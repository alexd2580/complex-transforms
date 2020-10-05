/* gcc optical_illusion.c -o optical_illusion `sdl2-config --cflags --libs` -lm -lSDL2_image -lGL -lGLU -Wall  */

#include <GL/gl.h>
#include <GL/glu.h>
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include <SDL2/SDL_keycode.h>
#include <SDL2/SDL_opengl.h>
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define WINDOW_WIDTH 1000
#define WINDOW_HEIGHT 1000

#define PI 3.14159265358979323846

int32_t check_gl_error(char const* activity) {
    GLenum error = glGetError();
    if(error == GL_NO_ERROR) {
        return 0;
    }
    printf("Error when %s: %s\n", activity, gluErrorString(error));
    return 1;
}

SDL_Window* init_sdl_window(void) {
    return SDL_CreateWindow("Complex transformation", 0, 0, WINDOW_WIDTH, WINDOW_HEIGHT, SDL_WINDOW_OPENGL);
}

int32_t init_opengl_projection(void) {
    glViewport(0, 0, WINDOW_WIDTH, WINDOW_HEIGHT);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if(check_gl_error("initializing projection matrix")) {
        return 0;
    }

    // TODO dynamic aspect ration.
    glOrtho(-1.2, 1.2, -1.2, 1.2, -1.0, 1.0);

    // Initialize Modelview Matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    if(check_gl_error("initializing model matrix")) {
        return 0;
    }

    return 1;
}

typedef struct {
    double x;
    double y;
} vec2;

vec2 mk_vec2(double x, double y) { return (vec2){.x = x, .y = y}; }

typedef struct {
    double r;
    double t;
} vec2_polar;

vec2_polar mk_vec2_polar(double r, double t) { return (vec2_polar){.r = r, .t = t}; }

vec2 add(vec2 a, vec2 b) { return mk_vec2(a.x + b.x, a.y + b.y); }
vec2 mul(vec2 a, vec2 b) { return mk_vec2(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x); }
vec2 mul_scalar(vec2 a, double b) { return mk_vec2(a.x * b, a.y * b); }

vec2_polar to_polar(vec2 a) {
    double radius = sqrt(pow(a.x, 2) + pow(a.y, 2));
    double angle = acos(a.x / radius);
    angle = a.y > 0 ? angle : -angle;
    return mk_vec2_polar(radius, angle);
}

vec2 to_cartesian(vec2_polar a) { return mk_vec2(a.r * cos(a.t), a.r * sin(a.t)); }

vec2 exponentiate(vec2 a, double p) {
    vec2_polar polar = to_polar(a);
    vec2_polar exponentiated = mk_vec2_polar(pow(polar.r, p), p * polar.t);
    return to_cartesian(exponentiated);
}

void render_line(vec2 a, vec2 b) {
    glBegin(GL_LINES);
    glVertex2d(a.x, a.y);
    glVertex2d(b.x, b.y);
    glEnd();
}

void render_point(vec2 a) {
    glBegin(GL_LINE_STRIP);
    for(int i = 0; i < 21; i++) {
        glVertex2d(a.x + 0.02 * cos(i * 2 * PI / 20), a.y + 0.02 * sin(i * 2 * PI / 20));
    }
    glEnd();
}

void render_coordinate_system(void) {
    render_line(mk_vec2(-1, 0), mk_vec2(1, 0));
    render_line(mk_vec2(0, -1), mk_vec2(0, 1));
}

vec2 jump(double f) { return f - floor(f) > 0.5 ? mk_vec2(1, 1) : mk_vec2(-1, -1); }

int map_range(double* x, double a1, double a2, double b1, double b2) {
    if(*x < a1 || *x > a2) {
        return 0;
    }

    // Shorten the interval to exclude the order values (often the computation is not defined for those).
    int positive_interval = b2 > b1;
    b1 += positive_interval ? 0.00000001 : -0.00000001;
    b2 -= positive_interval ? 0.00000001 : -0.00000001;

    *x = ((*x - a1) / (a2 - a1)) * (b2 - b1) + b1;
    return 1;
}

vec2 upper_wings_mask(double x) {
    double y = 1.5 * sqrt((-fabs(fabs(x) - 1)) * fabs(3 - fabs(x)) / ((fabs(x) - 1) * (3 - fabs(x)))) *
                   (1 + fabs(fabs(x) - 3) / (fabs(x) - 3)) * sqrt(1 - pow(x / 7, 2)) +
               (4.5 + 0.75 * (fabs(x - 0.5) + fabs(x + 0.5)) - 2.75 * (fabs(x - 0.75) + fabs(x + 0.75))) *
                   (1 + fabs(1 - fabs(x)) / (1 - fabs(x)));
    return mk_vec2(x, y);
}

vec2 shoulders(double x) {
    double y = (2.71052 + 1.5 - 0.5 * fabs(x) - 1.35526 * sqrt(4 - pow(fabs(x) - 1, 2))) *
               sqrt(fabs(fabs(x) - 1) / (fabs(x) - 1));
    return mk_vec2(x, y);
}

vec2 lower_wings(double x) {
    double y = (-3) * sqrt(1 - pow(x / 7, 2)) * sqrt(fabs(fabs(x) - 4) / (fabs(x) - 4));
    return mk_vec2(x, y);
}

vec2 lower_flaps(double x) {
    double y = fabs(x / 2) - 0.0913722 * pow(x, 2) - 3 + sqrt(1 - pow(fabs(fabs(x) - 2) - 1, 2));
    return mk_vec2(x, y);
}

/* vec2 (double x) { */
/*     if(x < 0.5) { */
/*         x = 2.0 * x * (7 - 4.00001) + 4.00001; */
/*     } else { */
/*         x = 2.0 * (x - 0.5) * (-4.00001 + 7) - 7; */
/*     } */
/*     double y = (-3) * sqrt(1 - pow(x / 7, 2)) * sqrt(fabs(fabs(x) - 4) / (fabs(x) - 4)); */
/*     return mk_vec2(x, y); */
/* } */
/*  */
/* vec2 pt3(double x) { */
/*     x = x * 8 - 4; */
/*     double y = fabs(x / 2) - 0.0913722 * pow(x, 2) - 3 + sqrt(1 - pow(fabs(fabs(x) - 2) - 1, 2)); */
/*     return mk_vec2(x, y); */
/* } */
/*  */
/* vec2 pt4(double x) { */
/*     if(x < 0.5) { */
/*         x = 2.0 * x * (2.99999 - 1.00001) + 1.00001; */
/*     } else { */
/*         x = 2.0 * (x - 0.5) * (-1.00001 + 2.99999) - 2.99999; */
/*     } */
/*     double y = (2.71052 + 1.5 - 0.5 * fabs(x) - 1.35526 * sqrt(4 - pow(fabs(x) - 1, 2))) * */
/*                sqrt(fabs(fabs(x) - 1) / (fabs(x) - 1)); */
/*     return mk_vec2(x, y); */
/* } */

vec2 batman(double x) {
    if(map_range(&x, 0, 0.125, -1, 1))
        return upper_wings_mask(x);
    if(map_range(&x, 0.125, 0.25, 1, 3))
        return shoulders(x);
    if(map_range(&x, 0.25, 0.375, 3, 7))
        return upper_wings_mask(x);
    if(map_range(&x, 0.375, 0.5, 7, 4))
        return lower_wings(x);
    if(map_range(&x, 0.5, 0.625, 4, -4))
        return lower_flaps(x);
    if(map_range(&x, 0.625, 0.75, -4, -7))
        return lower_wings(x);
    if(map_range(&x, 0.75, 0.875, -7, -3))
        return upper_wings_mask(x);
    if(map_range(&x, 0.875, 1, -3, -1))
        return shoulders(x);

    return mk_vec2(0, 0);
}

vec2 test(double x) {
    vec2 v = batman(x);
    // TODO assert not nan
    /* printf("%lf => %lf %lf\n", x, v.x, v.y); */
    return mul_scalar(v, 1.0 / 7);
    return jump(x);
}

vec2 compute_coefficient(int index, int sample_rate) {
    // We sample the interval [0, 1), meaning that
    // the sample rate equals the number of samples.
    double step_size = 1.0 / sample_rate;
    vec2 sum = mk_vec2(0, 0);
    for(int sample_index = 0; sample_index < sample_rate; sample_index++) {
        double parameter = sample_index * step_size;
        vec2 function_value = test(parameter);
        double angle = -index * 2 * PI * parameter;
        vec2 inverse_coefficient = to_cartesian(mk_vec2_polar(1, angle));
        sum = add(sum, mul(function_value, inverse_coefficient));
    }
    return mul_scalar(sum, step_size);
}

vec2* compute_coefficients(int num_coefficients, int sample_rate) {
    vec2* coefficients = (vec2*)malloc(num_coefficients * sizeof(vec2));
    int index = 0;
    for(int i = 0; i < num_coefficients; i++) {
        printf("%d\n", i);
        index += i % 2 == 0 ? -i : i;
        coefficients[i] = compute_coefficient(index, sample_rate);
    }
    return coefficients;
}

int zoom = 0;
int max_render_coeffifient_rank = 100000;

vec2 render_coefficients(vec2* coefficients, int num_coefficients, double t) {
    vec2 start = mk_vec2(0, 0);
    int index = 0;
    for(int i = 0; i < num_coefficients && i < max_render_coeffifient_rank; i++) {
        index += i % 2 == 0 ? -i : i;
        double angle = index * 2 * PI * t;
        vec2 base_arrow = to_cartesian(mk_vec2_polar(1, angle));
        vec2 arrow = mul(coefficients[i], base_arrow);
        vec2 end = add(start, arrow);
        render_line(start, end);
        start = end;
    }
    return start;
}

void free_coefficients(vec2* coefficients) { free(coefficients); }

GLuint create_empty_texture(void) {
    GLuint texture;
    glGenTextures(1, &texture);

    glBindTexture(GL_TEXTURE_2D, texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    unsigned char* data = (unsigned char*)malloc(3 * 1000 * 1000 * sizeof(unsigned char));
    memset(data, 0, 1000 * 1000 * 3);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1000, 1000, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
    glBindTexture(GL_TEXTURE_2D, 0);

    check_gl_error("creating empty texture");
    return texture;
}

int clamp(int min, int x, int max) { return x < min ? min : x > max ? max : x; }

void draw_to_texture_at(GLuint texture, vec2 pos) {
    glBindTexture(GL_TEXTURE_2D, texture);
    unsigned char five_pixels[5 * 3];
    for(int i = 0; i < 5; i++) {
        five_pixels[3 * i + 0] = 255;
        five_pixels[3 * i + 1] = 0;
        five_pixels[3 * i + 2] = 0;
    }
    int pos_x = clamp(0, (int)(1000 * (pos.x + 1) / 2) - 2, 995);
    int pos_y = clamp(0, (int)(1000 * (pos.y + 1) / 2) - 2, 995);
    glTexSubImage2D(GL_TEXTURE_2D, 0, pos_x, pos_y + 0, 5, 1, GL_RGB, GL_UNSIGNED_BYTE, five_pixels);
    glTexSubImage2D(GL_TEXTURE_2D, 0, pos_x, pos_y + 1, 5, 1, GL_RGB, GL_UNSIGNED_BYTE, five_pixels);
    glTexSubImage2D(GL_TEXTURE_2D, 0, pos_x, pos_y + 2, 5, 1, GL_RGB, GL_UNSIGNED_BYTE, five_pixels);
    glTexSubImage2D(GL_TEXTURE_2D, 0, pos_x, pos_y + 3, 5, 1, GL_RGB, GL_UNSIGNED_BYTE, five_pixels);
    glTexSubImage2D(GL_TEXTURE_2D, 0, pos_x, pos_y + 4, 5, 1, GL_RGB, GL_UNSIGNED_BYTE, five_pixels);
    glBindTexture(GL_TEXTURE_2D, 0);

    check_gl_error("drawing to texture");
}

void free_texture(GLuint texture) { glDeleteTextures(1, &texture); }

int main(int argc, char* argv[]) {
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 2);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);
    SDL_Window* window = init_sdl_window();
    SDL_GLContext context = SDL_GL_CreateContext(window);
    init_opengl_projection();
    glClearColor(0.f, 0.f, 0.f, 0.f);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_TEXTURE_2D);

    check_gl_error("initializing gl");

    GLuint tex = create_empty_texture();

    int num_coefficients = 10000;
    max_render_coeffifient_rank = num_coefficients;
    int sample_rate = 20000;
    vec2* coefficients = compute_coefficients(num_coefficients, sample_rate);
    vec2 last_pos = mk_vec2(0, 0);

    double t = 0.0;
    double step_size = 0.0001;
    int32_t running = 1;
    while(running) {
        // Handle events.
        SDL_Event event;
        while(SDL_PollEvent(&event)) {
            if(event.type == SDL_KEYDOWN) {
                switch(event.key.keysym.sym) {
                case SDLK_ESCAPE:
                    running = 0;
                    break;
                case SDLK_KP_MULTIPLY:
                    max_render_coeffifient_rank *= 2;
                    printf("max_render_coeffifient_rank %d\n", max_render_coeffifient_rank);
                    break;
                case SDLK_KP_DIVIDE:
                    max_render_coeffifient_rank /= 2;
                    printf("max_render_coeffifient_rank %d\n", max_render_coeffifient_rank);
                    break;
                case SDLK_KP_PLUS:
                    max_render_coeffifient_rank += 1;
                    printf("max_render_coeffifient_rank %d\n", max_render_coeffifient_rank);
                    break;
                case SDLK_KP_MINUS:
                    max_render_coeffifient_rank -= 1;
                    printf("max_render_coeffifient_rank %d\n", max_render_coeffifient_rank);
                    break;
                case SDLK_z:
                    zoom = 1 - zoom;
                    break;
                case SDLK_f:
                    step_size *= 2;
                    break;
                case SDLK_s:
                    step_size /= 2;
                    break;
                default:
                    break;
                }
            } else if(event.type == SDL_QUIT) {
                running = 0;
            }
        }

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        if (zoom) {
            glMatrixMode(GL_PROJECTION);
            glPushMatrix();
            glScaled(10.0, 10.0, 0.0);
            glTranslated(-last_pos.x, -last_pos.y, 0.0);
        }

        glColor3d(1.0, 1.0, 1.0);
        glBindTexture(GL_TEXTURE_2D, tex);
        glBegin(GL_QUADS);
        glTexCoord2i(0, 0);
        glVertex2i(-1, -1);
        glTexCoord2i(1, 0);
        glVertex2i(1, -1);
        glTexCoord2i(1, 1);
        glVertex2i(1, 1);
        glTexCoord2i(0, 1);
        glVertex2i(-1, 1);
        glEnd();
        glBindTexture(GL_TEXTURE_2D, 0);

        check_gl_error("rendering texture");

        glColor3d(0.7, 0.7, 0.7);
        render_coordinate_system();

        glColor3d(0.2, 0.2, 0.5);
        vec2 final_point = render_coefficients(coefficients, num_coefficients, t);
        last_pos = final_point;

        glColor3d(1.0, 0.0, 0.0);
        render_point(final_point);
        draw_to_texture_at(tex, final_point);

        SDL_GL_SwapWindow(window);
        SDL_Delay(13);

        if (zoom) {
            glMatrixMode(GL_PROJECTION);
            glPopMatrix();
        }

        t += step_size;
        if(t >= 1.0) {
            t -= 1.0;
        }
    }

    free_texture(tex);
    free_coefficients(coefficients);

    return 0;
}
