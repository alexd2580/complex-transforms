.PHONY: build
build:
	gcc main.c -o complex_transform `sdl2-config --cflags --libs` -lm -lSDL2_image -lGL -lGLU -Wall

run: build
	./complex_transform

