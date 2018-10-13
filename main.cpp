#include "upng/upng.h"
#include <cstdio>

upng_t* upng = NULL;
float* original_image_buffer;
float* gradient_buffer;

enum direction {
    EAST = 0,
    NORTHEAST = 45,
    NORTH = 90,
    NORTHWEST = 135
};

direction* direction_buffer;

int width;
int height;

inline int index(int i, int j) {
    return (i * width) + j;
}


// Load image into upng pointer
// Set width and height
void load_image() {
    printf("Loading image\n");
    upng = upng_new_from_file("img/gazelle.png");
    if (NULL != upng) {
        upng_decode(upng);
        if (UPNG_EOK == upng_get_error(upng)) {
            width = upng_get_width(upng);
            height = upng_get_height(upng);
            printf("Height: %d, Width: %d\n", height, width);
            printf("Bits per pixel: %d, Format: %d\n", upng_get_bpp(upng), upng_get_format(upng));
        }
    }
}

// Convert the upng image into an array of floats;
// Create gradient and direction buffers
void convert_image() {
    original_image_buffer = new float [width * height];
    gradient_buffer = new float [width * height];
    direction_buffer = new direction [width * height];

    printf("Converting image\n");
    const unsigned char* png_buffer = upng_get_buffer(upng);
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            original_image_buffer[index(i, j)] = (float) png_buffer[index(i, j)];
        }
    }
}

void process_image() {
    printf("Processing image\n");
    // Pad to nearest multiple of 15 AND CONVERT TO FLOAT ARRAY
    // Create array #2 (floats) for gradient and array #3 (chars) for direction
    // Step 1: 5 * 5 guassian blurring
    // Step 2: Get gradient and direction arrays
    // Step 3: Do the non-maximum suppression
    // Step 4: Thresholding (pixel by pixel)
    // Step 5: Write to file
}

void write_image() {
    printf("Writing image\n");
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            printf("%f ", original_image_buffer[index(i, j)]);
        }
        printf("\n");
    }
}

int main(int argc, char** argv) {
    load_image();
    convert_image();
    process_image();
    write_image();
    upng_free(upng);
    return 0;
}
