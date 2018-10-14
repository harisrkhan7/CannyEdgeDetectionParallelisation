#include "upng/upng.h"
#include <cmath>
#include <cstdio>

upng_t* upng = NULL;
float* original_image_buffer;
float* gaussian_filter_buffer;
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

float sobel_convolve_x [3][3] = {
    { 1.0f, 0.0f, -1.0f },
    { 2.0f, 0.0f, -2.0f },
    { 1.0f, 0.0f, -1.0f }
};

float sobel_convolve_y [3][3] = {
    { 1.0f, 2.0f, 1.0f },
    { 0.0f, 0.0f, 0.0f },
    { -1.0f, -2.0f, -1.0f }
};

// float gx_out [3][3] = {0.0f};
// float gy_out [3][3] = {0.0f};

void matrix_multiply_3x3(float left[3][3], float right[3][3], float out[3][3]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            out[i][j] = left[i][0] * right[0][j] + left[i][1] * right[1][j] + left[i][2] * right[2][j];
        }
    }
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
    gaussian_filter_buffer = new float [width*height];
    direction_buffer = new direction [width * height];

    printf("Converting image\n");
    const unsigned char* png_buffer = upng_get_buffer(upng);
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            original_image_buffer[index(i, j)] = (float) png_buffer[index(i, j)];
        }
    }
}
//Gaussian filter on 5 x 5 chunks of the image
void apply_gaussian_filter(){
    float gaussian_filter_matrix[5][5] ={
        {2,4,5,4,2},
        {4,9,12,9,4},
        {5,12,15,12,5},
        {4,9,12,9,4},
        {2,4,5,4,2}
    };
    for(int i=0;i<height;i++)
    {
        for(int j=0;j<width;j++)
        {
            gaussian_filter_buffer [index(i,j)] = 0;
            for(int k=0;k<5;k++)
            {
                int original_col_iterator = ((j/5)*5)+k;
                gaussian_filter_buffer[index(i,j)] += original_image_buffer[index(i,original_col_iterator)]*
                    gaussian_filter_matrix[k][j%5];
            }
            gaussian_filter_buffer[index(i,j)] /= 159;
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

void grad_dir() {
    float max_in_image = 0.0f;
    float image_piece [3][3];
    float out_x [3][3];
    float out_y [3][3];
    direction out_direction [3][3];
    for (int i_top = 0; i_top < height; i_top += 3) {
        for (int j_left = 0; j_left < width; j_left +=3) {

            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    image_piece[i][j] = original_image_buffer[index(i_top + i, j_left + j)];
                }
            }
            matrix_multiply_3x3(sobel_convolve_x, image_piece, out_x);
            matrix_multiply_3x3(sobel_convolve_y, image_piece, out_y);
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    gradient_buffer[index(i + i_top, j + j_left)] = sqrt(out_x[i][j] * out_x[i][j] + out_y[i][j] * out_y[i][j]);
                    if (gradient_buffer[index(i + i_top, j + j_left)] > max_in_image) {
                        max_in_image = gradient_buffer[index(i + i_top, j + j_left)];
                    }
                    double direction_d = 0.0;
                    double gx = out_x[i][j];
                    if (gx == 0.0) {
                        direction_d = 90.0;
                    } else if (gx == -0.0) {
                        direction_d = -90.0;
                    } else {
                        direction_d = atan(out_y[i][j] / out_x[i][j]);
                    }
                    int direction_rounded = ((int) round(direction_d / 45.0)) * 45;
                    direction rounded_direction = EAST;
                    if (direction_rounded == -90) {
                        rounded_direction = NORTH;
                    } else if (direction_rounded == -45) {
                        rounded_direction = NORTHWEST;
                    } else if (direction_rounded == 45) {
                        rounded_direction = NORTHEAST;
                    } else if (direction_rounded == 90) {
                        rounded_direction = NORTH;
                    }
                    direction_buffer[index(i + i_top, j + j_left)] = rounded_direction;
                }
            }

            //for (int i = i_top; i < i_top + 3; i++) {
            //    for (int j = j_left; j < j_left + 3; j++) {
            //        original_image_buffer[index(i, j)] = image_piece[i - i_top][j - j_left];
            //    }
            //}
        }
    }

    float scaling_factor = 1.0f;
    if (max_in_image != 0.0f) {
        scaling_factor = 255.0f / max_in_image;
    }

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            gradient_buffer[index(i, j)] *= scaling_factor;
        }
    }
}

unsigned char* out_buffer;

void write_image(float* input_buffer) {
    printf("Writing image\n");
    out_buffer = new unsigned char [width * height];
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            printf("%f ", input_buffer[index(i, j)]);
            out_buffer[index(i, j)] = (unsigned char) input_buffer[index(i, j)];
        }
        //printf("\n");
    }
    FILE* outfile = fopen("out/out.pgm", "wb");
    fprintf(outfile, "P5\n%d\n%d\n255\n", width, height);
    fwrite(out_buffer, 1, width * height, outfile);
}

void test_gaussian_filter(){
    load_image();
    convert_image();
    apply_gaussian_filter();
    write_image(gaussian_filter_buffer);
    for(int i=0;i<height;i++)
    {
        for(int j=0;j<width;j++)
        {
            printf("%f, ",gaussian_filter_buffer[index(i,j)]);
        }
    }
}
int main(int argc, char** argv) {
test_gaussian_filter();
//    load_image();
//    convert_image();
//    process_image();
//    write_image();
//    upng_free(upng);
    return 0;
}
