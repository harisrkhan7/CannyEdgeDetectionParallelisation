#include "upng/upng.h"
#include "omp.h"
#include <cmath>
#include <cstdio>
#include <algorithm>

upng_t* upng = NULL;
float* original_image_buffer;
float* gaussian_filter_buffer;
float* gradient_buffer;
float* suppression_buffer;
bool* strong_edge_buffer;
float* final_buffer;
float hysteresis_max = 0.0f;


float low_threshold = 0.12;
float high_threshold = 0.24;

enum Direction {
    EAST = 0,
    NORTHEAST = 45,
    NORTH = 90,
    NORTHWEST = 135
};

Direction* direction_buffer;

int width;
int height;

inline int index(int i, int j) {
    return (i * width) + j;
}

void matrix_multiply_3x3(float left[3][3], float right[3][3], float out[3][3]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            out[i][j] = left[i][0] * right[0][j] + left[i][1] * right[1][j] + left[i][2] * right[2][j];
        }
    }
}

void matrix_multiply_5x5(float left[5][5], float right[5][5], float out[5][5]) {
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            out[i][j] = left[i][0] * right[0][j] + left[i][1] * right[1][j] + left[i][2] * right[2][j] + left[i][3] * right[3][j] + left[i][4] * right[4][j];
        }
    }
}

// Load image into upng pointer
// Set width and height
void load_image() {
    printf("Loading image\n");
    upng = upng_new_from_file("img/lion.png");
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
    gaussian_filter_buffer = new float [width*height];
    gradient_buffer = new float [width * height];
    direction_buffer = new Direction [width * height];
    suppression_buffer = new float [width * height];
    strong_edge_buffer = new bool [width * height];
    final_buffer = new float [width * height];


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

void alt_gaussian() {
    // Precalculated Gaussian filter kernel
    float gaussian_filter_matrix[5][5] = {
        {2, 4, 5, 4, 2},
        {4, 9, 12, 9, 4},
        {5, 12, 15, 12, 5},
        {4, 9, 12, 9, 4},
        {2, 4, 5, 4, 2}
    };
    // Keep track of maximum brightness in output image so we can normalise
    float max_in_image = 0.0f;
    // Sliding window
    float image_piece [5][5];
    // Output of matrix multiplication
    float out [5][5];
    // For every top-left corner of a 5x5 sliding window over the input image
    #pragma omp parallel for shared(height, width, original_image_buffer, gaussian_filter_buffer, gaussian_filter_matrix) private(image_piece, out) reduction(max:max_in_image) num_threads(4)
    for (int i_top = 0; i_top < height - 4; i_top++) {
        for (int j_left = 0; j_left < width - 4; j_left++) {
            // Copy the 5x5 window into the image into a local matrix
            for (int i = 0; i < 5; i++) {
                for (int j = 0; j < 5; j++) {
                    image_piece[i][j] = original_image_buffer[index(i_top + i, j_left + j)];
                }
            }
            // Matrix multiply the gaussian kernel with the window of the image, into output matrix
            matrix_multiply_5x5(gaussian_filter_matrix, image_piece, out);
            // Divide by 159 and then copy the central pixel of the output matrix into the output buffer
            gaussian_filter_buffer[index(i_top + 2, j_left + 2)] = out[2][2] / 159.0f;
            // Update maxium brightness
            max_in_image = std::max(max_in_image, out[2][2]);
        }
    }
    // Normalise to avoid darkness
    // TODO: NORMALISE AT THE FINAL STEP ONLY
    max_in_image /= 159.0f;
    float scaling_factor = 1.0f;
    if (max_in_image != 0.0f) {
        scaling_factor = 255.0f / max_in_image;
    }
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            gaussian_filter_buffer[index(i, j)] *= scaling_factor;
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
    float max_in_image = 0.0f;
    float image_piece [3][3];
    float out_x [3][3];
    float out_y [3][3];
    #pragma omp parallel for shared(height, width, gaussian_filter_buffer, gradient_buffer, direction_buffer, sobel_convolve_x, sobel_convolve_y) private(image_piece, out_x, out_y) reduction(max:max_in_image) num_threads(4)
    for (int i_top = 0; i_top < height - 2; i_top++) {
        for (int j_left = 0; j_left < width - 2; j_left++) {

            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    image_piece[i][j] = gaussian_filter_buffer[index(i_top + i, j_left + j)];
                }
            }
            matrix_multiply_3x3(sobel_convolve_x, image_piece, out_x);
            matrix_multiply_3x3(sobel_convolve_y, image_piece, out_y);
            gradient_buffer[index(i_top + 1, j_left + 1)] = sqrt(out_x[1][1] * out_x[1][1] + out_y[1][1] * out_y[1][1]);
            max_in_image = std::max(max_in_image, gradient_buffer[index(i_top + 1, j_left + 1)]);
            double direction_d = 0.0;
            double gx = out_x[1][1];
            if (gx == 0.0) {
                direction_d = 90.0;
            } else if (gx == -0.0) {
                direction_d = -90.0;
            } else {
                direction_d = atan2(out_y[1][1], out_x[1][1]);
            }
            int direction_i = ((int) round(direction_d / 45.0)) * 45;
            Direction direction = EAST;
            if (direction_i == -90) {
                direction = NORTH;
            } else if (direction_i == -45) {
                direction = NORTHWEST;
            } else if (direction_i == 45) {
                direction = NORTHEAST;
            } else if (direction_i == 90) {
                direction = NORTH;
            }
            direction_buffer[index(i_top + 1, j_left + 1)] = direction;
            if (i_top == 4 && j_left == 4) {
                printf("OMP THREADS: %d\n", omp_get_num_threads());
            }
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

void suppress() {
    float gradient_piece [3][3];
    #pragma omp parallel for shared(height, width, gradient_buffer, direction_buffer, suppression_buffer) private(gradient_piece) reduction(max:hysteresis_max) num_threads(4)
    for (int i_top = 2; i_top < height - 2 - 2; i_top++) {
        for (int j_left = 2; j_left < width - 2 - 2; j_left++) {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    gradient_piece[i][j] = gradient_buffer[index(i_top + i, j_left + j)];
                }
            }
            float compare_one, compare_two;
            Direction direction = direction_buffer[index(i_top + 1, j_left + 1)];
            if (direction == EAST) {
                compare_one = gradient_piece[1][0];
                compare_two = gradient_piece[1][2];
            } else if (direction == NORTHEAST) {
                compare_one = gradient_piece[2][0];
                compare_two = gradient_piece[0][2];
            } else if (direction == NORTH) {
                compare_one = gradient_piece[0][1];
                compare_two = gradient_piece[2][1];
            } else if (direction == NORTHWEST) {
                compare_one = gradient_piece[0][0];
                compare_two = gradient_piece[2][2];
            }
            float current_gradient = gradient_piece[1][1];
            if (current_gradient >= compare_one && current_gradient >= compare_two) {
                suppression_buffer[index(i_top + 1, j_left + 1)] = current_gradient;
            } else {
                suppression_buffer[index(i_top + 1, j_left + 1)] = 0.0f;
            }

            hysteresis_max = std::max(hysteresis_max, current_gradient);
        }
    }
}

int safe_get_is_strong(int i, int j) {
    if (i < 0) {
        return false;
    }
    if (i >= height) {
        return false;
    }
    if (j < 0) {
        return false;
    }
    if (j >= width) {
        return false;
    }
    return strong_edge_buffer[index(i, j)];
}

bool has_strong_neighbour(int i, int j) {
    for (int q = -1; q < 2; q++) {
        for (int r = -1; r < 2; r++) {
            if (q == 0 && r == 0) {
                continue;
            }
            if (safe_get_is_strong(i + q, j + r)) {
                return true;
            }
        }
    }
    return false;
}

void hysteresis() {
    #pragma omp parallel for shared(height, width, strong_edge_buffer, suppression_buffer) num_threads(4)
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            float current_pixel = suppression_buffer[index(i, j)];
            if (current_pixel >= high_threshold * hysteresis_max) {
                strong_edge_buffer[index(i, j)] = true;
            } else {
                strong_edge_buffer[index(i, j)] = false;
            }
            if (current_pixel < low_threshold * hysteresis_max) {
                suppression_buffer[index(i, j)] = 0.0f;
            }
        }
    }
    #pragma omp parallel for shared(height, width, strong_edge_buffer, suppression_buffer, final_buffer) num_threads(4)
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            float current_pixel = suppression_buffer[index(i, j)];
            if (current_pixel > 0.0f) {
                if (strong_edge_buffer[index(i, j)]) {
                    final_buffer[index(i, j)] = 255.0f;
                } else {
                    if (has_strong_neighbour(i, j)) {
                        final_buffer[index(i, j)] = 255.0f;
                    } else {
                        final_buffer[index(i, j)] = 0.0f;
                    }
                }
            }

        }
    }
}

unsigned char* out_buffer;

void write_image(float* input_buffer) {
    printf("Writing image\n");
    out_buffer = new unsigned char [width * height];
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            // printf("%f ", input_buffer[index(i, j)]);
            out_buffer[index(i, j)] = (unsigned char) input_buffer[index(i, j)];
        }
        //printf("\n");
    }
    FILE* outfile = fopen("out/out.pgm", "wb");
    if (outfile == NULL) {
        fprintf(stderr, "Unable to open file for writing!\n");
    } else {
        fprintf(outfile, "P5\n%d\n%d\n255\n", width, height);
        fwrite(out_buffer, 1,  width * height, outfile);
    }
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
// test_gaussian_filter();
    load_image();
    convert_image();
    alt_gaussian();
    grad_dir();
    suppress();
    hysteresis();
    // process_image();
    write_image(final_buffer);
//    upng_free(upng);
    return 0;
}
