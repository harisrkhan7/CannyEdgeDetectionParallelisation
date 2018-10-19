#include "upng/upng.h"
#include "omp.h"
#include <cmath>
#include <cstdio>
#include <mpi.h>
#include <unistd.h>
#include <algorithm>

enum Direction {
    EAST = 0,
    NORTHEAST = 45,
    NORTH = 90,
    NORTHWEST = 135
};

upng_t* upng = NULL;
float* original_image_buffer;
float* gaussian_filter_buffer;
float* gradient_buffer;
float* suppression_buffer;
bool* strong_edge_buffer;
float* final_buffer;
float hysteresis_max = 0.0f;
Direction* direction_buffer;

char* input_filename;
int width;
int height;
float low_threshold = 0.12;
float high_threshold = 0.24;

int global_height;
int global_width;

#ifndef OUR_INDEX
#define OUR_INDEX
inline int index(int i, int j) {
    return (i * width) + j;
}
#endif

void matrix_multiply_3x3(float left[3][3], float right[3][3], float out[3][3]) {
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            out[i][j] = left[i][0] * right[0][j] + left[i][1] * right[1][j] + left[i][2] * right[2][j];
        }
    }
}

void matrix_multiply_5x5(float left[5][5], float right[5][5], float out[5][5]) {
    for (int i = 2; i < 3; i++) {
        for (int j = 2; j < 3; j++) {
            out[i][j] = left[i][0] * right[0][j] + left[i][1] * right[1][j] + left[i][2] * right[2][j] + left[i][3] * right[3][j] + left[i][4] * right[4][j];
        }
    }
}

// Load image into upng pointer
// Set width and height
void load_image() {
    fprintf(stderr, "Loading image\n");
    upng = upng_new_from_file(input_filename);
    if (NULL != upng) {
        upng_decode(upng);
        if (UPNG_EOK == upng_get_error(upng)) {
            width = upng_get_width(upng);
            height = upng_get_height(upng);
            fprintf(stderr, "Height: %d, Width: %d\n", height, width);
            fprintf(stderr, "Bits per pixel: %d, Format: %d\n", upng_get_bpp(upng), upng_get_format(upng));
        }
    }
}

// Convert the upng image into an array of floats;
// Create gradient and direction buffers

void prepare_memory() {
    original_image_buffer = new float [width * (height + 6)];
    gaussian_filter_buffer = new float [width * (height + 6)];
    gradient_buffer = new float [width * (height + 6)];
    direction_buffer = new Direction [width * (height + 6)];
    suppression_buffer = new float [width * (height + 6)];
    strong_edge_buffer = new bool [width * (height + 6)];
    final_buffer = new float [width * (height + 6)];
}

void convert_image() {
    


    fprintf(stderr, "Converting image\n");
    const unsigned char* png_buffer = upng_get_buffer(upng);
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            original_image_buffer[index(i, j)] = (float) png_buffer[index(i, j)];
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
    #pragma omp parallel for shared(height, width, original_image_buffer, gaussian_filter_buffer, gaussian_filter_matrix) private(image_piece, out) reduction(max:max_in_image) schedule(dynamic)
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
    #pragma omp parallel for shared(height, width, gaussian_filter_buffer, scaling_factor) schedule(dynamic)
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            gaussian_filter_buffer[index(i, j)] *= scaling_factor;
        }
    }
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
    #pragma omp parallel for shared(height, width, gaussian_filter_buffer, gradient_buffer, direction_buffer, sobel_convolve_x, sobel_convolve_y) private(image_piece, out_x, out_y) reduction(max:max_in_image) schedule(dynamic)
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
                fprintf(stderr, "OMP threads: %d\n", omp_get_num_threads());
            }
        }
    }

    float scaling_factor = 1.0f;
    if (max_in_image != 0.0f) {
        scaling_factor = 255.0f / max_in_image;
    }
    #pragma omp parallel for shared(height, width, gaussian_filter_buffer, scaling_factor) schedule(dynamic)
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            gradient_buffer[index(i, j)] *= scaling_factor;
        }
    }
}

void suppress() {
    float gradient_piece [3][3];
    #pragma omp parallel for shared(height, width, gradient_buffer, direction_buffer, suppression_buffer) private(gradient_piece) reduction(max:hysteresis_max) schedule(dynamic)
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
    #pragma omp parallel for shared(height, width, strong_edge_buffer, suppression_buffer) schedule(dynamic) 
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
    #pragma omp parallel for shared(height, width, strong_edge_buffer, suppression_buffer, final_buffer) schedule(dynamic)
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

void write_image(float* input_buffer, unsigned int leave_out, bool append) {
    printf("Writing image\n");
    out_buffer = new unsigned char [width * height];
    #pragma omp parallel for shared(height, width, out_buffer, input_buffer) schedule(dynamic)
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            // fprintf(stderr, "%f ", input_buffer[index(i, j)]);
            out_buffer[index(i, j)] = (unsigned char) input_buffer[index(i, j)];
        }
        //fprintf(stderr, "\n");
    }
    
    if (append == false) {
        FILE* outfile = fopen("out/out.pgm", "wb");
        fprintf(outfile, "P5\n%d\n%d\n255\n", global_width, global_height);
        fwrite(out_buffer, 1, width * (height - leave_out), outfile);
        fclose(outfile);
    } else {
        FILE* outfile = fopen("out/out.pgm", "ab");
        fwrite(out_buffer + (leave_out * width), 1, width * (height -  leave_out), outfile);
        fclose(outfile);
    }
    
}

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    input_filename = argv[1];
    
    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    //For termination
    int commSize = world_size;

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    MPI_Status status;
    
    //Extra rows for neighbourhood operation
    int num_extra_rows = 3;
    
    //Number of rows per process
    int num_rows_per_process;
    
    //Total Items to Recv 
    int total_items_to_recv;
    
    // If the world size is 1, just do everything in this process and then exit
    // This MUST be before the world size stuff below
    if (world_size == 1) {
        load_image();
        prepare_memory();
        convert_image();
        global_height = height;
        global_width = width;
        alt_gaussian();
        grad_dir();
        suppress();
        hysteresis();
        if (world_rank == 0) {
            write_image(final_buffer, 0, false);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        return 0;
    }
    
    if (world_rank == 0) {
        load_image();
        prepare_memory();
        convert_image();
        global_height = height;
        global_width = width;
        
        int resolution[2] = { width, height };
        
        MPI_Bcast(&resolution, 2, MPI_INT, world_rank, MPI_COMM_WORLD);
        
        num_rows_per_process = height/world_size;
        
        //Divide and send data to every processor 
        int num_rows_to_send = num_rows_per_process + (2 * num_extra_rows);
        
        int total_items_to_send = num_rows_to_send * width;
        
        int original_buffer_index = (num_rows_per_process * width) - (num_extra_rows * width);
        
        
        float* sending_buffer;
        
        int total_processors;
        
        if ((height % world_size) == 0) {
            total_processors = world_size;            
        } else {
            total_processors = world_size - 1;
        }

        //Start from rank 1 as rank 0 is master process itself
        //Leave the last process out in case not completely divisible
        for (int destination = 1; destination < total_processors; destination++) {
            if (destination == (world_size -1)) {
                total_items_to_send -= (2 * num_extra_rows * width);
            }
            sending_buffer = new float[total_items_to_send];
            for (int i = 0; i < total_items_to_send; i++, original_buffer_index++) {
                sending_buffer[i] = original_image_buffer[original_buffer_index];
            }
            
            MPI_Send(sending_buffer, total_items_to_send, MPI_FLOAT, destination, 0, MPI_COMM_WORLD);
            if (destination != (world_size -1)) {
                original_buffer_index -= 2 * num_extra_rows * width;     
            }
        }
        
        //Send last processor more work if it can't be evenly divided
        if ((height % world_size) != 0) {
            //Give the remaining work to the last processor
            int last_processor_rank = world_rank - 1;
        
            int num_rows_for_last_processor = height - (num_rows_per_process * (world_size -  1)) + num_extra_rows;
            total_items_to_send = num_rows_for_last_processor * width;
            sending_buffer = new float[total_items_to_send];
            
            for (int i = 0; i < total_items_to_send; i++, original_buffer_index++) {
                sending_buffer[i] = original_image_buffer[original_buffer_index];
            }
            
            int destination = world_size - 1;
            MPI_Send(sending_buffer, total_items_to_send, MPI_FLOAT, destination, 0, MPI_COMM_WORLD);
        }
        
//        height = num_rows_per_process;
        //Compute Master Chunk 
        alt_gaussian();
        grad_dir();
        suppress();
        hysteresis();
        //Write Results and Send Ack to Next Process
        
        write_image(final_buffer, num_extra_rows, false);
        
        
        MPI_Barrier(MPI_COMM_WORLD);
        int ack = 0;
        MPI_Send(&ack, 1, MPI_INT, 1, 1, MPI_COMM_WORLD);
        
        //Wait for all processes to finish and terminate
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
    } else if (world_rank == world_size - 1) {
        int resolution[2];
        MPI_Bcast(&resolution, 2, MPI_INT, 0, MPI_COMM_WORLD);
        width = resolution[0];
        height = resolution[1];
        num_rows_per_process = height/world_size;
        
        if ((height % world_size) != 0) {
            //Give the remaining work to the last processor
            int num_rows_for_last_processor = height - (num_rows_per_process * (world_size -  1));
            num_rows_for_last_processor += num_extra_rows;
            
            
            total_items_to_recv = num_rows_for_last_processor * width;
            
            height = num_rows_for_last_processor;
        } else {
            total_items_to_recv = ((height / world_size) + num_extra_rows) * width;
            
            height = num_rows_per_process + num_extra_rows;
        }
        
        prepare_memory();
        //original_image_buffer = new float[total_items_to_recv];
        MPI_Recv(original_image_buffer, total_items_to_recv, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);
        //Compute Image 
        
        alt_gaussian();
        grad_dir();
        suppress();
        hysteresis();
        
        //Wait for Ack from Previous Process
        MPI_Barrier(MPI_COMM_WORLD);
        int ack=0;
        MPI_Recv(&ack, 1, MPI_INT, (world_rank - 1), 1, MPI_COMM_WORLD, &status);
        
        //Write Results 
        
        write_image(final_buffer, num_extra_rows, true);
        
        //Wait for all processes to finish and terminate
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
    } else {
        int resolution[2];
        MPI_Bcast(&resolution, 2, MPI_INT ,0, MPI_COMM_WORLD);
        width = resolution[0];
        height = resolution[1];
        num_rows_per_process = height / world_size;
        
        total_items_to_recv = (height / world_size) * width;
        total_items_to_recv += (2 * num_extra_rows * width);
        height = num_rows_per_process;
        prepare_memory();
        //original_image_buffer = new float[total_items_to_recv];
        MPI_Recv(original_image_buffer, total_items_to_recv, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);
        
        //Compute Image
        alt_gaussian();
        grad_dir();
        suppress();
        hysteresis();
        
        
        //Wait for Ack from Previous Process
        MPI_Barrier(MPI_COMM_WORLD);
        int ack = 0;
        MPI_Recv(&ack, 1, MPI_INT, (world_rank - 1), 1, MPI_COMM_WORLD, &status);
        
        //Write the file
         write_image(final_buffer,num_extra_rows,true);
        //Send Ack to next process
        ack = 1;
        MPI_Send(&ack, 1, MPI_INT, (world_rank + 1), 1, MPI_COMM_WORLD);
        
        //Wait for all processes to finish and terminate
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
    }
    
    
    return 0;
}
