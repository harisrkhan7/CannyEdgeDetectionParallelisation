#ifndef OUR_INDEX
#define OUR_INDEX
//inline int index(int width, int i, int j) {
//    return (i * width) + j;
//}
#define INDEX(width, i, j) (i * width) + j
#endif

__global__
void __alt_gaussian__(int height, int width, float* raw_buffer, float* filter_buffer) {
    float left [5][5] = {
        {2, 4, 5, 4, 2},
        {4, 9, 12, 9, 4},
        {5, 12, 15, 12, 5},
        {4, 9, 12, 9, 4},
        {2, 4, 5, 4, 2}
    };
    float right [5][5];
    float out [5][5];
    int thread_idx = threadIdx.x;
    int stride = blockDim.x;
    for (int i_top = 0; i_top < height; i_top++) {
        for (int j_left = thread_idx; j_left < width; j_left += stride) {

            for (int i = 0; i < 5; i++) {
                for (int j = 0; j < 5; j++) {
                    right[i][j] = raw_buffer[INDEX(width, i_top + i, j_left + j)];
                }
            }

            for (int i = 2; i < 3; i++) {
                for (int j = 2; j < 3; j++) {
                    out[i][j] = left[i][0] * right[0][j] + left[i][1] * right[1][j] + left[i][2] * right[2][j] + left[i][3] * right[3][j] + left[i][4] * right[4][j];
                }
            }

            filter_buffer[INDEX(width, i_top + 2, j_left + 2)] = out[2][2];
        }
    }
}

extern "C" alt_gaussio(int height, int width, float* raw_buffer, float* filter_buffer) {
    float* cuda_raw;
    float* cuda_filter;
    cudaMalloc(&cuda_raw, height * width * sizeof(float));
    cudaMalloc(&cuda_filter, height * width * sizeof(float));
    __alt_guassian__<<<width * height / 256, 256>>>(height, width, cuda_raw, cuda_filter);
    cudaMemcpy(filter_buffer, cuda_filter, height * width * sizeof(float), cudaMemcpyDeviceToHost);
}
