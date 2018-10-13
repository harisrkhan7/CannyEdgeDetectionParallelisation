#include "upng/upng.h"
#include <stdio.h>

int main(int argc, char** argv) {
	printf("Test!\n");
	upng_t* upng = NULL;
	upng = upng_new_from_file("img/bzaZ1.png");
	if (NULL != upng) {
		upng_decode(upng);
		if (UPNG_EOK == upng_get_error(upng)) {
			printf("Height: %d, Width: %d\n", upng_get_height(upng), upng_get_width(upng));
		}
		upng_free(upng);
	}
	return 0;
}
