#include <stdlib.h>
#include <stdio.h>

int main() {
    float *heap_array = (float *)malloc(sizeof(float) * 1e10);

    if(heap_array == NULL) {
	    perror("malloc failed");
	    exit(1);
    }
    printf("%i\n", sizeof(float));
    printf("%i\n", sizeof(double));

    long int m = 1000005;
    heap_array[m] = m;
    printf("%f\n", heap_array[m]);

    free(heap_array);
    return 0;
}
