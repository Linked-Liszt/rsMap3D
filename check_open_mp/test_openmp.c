#include <stdio.h>
#include <omp.h>

int main(int argc, char** argv){

    printf("Max Threads %d\n", omp_get_max_threads());

    return 0;
}