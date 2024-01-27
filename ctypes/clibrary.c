#include <stdio.h>
#include <math.h>
#include "middle.c"

char* display(char* str, int age){
    printf("my name is %s, and my age is %d\n", str, age);
    middlefun();
    return "completed";
}

float styblinski_tang(float x, float y){
    return 0.5 * (pow(x, 4) - 16*pow(x, 2) + 5*x + pow(y, 4) - 16*pow(y, 2) + 5*y);
}