#include <stdio.h>
#include "Vector3.h"

int main() {
    Vector3 v,w;
    v.x = 1.0;
    v.y = 2.0;
    v.z = 3.0;
    printf("hello world\n");
    printf("< %g %g %g >\n", v.x, v.y, v.z);
    w = v;
    printf("< %g %g %g >\n", w.x, w.y, w.z);
}
