  #include "stdio.h"
  #include "linearalgebraic.h"
int main()
{
    int i;
    static double a[16] =
        {0.2368, 0.2471, 0.2568, 1.2671,
         0.1968, 0.2071, 1.2168, 0.2271,
         0.1581, 1.1675, 0.1768, 0.1871,
         1.1161, 0.1254, 0.1397, 0.1490};
    static double b[4] = {1.8471, 1.7471, 1.6471, 1.5471};
    if (Gauss(a, b, 4) != 0)
        for (i = 0; i <= 3; i++)
            printf("x(%d)=%e\n", i, b[i]);
    printf("\n");
    static double a1[16] = {1.0, 3.0, 2.0, 13.0,
                            7.0, 2.0, 1.0, -2.0,
                            9.0, 15.0, 3.0, -2.0,
                            -2.0, -2.0, 11.0, 5.0};
    static double b1[8] = {9.0, 0.0, 6.0, 4.0,
                           11.0, 7.0, -2.0, -1.0};
    if (GaussJordan(a1, b1, 4, 2) != 0)
        for (i = 0; i <= 3; i++)
            printf("x(%d)=%13.7e,  %13.7e\n", i, b1[i * 2], b1[i * 2 + 1]);
    printf("\n");

    static double b2[13] = {13.0, 12.0, 11.0, 10.0, 9.0, 8.0, 7.0,
                            6.0, 5.0, 4.0, 3.0, 2.0, 1.0};
    static double d[5] = {3.0, 0.0, -2.0, 6.0, 8.0};
    if (TriMatrix(b2, 5, 13, d) > 0)
        for (i = 0; i <= 4; i++)
            printf("x(%d)=%13.7e\n", i, d[i]);

    return 1;
}