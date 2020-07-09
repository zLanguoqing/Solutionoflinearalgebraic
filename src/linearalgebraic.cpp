#include "linearalgebraic.h"
/*
a 系数矩阵，返回时被破坏
b 右边的常量向量，返回解向量
n 阶数
   */
int Gauss(double a[], double b[], int n)
{
    int *js, l, k, i, j, is, p, q;
    double d, t;
    js = (int *)malloc(n * sizeof(int));
    l = 1;
    for (k = 0; k <= n - 2; k++)
    {
        d = 0.0;
        for (i = k; i <= n - 1; i++)
            for (j = k; j <= n - 1; j++)
            {
                t = fabs(a[i * n + j]);
                if (t > d)
                {
                    d = t;
                    js[k] = j;
                    is = i;
                }
            }
        if (d < MIN)
            l = 0;
        else
        {
            if (js[k] != k)
                for (i = 0; i <= n - 1; i++)
                {
                    p = i * n + k;
                    q = i * n + js[k];
                    t = a[p];
                    a[p] = a[q];
                    a[q] = t;
                }
            if (is != k)
            {
                for (j = k; j <= n - 1; j++)
                {
                    p = k * n + j;
                    q = is * n + j;
                    t = a[p];
                    a[p] = a[q];
                    a[q] = t;
                }
                t = b[k];
                b[k] = b[is];
                b[is] = t;
            }
        }
        if (l == 0)
        {
            free(js);
            printf("fail\n");
            return 0;
        }
        d = a[k * n + k];
        for (j = k + 1; j <= n - 1; j++)
        {
            p = k * n + j;
            a[p] = a[p] / d;
        }
        b[k] = b[k] / d;
        for (i = k + 1; i <= n - 1; i++)
        {
            for (j = k + 1; j <= n - 1; j++)
            {
                p = i * n + j;
                a[p] = a[p] - a[i * n + k] * a[k * n + j];
            }
            b[i] = b[i] - a[i * n + k] * b[k];
        }
    }
    d = a[(n - 1) * n + n - 1];
    if (fabs(d) < MIN)
    {
        free(js);
        printf("fail\n");
        return 0;
    }
    b[n - 1] = b[n - 1] / d;
    for (i = n - 2; i >= 0; i--)
    {
        t = 0.0;
        for (j = i + 1; j <= n - 1; j++)
            t = t + a[i * n + j] * b[j];
        b[i] = b[i] - t;
    }
    js[n - 1] = n - 1;
    for (k = n - 1; k >= 0; k--)
        if (js[k] != k)
        {
            t = b[k];
            b[k] = b[js[k]];
            b[js[k]] = t;
        }
    free(js);
    return 1;
}
/*  
a 稀疏矩阵
b 右端的m组常量向量
n 方程组的阶数
m 右端向量的组数
 */
int GaussJordan(double a[], double b[],int n,int m)
{
    int *js, l, k, i, j, is, p, q;
    double d, t;
    js =(int *) malloc(n * sizeof(int));
    l = 1;
    for (k = 0; k <= n - 1; k++)
    {
        d = 0.0;
        for (i = k; i <= n - 1; i++)
            for (j = k; j <= n - 1; j++)
            {
                t = fabs(a[i * n + j]);
                if (t > d)
                {
                    d = t;
                    js[k] = j;
                    is = i;
                }
            }
        if (d + 1.0 == 1.0)
            l = 0;
        else
        {
            if (js[k] != k)
                for (i = 0; i <= n - 1; i++)
                {
                    p = i * n + k;
                    q = i * n + js[k];
                    t = a[p];
                    a[p] = a[q];
                    a[q] = t;
                }
            if (is != k)
            {
                for (j = k; j <= n - 1; j++)
                {
                    p = k * n + j;
                    q = is * n + j;
                    t = a[p];
                    a[p] = a[q];
                    a[q] = t;
                }
                for (j = 0; j <= m - 1; j++)
                {
                    p = k * m + j;
                    q = is * m + j;
                    t = b[p];
                    b[p] = b[q];
                    b[q] = t;
                }
            }
        }
        if (l == 0)
        {
            free(js);
            printf("fail\n");
            return 0;
        }
        d = a[k * n + k];
        for (j = k + 1; j <= n - 1; j++)
        {
            p = k * n + j;
            a[p] = a[p] / d;
        }
        for (j = 0; j <= m - 1; j++)
        {
            p = k * m + j;
            b[p] = b[p] / d;
        }
        for (j = k + 1; j <= n - 1; j++)
            for (i = 0; i <= n - 1; i++)
            {
                p = i * n + j;
                if (i != k)
                    a[p] = a[p] - a[i * n + k] * a[k * n + j];
            }
        for (j = 0; j <= m - 1; j++)
            for (i = 0; i <= n - 1; i++)
            {
                p = i * m + j;
                if (i != k)
                    b[p] = b[p] - a[i * n + k] * b[k * m + j];
            }
    }
    for (k = n - 1; k >= 0; k--)
        if (js[k] != k)
            for (j = 0; j <= m - 1; j++)
            {
                p = k * m + j;
                q = js[k] * m + j;
                t = b[p];
                b[p] = b[q];
                b[q] = t;
            }
    free(js);
    return 1;
}
/*
三对角线方程组的追赶法求解
b 放置三对角线的矩阵
n 矩阵的阶数
m 三对角线上的元素个数 m = 3n-2
d 右端常数向量，返回是方程组的解
返回值 0 工作失败， -1 m的值不对， 1 正常返回  
*/
int TriMatrix(double b[], int n, int m, double d[])
{
    int k, j;
    double s;
    if (m != (3 * n - 2))
    {
        printf("err\n");
        return -1;
    }
    for (k = 0; k <= n - 2; k++)
    {
        j = 3 * k;
        s = b[j];
        if (fabs(s) + 1.0 == 1.0)
        {
            printf("fail\n");
            return 0;
        }
        b[j + 1] = b[j + 1] / s;
        d[k] = d[k] / s;
        b[j + 3] = b[j + 3] - b[j + 2] * b[j + 1];
        d[k + 1] = d[k + 1] - b[j + 2] * d[k];
    }
    s = b[3 * n - 3];
    if (fabs(s) + 1.0 == 1.0)
    {
        printf("fail\n");
        return 0;
    }
    d[n - 1] = d[n - 1] / s;
    for (k = n - 2; k >= 0; k--)
        d[k] = d[k] - b[3 * k + 1] * d[k + 1];
    return 1;
}
