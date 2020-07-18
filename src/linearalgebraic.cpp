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
int GaussJordan(double a[], double b[], int n, int m)
{
    int *js, l, k, i, j, is, p, q;
    double d, t;
    js = (int *)malloc(n * sizeof(int));
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
int ChaseTridiagonal(double b[], int n, int m, double d[])
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

/*
一般带型方程组求解
b 存放带区的元素
d 右端m组常数向量
n 方程组阶数
l 矩阵的半带宽
il 系数矩阵的带宽 il = 2l+1
m 方程组右端常量向量组数 
*/
int Band(double b[], double d[], int n, int l, int il, int m)
{
    int ls, k, i, j, is, u, v;
    double p, t;
    if (il != (2 * l + 1))
    {
        printf("fail\n");
        return -1;
    }
    ls = l;
    for (k = 0; k <= n - 2; k++)
    {
        p = 0.0;
        for (i = k; i <= ls; i++)
        {
            t = fabs(b[i * il]);
            if (t > p)
            {
                p = t;
                is = i;
            }
        }
        if (p < MIN)
        {
            printf("fail\n");
            return (0);
        }
        for (j = 0; j <= m - 1; j++)
        {
            u = k * m + j;
            v = is * m + j;
            t = d[u];
            d[u] = d[v];
            d[v] = t;
        }
        for (j = 0; j <= il - 1; j++)
        {
            u = k * il + j;
            v = is * il + j;
            t = b[u];
            b[u] = b[v];
            b[v] = t;
        }
        for (j = 0; j <= m - 1; j++)
        {
            u = k * m + j;
            d[u] = d[u] / b[k * il];
        }
        for (j = 1; j <= il - 1; j++)
        {
            u = k * il + j;
            b[u] = b[u] / b[k * il];
        }
        for (i = k + 1; i <= ls; i++)
        {
            t = b[i * il];
            for (j = 0; j <= m - 1; j++)
            {
                u = i * m + j;
                v = k * m + j;
                d[u] = d[u] - t * d[v];
            }
            for (j = 1; j <= il - 1; j++)
            {
                u = i * il + j;
                v = k * il + j;
                b[u - 1] = b[u] - t * b[v];
            }
            u = i * il + il - 1;
            b[u] = 0.0;
        }
        if (ls != (n - 1))
            ls = ls + 1;
    }
    p = b[(n - 1) * il];
    if (fabs(p) < MIN)
    {
        printf("fail\n");
        return (0);
    }
    for (j = 0; j <= m - 1; j++)
    {
        u = (n - 1) * m + j;
        d[u] = d[u] / p;
    }
    ls = 1;
    for (i = n - 2; i >= 0; i--)
    {
        for (k = 0; k <= m - 1; k++)
        {
            u = i * m + k;
            for (j = 1; j <= ls; j++)
            {
                v = i * il + j;
                is = (i + j) * m + k;
                d[u] = d[u] - b[v] * d[is];
            }
        }
        if (ls != (il - 1))
            ls = ls + 1;
    }
    return 1;
}

/*
ldlt分解
a 存放对称系数矩阵
n 阶数
m 右端常数向量的组数
c 返回m组解向量
*/
int LDLT(double a[], int n, int m, double c[])
{
    int i, j, l, k, u, v, w, k1, k2, k3;
    double p;
    if (fabs(a[0]) < MIN)
    {
        printf("fail \n");
        return -1;
    }
    for (i = 1; i <= n - 1; i++)
    {
        u = i * n;
        a[u] = a[u] / a[0];
    }
    for (i = 1; i <= n - 2; i++)
    {
        u = i * n + i;
        for (j = 1; j <= i; j++)
        {
            v = n * i + j - 1;
            l = (j - 1) * n + j - 1;
            a[u] = a[u] - a[v] * a[v] * a[l];
        }
        p = a[u];
        if (fabs(p) < MIN)
        {
            printf("1 fail\n");
            return -1;
        }
        for (k = i + 1; k <= n - 1; k++)
        {
            u = k * n + i;
            for (j = 1; j <= i; j++)
            {
                v = k * n + j - 1;
                l = i * n + j - 1;
                w = (j - 1) * n + j - 1;
                a[u] = a[u] - a[v] * a[l] * a[w];
            }
            a[u] = a[u] / p;
        }
    }
    u = n * n - 1;
    for (j = 1; j <= n - 1; j++)
    {
        v = (n - 1) * n + j - 1;
        w = (j - 1) * n + j - 1;
        a[u] = a[u] - a[v] * a[v] * a[w];
    }
    p = a[u];
    if (fabs(p)+1.0 == 1.0)
    {
        printf("2 fail\n");
        return -1;
    }
    for (j = 0; j <= m - 1; j++)
    {
        for (i = 1; i <= n - 1; i++)
        {
            u = i * m + j;
            for (k = 1; k <= i; k++)
            {
                v = i * n + k - 1;
                w = (k - 1) * m + j;
                c[u] = c[u] - a[v] * c[w];
            }
        }
    }
    for (i = 1; i <= n - 1; i++)
    {
        u = (i - 1) * n + i - 1;
        for (j = i; j <= n - 1; j++)
        {
            v = (i - 1) * n + j;
            w = j * n + i - 1;
            a[v] = a[u] * a[w];
        }
    }
    for (j = 0; j <= m - 1; j++)
    {
        u = (n - 1) * m + j;
        c[u] = c[u] / p;
        for (k = 1; k <= n - 1; k++)
        {
            k1 = n - k;
            k3 = k1 - 1;
            u = k3 * m + j;
            for (k2 = k1; k2 <= n - 1; k2++)
            {
                v = k3 * n + k2;
                w = k2 * m + j;
                c[u] = c[u] - a[v] * c[w];
            }
            c[u] = c[u] / a[k3 * n + k3];
        }
    }
    return 1;
}
/*
*Cholesky分解
a 对称正定矩阵，返回时其上三角存放着分解后的U矩阵
n 阶数
m 右端组数
d 返回m组解向量
*/
int Cholesky(double a[],int n,int m,double d[])
{
    int i,j,k,u,v;
    if(a[0] <MIN || a[0]<0.0)
    {
        printf(" 1 fail \n");
        return -1;
    }
    a[0]=sqrt(a[0]);
    for(j =1;j <= n-1;j++)
    {
         a[j] = a[j]/a[0];
    }

    for (i = 1; i <= n - 1; i++)
    {
        u = i * n + i;
        for (j = 1; j <= i; j++)
        {
            v = (j - 1) * n + i;
            a[u] = a[u] - a[v] * a[v];
        }
        printf("a[u] %f \n", a[u]);
        if(a[u] < 0.0 )
        {
            printf("2 fail \n");
            return -1;
        }
        a[u] = sqrt(a[u]);
        if(i != (n-1))
        {
            for(j = i+1; j <= n-1;j++)
            {
                v = i*n+j;
                for(k =1;k<=i;k++)
                {
                    a[v] = a[v] - a[(k - 1) * n + i] * a[(k - 1) * n + j];
                }
                a[v] = a[v] / a[u];
            }
        }
    }
    for(j =0; j <=m-1;j++)
    {
        d[j] = d[j]/a[0];
        for(i =1;i <=n-1;i++)
        {
            u = i*n+i;
            v = i*m+j;
            for(k = 1; k <= i; k++)
            {
                d[v] = d[v]-a[(k-1)*n+i]*d[(k-1)*m+j];
            }
            d[v] =d[v]/a[u];
        }
    }
    for(j = 0; j <=m-1; j++)
    {
        u = (n-1)*m+j;
        d[u] = d[u] / a[n * n - 1];
        for (k = n - 1; k >= 1; k--)
        {
            u = (k-1)*m+j;
            for (i = k; i <= n - 1; i++)
            {
                v = (k - 1) * n + i;
                d[u] = d[u] - a[v] * d[i * m + j];
            }
            v = (k - 1) * n + k - 1;
            d[u] = d[u] / a[v];
        }

    }
    return 1;
}
/*
高斯赛尔德迭代法
a 存放系数矩阵
b 右端常量向量
n 阶数
x 接向量
eps 精度要求
*/
int GaussSeidel(double a[], double b[], int n, double x[], double eps)
{
    int i, j, u, v;
    double p, t, s, q;
    for (i = 0; i <= n - 1; i++)
    {
        u = i * n + i;
        p = 0.0;
        x[i] = 0.0;
        for (j = 0; j <= n - 1; j++)
            if (i != j)
            {
                v = i * n + j;
                p = p + fabs(a[v]);
            }
        if (p >= fabs(a[u]))
        {
            printf("fail\n");
            return -1;
        }
    }
    p = eps + 1.0;
    while (p >= eps)
    {
        p = 0.0;
        for (i = 0; i <= n - 1; i++)
        {
            t = x[i];
            s = 0.0;
            for (j = 0; j <= n - 1; j++)
                if (j != i)
                    s = s + a[i * n + j] * x[j];
            x[i] = (b[i] - s) / a[i * n + i];
            q = fabs(x[i] - t) / (1.0 + fabs(x[i]));
            if (q > p)
                p = q;
        }
    }
    return 1;
}