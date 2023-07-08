#include <float.h>
#include <quadmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void Uf (__float128** U, __float128 h, __float128 t);
void Us (__float128** U, __float128 h, __float128 t);
void UUf (__float128** U, __float128 h, __float128 t);
void UUs (__float128** U, __float128 h, __float128 t);

int main (int argc, char* argv[])
{
    clock_t timer;
    timer = clock();

    long double  td = 0.1L, hd = 0.1L;
    char         prt  = 0;
    unsigned int mode = 0;
    if (argc > 1) sscanf (argv[1], "%c", &prt);
    if (argc > 2)
    {
        sscanf (argv[2], "%u", &mode);
        mode = mode % 4;
        if (argc > 3) sscanf (argv[3], "%Lf", &td);
        if (argc > 4) sscanf (argv[4], "%Lf", &hd);
    }
    __float128 t, h;
    t = fabsq ((__float128) td);
    h = fabsq ((__float128) hd);

    unsigned int l = (unsigned int) (1.Q / t) + 1;
    unsigned int k = (unsigned int) (2.Q / h) + 1;
    __float128** U;

    U = (__float128**) malloc (sizeof (__float128*) * (l + 1));
    for (unsigned int i = 0; i < l + 1; i++)
    {
        U[i] = (__float128*) malloc (sizeof (__float128) * k);
        for (unsigned int j = 0; j < k; j++) U[i][j] = (FLT128_MAX - 1);
        U[i][0]     = 0;
        U[i][k - 1] = 1;
    }
    for (unsigned int j = 0; j < k; j++)
    {
        if (j * h < 1.Q + h / 4.Q)
        {
            U[0][j] = 0;
            continue;
        }
        if (j * h > 1.25Q + h / 4.Q)
        {
            U[0][j] = 1;
            continue;
        }
        U[0][j] = 4.Q * (j * h - 1.Q);
    }
    switch (mode)
    {
        case (0): Uf (U, h, t); break;
        case (1): UUf (U, h, t); break;
        case (2):
            // Us (U, h, t);
            break;
        case (3):
            // UUs (U, h, t);
            break;
        default: return -1;
    }
    if (prt == '0')
    {
        FILE* outf = NULL;
        outf       = fopen ("trajectory.txt", "w");
        if (outf == NULL)
        {
            printf ("File for trajectories can`t be opened!\n");
            return -2;
        }
        for (unsigned int j = 0; j < k; j++)
            fprintf (outf, "%.40Lf|%.40Lf|%.40Lf\n",
                     (long double) (j * h - 1.Q), (long double) U[l - 1][j],
                     (long double) U[l][j]);
        fclose (outf);
    }
    else
    {
        /*  printf ("t = %lf\n", (__float128) (l - 1) * t);
          for (unsigned int j = 0; j < k; j++)
              printf ("%20.17lf  |%22.17lf  |%22.17lf\n", (__float128) j * h
          - 1., U[l-1][j], U[l][j]);
      */
        for (unsigned int i = 0; i < l + 1; i++)
        {
            for (unsigned int j = 0; j < k; j++)
                printf ("%8.5Lf  ", (long double) U[i][j]);
            printf ("\n");
        }
    }
    timer -= clock();
    printf ("%.5f seconds\n", ((double) -timer) / CLOCKS_PER_SEC);

    for (unsigned int i = 0; i <= l; i++) free (U[i]);
    free (U);

    return 0;
}

// __float128 fabs (__float128 x) { return x < 0 ? -x : x; }

void Uf (__float128** U, __float128 h, __float128 t)
{
    unsigned int l = (unsigned int) (1.Q / t) + 1;
    unsigned int k = (unsigned int) (2.Q / h) + 1;
    unsigned int i, j = 1;

    for (j = 1; j < k - 1; j++)
    {
        U[1][j] =
            U[0][j] + t * (U[0][j + 1] - U[0][j - 1]) / (4.Q * h) +
            t * t * (U[0][j + 1] - 2.Q * U[0][j] + U[0][j - 1]) / (8.Q * h * h);
        if (j * h < 0.5Q + h / 4.Q)
        {
            U[l][j] = 0;
            continue;
        }
        if (j * h > 0.75Q + h / 4.Q)
        {
            U[l][j] = 1;
            continue;
        }
        U[l][j] = 4.Q * (j * h - 1.Q) + 2.Q;
    }

    for (i = 2; i < l; i++)
        for (j = 1; j < k - 1; j++)
            U[i][j] = t * (U[i - 1][j + 1] - U[i - 1][j - 1]) / (2.Q * h) +
                      U[i - 2][j];
}

void UUf (__float128** U, __float128 h, __float128 t)
{
    unsigned int l = (unsigned int) (1.Q / t) + 1;
    unsigned int k = (unsigned int) (2.Q / h) + 1;
    unsigned int i, j = 1;
    __float128   dev = t / h;

    for (j = 1; j < k - 1; j++)
    {
        U[1][j] =
            U[0][j] * (1.Q + dev * (U[0][j + 1] - U[0][j - 1]) / 2.Q -
                       powq (dev, 2) * U[0][j] *
                           (U[0][j + 1] - 2.Q * U[0][j] + U[0][j - 1]) / 2.Q);
        if (j * h < 0.625Q + h / 4.Q)
        {
            U[l][j] = 0;
            continue;
        }
        U[l][j] = 1;
    }

    for (i = 2; i < l; i++)
        for (j = 1; j < k - 1; j++)
            U[i][j] =
                dev * (powq (U[i - 1][j + 1], 2) - powq (U[i - 1][j - 1], 2)) /
                    2.Q +
                U[i - 2][j];
}
