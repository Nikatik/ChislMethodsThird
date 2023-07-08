#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double fabs (double x);
void   Uf (double** U, double h, double t);
void   Us (double** U, double h, double t);
void   UUf (double** U, double h, double t);
void   UUs (double** U, double h, double t);

int main (int argc, char* argv[])
{
    clock_t timer;
    timer = clock();

    double       t = 0.1, h = 0.1;
    char         prt  = 0;
    unsigned int mode = 0;
    if (argc > 1) sscanf (argv[1], "%c", &prt);
    if (argc > 2)
    {
        sscanf (argv[2], "%u", &mode);
        mode = mode % 4;
        if (argc > 3) sscanf (argv[3], "%lf", &t);
        if (argc > 4) sscanf (argv[4], "%lf", &h);
    }
    t = fabs (t);
    h = fabs (h);

    unsigned int l = (unsigned int) (1. / t) + 1;
    unsigned int k = (unsigned int) (2. / h) + 1;
    double**     U;

    U = (double**) malloc (sizeof (double*) * (l + 1));
    for (unsigned int i = 0; i < l + 1; i++)
    {
        U[i] = (double*) malloc (sizeof (double) * k);
        for (unsigned int j = 0; j < k; j++) U[i][j] = (DBL_MAX - 1);
        U[i][0]     = 0;
        U[i][k - 1] = 1;
    }
    for (unsigned int j = 0; j < k; j++)
    {
        if (j * h < 1. + h / 4)
        {
            U[0][j] = 0;
            continue;
        }
        if (j * h > 1.25 + h / 4)
        {
            U[0][j] = 1;
            continue;
        }
        U[0][j] = 4. * (j * h - 1.);
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
            fprintf (outf, "%.17lf|%.17lf|%.17lf\n", (double) j * h - 1.,
                     U[l - 1][j], U[l][j]);
    }
    else
    {
        /*  printf ("t = %lf\n", (double) (l - 1) * t);
          for (unsigned int j = 0; j < k; j++)
              printf ("%20.17lf  |%22.17lf  |%22.17lf\n", (double) j * h - 1.,
                      U[l-1][j], U[l][j]);
      */
        for (unsigned int i = 0; i < l + 1; i++)
        {
            for (unsigned int j = 0; j < k; j++) printf ("%8.5lf  ", U[i][j]);
            printf ("\n");
        }
    }
    timer -= clock();
    printf ("%.5f seconds\n", ((double) -timer) / CLOCKS_PER_SEC);

    for (unsigned int i = 0; i <= l; i++) free (U[i]);
    free (U);

    return 0;
}

double fabs (double x) { return x < 0 ? -x : x; }

void Uf (double** U, double h, double t)
{
    unsigned int l = (unsigned int) (1. / t) + 1;
    unsigned int k = (unsigned int) (2. / h) + 1;
    unsigned int i, j = 1;

    for (j = 1; j < k - 1; j++)
    {
        U[1][j] =
            U[0][j] + t * (U[0][j + 1] - U[0][j - 1]) / (4. * h) +
            t * t * (U[0][j + 1] - 2 * U[0][j] + U[0][j - 1]) / (8. * h * h);
        if (j * h < 0.5 + h / 4)
        {
            U[l][j] = 0;
            continue;
        }
        if (j * h > 0.75 + h / 4)
        {
            U[l][j] = 1;
            continue;
        }
        U[l][j] = 4. * (j * h - 1.) + 2.;
    }

    for (i = 2; i < l; i++)
        for (j = 1; j < k - 1; j++)
            U[i][j] = t * (U[i - 1][j + 1] - U[i - 1][j - 1]) / (2. * h) +
                      U[i - 2][j];
}

void UUf (double** U, double h, double t)
{
    unsigned int l = (unsigned int) (1. / t) + 1;
    unsigned int k = (unsigned int) (2. / h) + 1;
    unsigned int i, j = 1;

    for (j = 1; j < k - 1; j++)
    {
        U[1][j] = U[0][j] *
                  (1. + t * (U[0][j + 1] - U[0][j - 1]) / (2. * h) -
                   t * t * U[0][j] * (U[0][j + 1] - 2 * U[0][j] + U[0][j - 1]) /
                       (2. * h * h));
        if (j * h < 0.625 + h / 4)
        {
            U[l][j] = 0;
            continue;
        }
        U[l][j] = 1;
    }

    for (i = 2; i < l; i++)
        for (j = 1; j < k - 1; j++)
            U[i][j] = t *
                          (U[i - 1][j + 1] * U[i - 1][j + 1] -
                           U[i - 1][j - 1] * U[i - 1][j - 1]) /
                          (2. * h) +
                      U[i - 2][j];
}
