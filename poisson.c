#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 100
#define NY 100
#define NT 200000
#define PI 3.14159265979323

int main() 
{
    // Parameters
    int nx = NX;
    int ny = NY;
    int nt = NT;
    double xmin = 0;
    double xmax = 4*PI;
    double ymin = 0;
    double ymax = 4*PI;

    double dx = (xmax - xmin) / (nx - 1);
    double dy = (ymax - ymin) / (ny - 1);

    // Initialization
    double p[NY][NX] = {0};
    double pd[NY][NX] = {0};
    double b[NY][NX] = {0};
    double x[NX], y[NY];
    double X[NY][NX], Y[NY][NX];
    // Error calculation
    double err[NY][NX];
    double totalerr = 0;
    double meanerr = 0;
    double maxerr = 0;
    double maxp = 0;
    double maxpexact = 0;

    for (int i = 0; i < nx; i++) {
        x[i] = xmin + i * dx;
    }
    for (int j = 0; j < ny; j++) {
        y[j] = ymin + j * dy;
    }
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            X[j][i] = x[i];
            Y[j][i] = y[j];
        }
    }

    // Source
    int cx = 16;
    int cy = 16;
    b[ny / 4][nx / 4] = 100;
    b[3 * ny / 4][3 * nx / 4] = -100;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            b[j][i] = 2. * ((Y[j][i] * Y[j][i] - Y[j][i]) + (X[j][i] * X[j][i] - X[j][i]));
        }
    }
    double pexact[NY][NX];
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            pexact[j][i] = sin(cx * X[j][i]) * sin(cy * Y[j][i]);
        }
    }
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            b[j][i] = -1. * (cx * cx * pexact[j][i] + cy * cy * pexact[j][i]);
        }
    }

    for (int it = 0; it < nt; it++) {
        // Copy p to pd
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                pd[j][i] = p[j][i];
            }
        }

        for (int j = 1; j < ny - 1; j++) {
            for (int i = 1; i < nx - 1; i++) {
                p[j][i] = ((pd[j][i + 1] + pd[j][i - 1]) * dy * dy +
                            (pd[j + 1][i] + pd[j - 1][i]) * dx * dx -
                            b[j][i] * dx * dx * dy * dy) /
                           (2 * (dx * dx + dy * dy));
            }
        }

        for (int i = 0; i < nx; i++) {
            p[0][i] = 0;
            p[ny - 1][i] = 0;
        }
        for (int j = 0; j < ny; j++) {
            p[j][0] = 0;
            p[j][nx - 1] = 0;
        }
    }


    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            err[j][i] = p[j][i] - pexact[j][i];
            totalerr += abs(err[j][i])*abs(err[j][i]);
            if (fabs(err[j][i]) > maxerr) {
                maxerr = fabs(err[j][i]);
            }
            if (p[j][i] > maxp) {
                maxp = p[j][i];
            }
            if (pexact[j][i] > maxpexact) {
                maxpexact = pexact[j][i];
            }
        }
    }

    meanerr = totalerr / (nx * ny);
    printf("%d %d %f %f %f %f\n", nx,ny,maxerr, meanerr, maxp, maxpexact);

    return 0;
}
