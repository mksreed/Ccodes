#include<stdio.h>
void cyclic_thomas(const int X, double x[restrict X], const double a[restrict X], 
const double b[restrict X], const double c[restrict X], double tempcmod[restrict X], 
double v[restrict X]) 
{
    /* first solve a system of length X - 1 for two right hand sides, ignoring ix == 0 */
    tempcmod[1] = c[1] / b[1];
    v[1] = -a[1] / b[1];
    x[1] = x[1] / b[1];

    /* loop from 2 to X - 1 inclusive */
    for (int ix = 2; ix < X - 1; ix++) {
        const double m = 1.0 / (b[ix] - a[ix] * tempcmod[ix - 1]);
        tempcmod[ix] = c[ix] * m;
        v[ix] = (0.0f  - a[ix] * v[ix - 1]) * m;
        x[ix] = (x[ix] - a[ix] * x[ix - 1]) * m;
    }

    /* handle X - 1 */
    const double m = 1.0 / (b[X - 1] - a[X - 1] * tempcmod[X - 2]);
    tempcmod[X - 1] = c[X - 1] * m;
    v[X - 1] = (-c[0]    - a[X - 1] * v[X - 2]) * m;
    x[X - 1] = (x[X - 1] - a[X - 1] * x[X - 2]) * m;

    /* loop from X - 2 to 1 inclusive */
    for (int ix = X - 2; ix >= 1; ix--) {
        v[ix] -= tempcmod[ix] * v[ix + 1];
        x[ix] -= tempcmod[ix] * x[ix + 1];
    }

    x[0] = (x[0] - a[0] * x[X - 1] - c[0] * x[1]) / (b[0] + a[0] * v[X - 1] + c[0] * v[1]);

    /* loop from 1 to X - 1 inclusive */
    for (int ix = 1; ix < X; ix++)
        x[ix] += x[0] * v[ix];
}
int main()
{
  printf("\n Started thomas algorithm");
}
