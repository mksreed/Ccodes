#include<stdio.h>
#include<math.h>
#define nx 40
//--------------------------------------------------------------------------------------
void cyclic_thomas_2(const int X, double x[restrict X], const double a[restrict X], 
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
//--------------------------------------------------------------------------------
void cyclic_thomas(const int X, double x[restrict X], const double a[restrict X], const double b[restrict X], const double c[restrict X], double cmod[restrict X], double u[restrict X]) {
    /*
     solves Ax = v, where A is a cyclic tridiagonal matrix consisting of vectors a, b, c
     X = number of equations
     x[] = initially contains the input v, and returns x. indexed from [0, ..., X - 1]
     a[] = subdiagonal, regularly indexed from [1, ..., X - 1], a[0] is lower left corner
     b[] = main diagonal, indexed from [0, ..., X - 1]
     c[] = superdiagonal, regularly indexed from [0, ..., X - 2], c[X - 1] is upper right
     cmod[], u[] = scratch vectors each of length X
     */

    /* lower left and upper right corners of the cyclic tridiagonal system respectively */
    const double alpha = a[0];
    const double beta = c[X - 1];

    /* arbitrary, but chosen such that division by zero is avoided */
    const double gamma = -b[0];

    cmod[0] = alpha / (b[0] - gamma);
    u[0] = gamma / (b[0] - gamma);
    x[0] /= (b[0] - gamma);

    /* loop from 1 to X - 2 inclusive */
    for (int ix = 1; ix + 1 < X; ix++) {
        const double m = 1.0 / (b[ix] - a[ix] * cmod[ix - 1]);
        cmod[ix] = c[ix] * m;
        u[ix] = (0.0f  - a[ix] * u[ix - 1]) * m;
        x[ix] = (x[ix] - a[ix] * x[ix - 1]) * m;
    }

    /* handle X - 1 */
    const double m = 1.0 / (b[X - 1] - alpha * beta / gamma - beta * cmod[X - 2]);
    u[X - 1] = (alpha    - a[X - 1] * u[X - 2]) * m;
    x[X - 1] = (x[X - 1] - a[X - 1] * x[X - 2]) * m;

    /* loop from X - 2 to 0 inclusive */
    for (int ix = X - 2; ix >= 0; ix--) {
        u[ix] -= cmod[ix] * u[ix + 1];
        x[ix] -= cmod[ix] * x[ix + 1];
    }

    const double fact = (x[0] + x[X - 1] * beta / gamma) / (1.0 + u[0] + u[X - 1] * beta / gamma);

    /* loop from 0 to X - 1 inclusive */
    for (int ix = 0; ix < X; ix++)
        x[ix] -= fact * u[ix];
}//-------------------------------------------------------------------------------
void thomas(const int X, double x[restrict X],
            const double a[restrict X], const double b[restrict X],
            const double c[restrict X], double scratch[restrict X]) {
    /*
     solves Ax = d, where A is a tridiagonal matrix consisting of vectors a, b, c
     X = number of equations
     x[] = initially contains the input v, and returns x. indexed from [0, ..., X - 1]
     a[] = subdiagonal, indexed from [1, ..., X - 1]
     b[] = main diagonal, indexed from [0, ..., X - 1]
     c[] = superdiagonal, indexed from [0, ..., X - 2]
     scratch[] = scratch space of length X, provided by caller, allowing a, b, c to be const
     not performed in this example: manual expensive common subexpression elimination
     */
    scratch[0] = c[0] / b[0];
    x[0] = x[0] / b[0];

    /* loop from 1 to X - 1 inclusive */
    for (int ix = 1; ix < X; ix++) {
        if (ix < X-1){
        scratch[ix] = c[ix] / (b[ix] - a[ix] * scratch[ix - 1]);
        }
        x[ix] = (x[ix] - a[ix] * x[ix - 1]) / (b[ix] - a[ix] * scratch[ix - 1]);
    }

    /* loop from X - 2 to 0 inclusive */
    for (int ix = X - 2; ix >= 0; ix--)
        x[ix] -= scratch[ix] * x[ix + 1];
}
//--------------------------------------------------------------------------------------
int main()
{
    FILE * fpout;
    int i=0,j=0;
    double pi=3.14159265358979323;
    double xmax=2*pi;
    double dx=xmax/(nx-1);
    double s1=8;
    printf("\n Started thomas algorithm");
    double a[nx],b[nx],c[nx],x[nx],d[nx];
    double ap[nx],bp[nx],cp[nx],xp[nx],dp[nx];
    double aas[nx],bas[nx],cas[nx],xas[nx],das[nx];
    double ss1[nx],ss2[nx];
    double u[nx],ux[nx],uxe[nx];
    for(i=0;i<nx;i++)
    {
        x[i]=i*dx;
        u[i]=cos(s1*x[i]);
        uxe[i]=-s1*sin(s1*x[i]);
    }
    fpout=fopen("outfile.txt","w");
    for (int i=0;i<nx;i++)
    {
        a[i]=1/3.;
        b[i]=1;
        c[i]=1/3.;
        //printf("\n%d,%f,%f,%f,%f,%f",i,x[i],a[i],b[i],c[i],d[i]);
        ap[i]=a[i];
        bp[i]=b[i];
        cp[i]=c[i];
        aas[i]=a[i];
        bas[i]=b[i];
        cas[i]=c[i];
    }
    for(i=2;i<nx-2;i++)
    {
        d[i]=14./9.*(u[i+1]-u[i-1])/2./dx+1./9.*(u[i+2]-u[i-2])/4./dx;
    }
    i=0;
    d[i]=14./9.*(u[i+1]-u[nx-2])/2./dx+1./9.*(u[i+2]-u[nx-3])/4./dx;
    i=1;
    d[i]=14./9.*(u[i+1]-u[nx-1])/2./dx+1./9.*(u[i+2]-u[nx-1])/4./dx;
    i=nx-1;
    d[i]=14./9.*(u[1]-u[i-1])/2./dx+1./9.*(u[2]-u[i-2])/4./dx;
    i=nx-2;
    d[i]=14./9.*(u[0]-u[i-1])/2./dx+1./9.*(u[1]-u[i-2])/4./dx;

    cyclic_thomas(nx-1,d,a,b,c,ss1,ss2);
    dp[nx-1]=dp[0];
    //thomas(nx,dp,ap,bp,cp,ss1);
   // thomas(nx,d,a,b,c,ss1);
   // thomas(nx,das,aas,bas,cas,ss1);
    fprintf(fpout,"i,x,a,b,c,ux,u,uxe");
    for (int i=0;i<=nx-1;i++)
    {
        printf("\n%d,%f,%f,%f,%f,%f,%f,%f",
                    i,x[i],a[i],b[i],c[i],d[i],u[i],uxe[i]);
        fprintf(fpout,"\n%d,%f,%f,%f,%f,%f,%f,%f",
                    i,x[i],a[i],b[i],c[i],d[i],u[i],uxe[i]);
    }
}
