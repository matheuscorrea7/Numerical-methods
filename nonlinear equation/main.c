#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f(double x, double y, double y1)
{
    return -pow(y1,2) - y + log(y);
}

double f1(double x, double y, double y1)
{
    return -1 + 1.0/y;
}

double f11(double x, double y, double y1)
{   
    return -2*y1;
}

double solexata(double x)
{
    return log(x);
}

double *DiferencaFinitaNL(double a, double b, int N, double h, double alpha, double beta){
    
    double x, t, h2, y1, norma;
    double *D, *U, *L, *R, *y, *l, *u, *v, *z;
    int i, itera;
    
    h2 = h*h;
    
    D = malloc((N+2)*sizeof(double));
    U = malloc((N+2)*sizeof(double));
    L = malloc((N+2)*sizeof(double));
    R = malloc((N+2)*sizeof(double));
    
    l = malloc((N+2)*sizeof(double));
    u = malloc((N+2)*sizeof(double));
    v = malloc((N+2)*sizeof(double));
    y = malloc((N+2)*sizeof(double));
    z = malloc((N+2)*sizeof(double));
    
    for(i=1; i<=N; i++) 
    	y[i] = alpha + i*h*(beta - alpha)/(b-a);
    
    y[0] = alpha;
    y[N+1] = beta;
    
    itera = 0;
    do
    {
     
        x = a + h;
        y1 = (y[2] - alpha)/(2*h);
        D[1] = 2 + h2*f1(x, y[1], y1);
        U[1] = -1 + (h/2)*f11(x, y[1], y1);
        R[1] = -(2*y[1] - y[2] - alpha + h2*f(x,y[1],y1));
            
        for(i=2; i<N; i++)
        {    
            x = a +i*h;
            y1 = (y[i+1] - y[i-1])/(2*h);
            D[i] = 2 + h2*f1(x, y[i], y1);
            L[i] = -1 - (h/2)*f11(x, y[i], y1);
            U[i] = -1 + (h/2)*f11(x, y[i], y1);
            R[i] = -(2*y[i] - y[i+1] - y[i-1] + h2*f(x, y[i], y1));
        }
        
        
        x = b - h;
        y1 = (beta - y[N-1])/(2*h);
        D[N] = 2 + h2*f1(x, y[N], y1);
        L[N] = -1 - (h/2)*f11(x, y[N], y1);
        R[N] = -(2*y[N] - y[N-1] - beta + h2*f(x, y[N], y1));
        
        //Algoritmo de Crout
        
        l[1] = D[1];
        u[1] = U[1]/D[1];
        z[1] = R[1]/l[1];
        
        for(i=2; i<N; i++)
        {    
            l[i] = D[i] - L[i]*u[i-1];
            u[i] = U[i]/l[i];
            z[i] = (R[i] - L[i]*z[i-1])/l[i];
        }
        
        l[N] = D[N] - L[N]*u[N-1];
        z[N] = (R[N] - L[N]*z[N-1])/l[N];
        
        v[N] = z[N];
        y[N] += v[N];
        norma = 0;
        
        for(i=N-1; i>0; i--)
        {    
            v[i] = z[i] - u[i]*v[i+1];
            if(fabs(v[i])>norma) norma = fabs(v[i]);
            y[i] += v[i];
        }
        itera++;
        
    }while((norma > 1e-8));
    
    return y;
}


int main(int argc, char **argv)
{
    
    int N, i, j;
    
    double **M, *x, a, b, h, t, alpha, beta;
    
    a = atof(argv[1]);
    b = atof(argv[2]);
    N = atof(argv[3]);
    alpha = atof(argv[4]);
    beta = log(atof(argv[5]));
    h = (b-a)/(N+1);
    
    x = DiferencaFinitaNL(a, b, N, h, alpha, beta);
    
    for(i=0; i<=N+1; i++)
    {
            printf("%lf %lf\n", a+i*h, x[i]-solexata(a+i*h));
            //fprintf(filedif, "%lf %lf\n", a+i*h, x[i] - solexata(a+i*h));
    }
        
    return 0;
}
