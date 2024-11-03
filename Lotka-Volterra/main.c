#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 2 

typedef double (*equacoes)();
double f0(double x, double y[N], double alpha, double beta)
{
    
    return alpha*y[0] - beta*y[0]*y[1];
}

double f1(double x, double y[N], double delta, double gamma)
{
    
    return delta*y[0]*y[1] - gamma *y[1];
}

void rk4(double x, double y[N], double h, double alpha, double beta, double delta, double gamma)
{
    
    double k1[N], k2[N], k3[N], k4[N], yp[N], ypp[N], a[2][2] = {{alpha, beta}, {delta, gamma}};
    equacoes f[] = {f0,f1};
    int i;
    
    for(i = 0; i < N; i++)
    {    
        k1[i] = f[i](x, y, a[i][0], a[i][1]);
        yp[i] = y[i] + h/2*k1[i];
    }
    
    for(i = 0; i < N; i++)
    {    
        k2[i] = f[i](x+h/2, yp, a[i][0], a[i][1]);
        ypp[i] = y[i] + h/2*k2[i];
    }
    
    for(i = 0; i < N; i++)
    {
        k3[i] = f[i](x + h/2, ypp, a[i][0], a[i][1]);
        yp[i] = y[i] + h*k3[i];
    }
    
    for(i = 0; i < N; i++)
    {
        k4[i] = f[i](x+h, yp, a[i][0], a[i][1]);
        y[i] += h/6*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
    }
}

int main(int argc, char **argv)
{
    
    double a, b, h, x, y[N];
    int i;
    
    y[0] = atof(argv[1]);
    y[1] = atof(argv[2]);
    
    a = 0;
    b = 100;
    h = 0.1;
    x = a;
        
    for(x = a; x <= b; x+=h)
    {    
      printf("%g ", x);

      for(i = 0; i < N; i++) 
          printf("%g ", y[i]);
            
      rk4(x, y, h, 2./3., 4./3., 1., 1.0);
      puts("");
    }
      
}