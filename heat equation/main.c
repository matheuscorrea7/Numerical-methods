#include <stdio.h>
#include <stdlib.h>

double p(double x)
{   
    return 0;
}

double q(double x, double gamma)
{    
    return gamma;
}

double r(double x, double gamma, double TA)
{    
    return -gamma*TA;
}

double *Diferencafinita(double a, int N, double h, double alpha, double beta,double gamma, double TA)
{
  double h2, t;
  double *D, *U, *L, *R, *y, *l, *u, *v, *z;
    int i;
    
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
    
    t = a+h;
    D[1] = 2 + h2*q(t, gamma);
    U[1] = -1 + (h/2)*p(t);
    R[1] = -h2*r(t, gamma, TA) + (1 + h/2*p(t))*alpha;
    
    for(i = 2; i<N; i++){
        
        t += h;
        
        D[i] = 2+ h2*q(t, gamma);
        U[i] = -1+(h/2)*p(t);
        L[i] = -1-(h/2)*p(t);
        R[i] = -h2*r(t, gamma, TA );
    }
    
    t += h;
    D[N] = 2 + h2*q(t, gamma);
    L[N] = -1 - (h/2)*p(t);
    R[N] = -h2*r(t, gamma, TA) + (1 - h/2*p(t))*beta;
    
    //Crout
    
    l[1] = D[1];
    u[1] = U[1]/D[1];
    z[1] = R[1]/l[1];
    
    for(i=2; i<N; i++){
        
        l[i] = D[i] - L[i]*u[i-1];
        u[i] = U[i]/l[i];
        z[i] = (R[i] - L[i]*z[i-1])/l[i];
    }
    
    l[N] = D[N] - L[N]*u[N-1];
    z[N] = (R[N] - L[N]*z[N-1])/l[N];
    
    y[N] = z[N];
    y[0] = alpha;
    y[N+1] = beta;
    for(i=N-1; i>0; i--){
        y[i] = z[i] - u[i]*y[i+1];
    }
    
    return y;

} 

int main(int argc, char **argv) 
{
  int N, i;
  double **M, *x, a, b, h, t, TA, gamma, alpha, beta;

  a = atof(argv[1]);
  b = atof(argv[2]);
  N = atof(argv[3]); //Terceiro item no arquivo .replit
  alpha = atof(argv[4]);
  beta = atof(argv[5]);
  gamma = atof(argv[6]);
  TA = atof(argv[7]);
  h = (b-a)/(N+1);

  x = Diferencafinita(a, N, h, alpha, beta, gamma, TA);

  for(i=0;i<= N+1;i++)
    printf("%lf %lf\n", a+i*h, x[i]);

  return 0;
}