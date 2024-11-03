#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#define N 3

typedef double(*sistema)();

double *LeVetor(char *nome, int *m)
{
  FILE *fp = fopen(nome, "r");
  int i;
  double *v;

  fscanf(fp, "%d", m);
  v = (double *) malloc(*m *sizeof(double));

  for(i=0; i<*m; i++)
    fscanf(fp, "%lf", &v[i]);
  
  puts("");
  return v;
}

void ImprimeVetor(double *v, int j)
{
  int i;
  for(i=0; i<j; i++) printf("x[%d] = %g\t", i, v[i]);
}
double NormaVetor(double *vec, int c, int p)
{
  int i;
  double max = 0;

  if(p == 0)
  {
    for(i = 0;i<c;i++) 
    {
      if(fabs(vec[i])>max) 
        max = fabs(vec[i]);
    }
  }
  return max;
}

double func1(double *x)
{
  return (cos(x[1]*x[2])/3.0)+(1.0/6.0);
}

double func2(double *x)
{
  return (sqrt(pow(x[0],2)+sin(x[2])+1.06)/9) -0.1;
}

double func3(double *x)
{
  //double pi = 4*atan(1.0) //Usei o M_PI;

  return -(1.0/20.0)*exp(-x[0]*x[1])-((10*M_PI-3)/60);
}

double PontoFixoNL(sistema f[],double x[])
{
  int i;
  double xa[N];

  for(i=0;i<N;i++)
  {
    xa[i] = x[i];
    x[i] = f[i](x);
    xa[i] = fabs(xa[i]-x[i]);
  }
  return NormaVetor(xa,N,0);
}

int main(int argc, char **argv) 
{
  int m, itera=0;
  double norma, tol=1e-6;
  double *x;
  sistema equacao[N]={func1,func2,func3};

  x = LeVetor(argv[1], &m);
  puts("Vetor chute:");
  ImprimeVetor(x,N);
  puts("");
  
  do
  {
    printf("\nIteração:[%d]\n ", itera);
    norma = PontoFixoNL(equacao, x);
    itera++;
    ImprimeVetor(x, m);
    printf("Norma: %g", norma);

  }while (norma > tol);
  
  puts("");
  return 0;
}