#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double **LeMatriz(char *nome, int *m, int *n){
  FILE *fp = fopen(nome, "r");
  int i, j;
  double **Matriz;

  fscanf(fp, "%d %d", m, n);
  Matriz = (double **) malloc(*m *sizeof(double *));

  for (i=0; i<*m; i++)
      Matriz[i] = (double *) malloc(*n *sizeof(double));
  
  for (i=0; i<*m; i++)
  {
    for (j=0; j<*n; j++)
      fscanf(fp, "%lf", &Matriz[i][j]);
    
  }
  puts("");
  return Matriz;
}
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
double NormaVetor(double *v, int m, int p)
{
  int i;
  double max = 0;

  if(p == 0)
  {
    for(i = 0;i<m;i++)
    {
      if(fabs(v[i])>max) 
        max = fabs(v[i]);
    }
  }
  else
  {
    for(i =0 ;i<m;i++) 
      max += pow(fabs(v[i]),p); 
    max = pow(max,1.0/p); 
  }
  return max;
}

double NormaMatriz(double **v, int m, int n, int p)
{
  double aux, max = 0;
  int i,j;

  if(p == 0)
  {
    for(i = 0;i<m; i++)
    {
      aux = 0;
      for(j = 0; j<n; j++)
        aux += fabs(v[i][j]);
      
      if(aux > max) 
        max = aux;
    }
  }
  else if(p == 1)
  {
    for(i = 0;i<n; i++)
    {
      aux = 0;
      for(j = 0;j<m;j++)
        aux += fabs(v[j][i]);
      
      if(aux > max) 
        max = aux;
    }
  }
  else{
    for(i = 0;i < m;i++)
    {
        for(j = 0;j < n;j++)
          max += pow(v[i][j],2);  
    }
    max = sqrt(max);
  }
  return max;
}

double Jacobi(double **Matriz, int m, int n, double *x0, int p)
{ 
  double *v, *r, sum;
  int i,j;

  v = (double *)malloc(m *sizeof(double));
  r = (double *) malloc(m *sizeof(double));

  for(i=0;i<m;i++)
  {
    sum = 0;
    for(j=0;j<n-1;j++)
    {
      if(i!=j) 
        sum += Matriz[i][j]*x0[j];
    }
    v[i]=(Matriz[i][n-1]-sum)/Matriz[i][i];
  }
  
  for(i=0;i<m;i++)
    r[i] = v[i] - x0[i]; 
  
  memcpy(x0, v, m*sizeof(double)); //Mudando a solução
  
  return NormaVetor(r, m, p);
}

double Gauss(double **Matriz, int m, int n, double *x0, int p)
{
  double *a, *b, sum;
  int i,j;

  a = (double *)malloc(m *sizeof(double));
  b = (double *) malloc(m *sizeof(double));
  memcpy(a, x0, m*sizeof(double));
  
  for(i = 0;i<m;i++)
  {
    sum = 0;
    for(j=0;j<n-1;j++)
    {
      if(i!=j) 
        sum += Matriz[i][j]*x0[j];
    }
    x0[i]=(Matriz[i][n-1]-sum)/Matriz[i][i]; 
  }
    
  for(i=0;i<m;i++)
    b[i] = a[i] - x0[i]; 
  
  return NormaVetor(b, m, p);
}

double Relaxacao(double **Matriz, int m, int n,double *x0, double omega,int p)
{
  double *a, *b, sum;
  int i,j;

  a = (double *)malloc(m *sizeof(double));
  b = (double *) malloc(m *sizeof(double));
  memcpy(a, x0, m*sizeof(double));
  
  if(omega >= 1 && omega < 2)
  {
    for(i=0;i<m;i++)
    {
      sum = 0;
      for(j=0;j<n-1;j++)
      {
        if(i!=j) 
          sum += Matriz[i][j]*x0[j];
      }
      x0[i] = ((1-omega)*x0[i]) + ((omega/Matriz[i][i])*(Matriz[i][n-1]-sum)); 
    }
  }
  else 
    printf("\nRaiz não foi encontrada\n");
  
  for(i=0;i<m;i++)
    b[i] = a[i] - x0[i];
  
  return NormaVetor(b, m, p);
}

 int main(int argc, char **argv) 
 {
    double **M, *v, dx, tolerance=1e-7;
    int m, n, l, i, it=0, p=0;

    M = LeMatriz(argv[1], &m, &n);//Matriz teste é a do exercicio 2
    v = LeVetor(argv[2], &l);
    for( i=0; i<m; i++) printf("%11.6g ", v[i]);
    puts("");

    do
    {
      it++;
      dx = Relaxacao(M, m, n, v,1.5, p);
      printf("%d  %8.4g ", it,dx);
      for( i=0; i<m; i++) printf("%11.6g", v[i]);
      puts("");
    } while (dx > tolerance);
 
    return 0;
 }