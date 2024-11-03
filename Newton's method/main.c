#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
double Jacobiana(double **J,  double *x)
{
  //Primeira linha
    J[0][0] = 4-x[3]; // derivada parcial em relação a x_{1}
    J[0][1] = -1; // derivada parcial em relação a x_{2}
    J[0][2] = 1; // derivada parcial em relação a x_{3}
    J[0][3] = -x[0]; // derivada parcial em relação a x_{4}
  //Segunda linha
    J[1][0] = -1; // derivada parcial em relação a x_{1}
    J[1][1] = 3-x[3]; // derivada parcial em relação a x_{2}
    J[1][2] = -2; // derivada parcial em relação a x_{3}
    J[1][3] = -x[1]; // derivada parcial em relação a x_{4}
  //Terceira linha
    J[2][0] = 1; // derivada parcial em relação a x_{1}
    J[2][1] = -2; // derivada parcial em relação a x_{2}
    J[2][2] = 3-x[3]; // derivada parcial em relação a x_{3}
    J[2][3] = -x[2]; // derivada parcial em relação a x_{4}
    //Quarta linha
    J[3][0] = 2*x[0]; // derivada parcial em relação a x_{1}
    J[3][1] = 2*x[1]; // derivada parcial em relação a x_{2}
    J[3][2] = 2*x[2]; // derivada parcial em relação a x_{3}
    J[3][3] = 0; // derivada parcial em relação a x_{4}
}

void TrocaLinha(double **k, double **l)
{
  double *temp;

  temp = *k;
  *k = *l;
  *l = temp;
}
void **Lusist(double **M, double ***L, double ***U, double **b, int m, int n)
{
  int i, j, k, p;
  double lbd, pivot;

  *L = calloc(m, sizeof(double *));
  for(i=0; i<m; i++) (*L)[i] = calloc(m, sizeof(double));

  *U = malloc(m*sizeof(double *));
  for(i=0; i<m; i++) (*U)[i] = malloc(m*sizeof(double));

  *b = malloc(m*sizeof(double));
  for(i=0; i<m; i++) (*b)[i] = M[i][n-1];

  for(i=0; i<m; i++)
  {
    (*L)[i][i] = 1;
    for(j=0; j<m; j++)
    {
      (*U)[i][j] = M[i][j];
    }
  }

  for(j=0; j<m-1; j++)
  {
    pivot = (*U)[j][j];
    p=0;

    for(i=j+1; i<m; i++){
      if(fabs((*U)[i][j])>pivot)
      {
        pivot=fabs((*U)[i][j]);
        p=i;
      }
    }

    if(p){
      TrocaLinha(&(*U)[p], &(*U)[j]);
      lbd=(*b)[p];
      (*b)[p]=(*b)[j];
      (*b)[j]=lbd;
    }

    for(i=j+1; i<m; i++){
      lbd=(*L)[i][j] = (*U)[i][j]/(*U)[j][j];
      for(k=j; k<m; k++) (*U)[i][k] -= lbd*(*U)[j][k];
    }

  }
}

double *SubstituicaoDireta(double **M, int m, double *b){
  double *y, sum;
  int i, j;

  y=malloc(m*sizeof(double));

  for(i=0; i<m; i++)
  {
    sum=0;
    for(j=i-1; j>=0; j--)
      sum += M[i][j]*y[j];
    
    y[i] = (b[i]-sum)/M[i][i];
  }
  return y;
}
double *SubstituicaoReversa(double **M, int m, double *b){
  double *y, sum;
  int i, j;

  y=malloc(m*sizeof(double));

  for(i=m-1; i>=0; i--){
    sum=0;
    for(j=i+1; j<=m-1; j++){
      sum += M[i][j]*y[j];
    }
    y[i] = (b[i]-sum)/M[i][i];
  }
  return y;
}

void ImprimeVetor(double *v, int j)
{
  int i;
  for(i=0; i<j; i++) printf("x[%d] = %g\n", i, v[i]);
}
double func1(double *x)
{
  return (4*x[0])-x[1]+x[2]-(x[0]*x[3]);
}

double func2(double *x)
{
  return -x[0]+(3*x[1])-(2*x[2])-(x[1]*x[3]);
}

double func3(double *x)
{
  return x[0]-(2*x[1])+(3*x[2])-(x[2]*x[3]);
}

double func4(double *x)
{
  return pow(x[0],2)+pow(x[1],2)+pow(x[2],2)-1;
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

int main(int argc, char **argv) 
{
  double **J, **L, **U, *x, *b, norma, tol=1e-7;
  int i, itera=0, N;
  double (*f[])(double*) = {func1,func2,func3,func4};//Isso é coisa do mal

  x = LeVetor(argv[1], &N);

  J = (double**)malloc(N*sizeof(double*));

  for(i=0;i<N;i++)
    J[i] = (double*)malloc((N+1)*sizeof(double));

  do
  {
    Jacobiana(J,x);

    for(i=0;i<N;i++)
      J[i][N] = -f[i](x);

    Lusist(J, &L, &U, &b, N, N+1);

    b = SubstituicaoDireta(L, N, b);
    b = SubstituicaoReversa(U, N, b);
    itera++;

    norma = NormaVetor(b, N, 0);

    for(i=0;i<N;i++)
      x[i]+=b[i];

    printf("\n%2d ||b||=%g\t x= ",itera, norma);
    ImprimeVetor(x, N);

  }while((norma>tol)&& (itera<20));

}