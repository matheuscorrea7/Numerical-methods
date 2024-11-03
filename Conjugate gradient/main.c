#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double **LeMatriz(char *nome, int *m, int *n)
{
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
double NormaVetor(double *r, int m, int p)
{
  int i;
  double max = 0;

  if(p == 0)
  {
    for(i = 0;i<m;i++)
    {
      if(fabs(r[i]) > max) 
        max = fabs(r[i]);
    }
  }
  else
  {
    for(i =0 ;i<m;i++) 
      max += pow(fabs(r[i]),p); 
    
    max = pow(max,1.0/p); 
  }

  return max;
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
//Função produto de Matriz por vetor
double *MatrizxVetor(double **A, double *x, int nl, int nc)
{
  int i, j;
  double *vec;

  vec = (double*)malloc(nl*sizeof(double *));

  for(i=0; i<nl; i++)
  {
    for(j=0; j<nc; j++)
      vec[i]+= A[i][j]*(x[j]);
    
  }
  return vec;
}
//Fução produto interno
double prodint(double *x, double *y, int j)
{
  int i;
  double soma = 0;

  for(i=0;i<j;i++)
     soma += x[i] * y[i];
  
  return soma;
}
double *Residuo(double **A, double *x, int nl, int nc)
{
   int i;
   double *r, *Ax;

  //Calculando o residuo 
  r = (double*)malloc(nl*sizeof(double));
  Ax = MatrizxVetor(A,x,nl,nc);
  for(i=0;i<nl;i++)
    r[i] = A[i][nc - 1] - Ax[i];     
  
  //fim
  return r;
}
double gradienteconjugado(double **A, double *x, double *d,int nl, int nc, int p)
{
  int i;
  double *r, pr, *Ad, Add, alfa, beta=0, *rnovo;

 //Função que calcula resíduo
  r = rnovo = Residuo(A,x,nl,nc);
  //calculo do alfa em partes
  pr = prodint(r,r,nl);// produto de r com r
  Ad = MatrizxVetor(A,d,nl,nc); //produto de A com d
  Add = prodint(Ad,d,nl); //produto de Ad com d
  alfa = pr/(Add);
  /////fim

  for(i=0;i<nl;i++)
  {
    x[i]= x[i] + alfa*d[i];
    r[i]= r[i] - alfa*Ad[i];
  }

  //Calculando Beta
  beta = prodint(r,r,nl)/(pr);
  //fim

  //Calculando d^{i+1}
  for(i=0;i<nl;i++)
    d[i] = r[i] + beta*d[i];

  //printf("Produto interno entre rnovo e r é: %g", prodint(rnovo,r,nl));
  //exit(0);
  
  return NormaVetor(r, nl, p);
}

double gradiente(double **A, int nl, int nc, double *x, int p)
{
  int i,j;
  double *r, lambda, soma, pr, prr, *Ar;
 
  //Calculando o residuo 
  //r = malloc(nl*sizeof(double *));
  r = Residuo(A,x,nl,nc);

  //Calculando lambda separadamente
  pr = prodint(r,r,nl);
  Ar = MatrizxVetor(A, r, nl, nc);
  lambda = pr/prodint(r, Ar, nl);
  //fim do labda

  for(i=0;i<nl;i++)
    x[i]+=lambda*r[i];
  

  return NormaVetor(r, nl, p);
}



int main(int argc, char **argv) 
{
  double **A, *x, *d ,dx, tol=1e-8;
  int i,j, nl, nc, c, it = 0, w = 0;

  A = LeMatriz(argv[1], &nl, &nc);
  x = LeVetor(argv[2], &c);

  //Primeira condição
  d = Residuo(A, x, nl, nc);
  printf("it\t\tdx\n");
  do
  {
    it++;
    dx = gradienteconjugado(A, x, d, nl, nc, w);
    // dx = gradiente(A, nl, nc, x, w);
    //dx = Gauss(A, nl, nc, x, w);
    printf("%d %11.6g\n", it,dx);
  }while (dx > tol);

  return 0;
}