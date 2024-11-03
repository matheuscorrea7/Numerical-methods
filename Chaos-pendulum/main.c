#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#define N 2

//Parâmetros definidos
#define g 9.8 
#define L 9.8 
#define b 0.1 
#define fext 0.6
#define omega 0.66666

typedef double (*equacoes)();

//Função f0
double f0(double x, double y[N])
{ 
  return y[1];
}
//Função f1
double f1(double x, double y[N])
{
  return -(g/L)*sin(y[0])-b*y[1]-fext*cos(omega*x);
}

//Runge-Kutta
void rk4(double x, double y[N], equacoes f[N],double h)
{
  double k[4][N], a[N],a1[N];
  int i;

  for(i=0;i<N;i++)
  {
    k[0][i]=f[i](x,y);
    a[i]=y[i]+(h/2.0)*k[0][i];
  }
  for(i=0;i<N;i++)
  {
    k[1][i]=f[i](x+(h/2.0),a);
    a1[i]=y[i]+(h/2.0)*k[1][i];
  }
  for(i=0;i<N;i++)
  {
    k[2][i]=f[i](x+(h/2.0),a1);
    a[i]=y[i]+h*k[2][i];
  }
  for(i=0;i<N;i++)
  {
    k[3][i]=f[i](x+h,a);
    y[i]+=(h/6.0)*(k[0][i]+2.0*k[1][i]+2.0*k[2][i]+k[3][i]);
  }
}

int main() 
{
  double h=0.1,t, y1[N], y2[N];
  equacoes f[N]={f0,f1};

  //30 Graus
  y1[0]=(30*(4*atan(1.0))/180.0);
  y1[1]=0.0;

  //31 Graus
  y2[0]=(31*(4*atan(1.0))/180.0); 
  y2[1]=0.0;
  
  //|theta(30)-theta(31)| = fabs(yrk1[0]-yrk2[0])
  for(t=0;t<=10.0;t+=h) 
  {
    printf("%g %6g\n", t, fabs(y1[0]-y2[0]));
    rk4(t,y1,f,h);
    rk4(t,y2,f,h);
  }
  return 0;
}