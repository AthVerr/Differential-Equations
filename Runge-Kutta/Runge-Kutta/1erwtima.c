#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Pi 3.14
#define N 100

FILE *output; //orismos arxeiou gia ektypwsh grafikwn mesw gnuplot

int main(int argc, char *argv[])
{
  int j;
  double a, b, h, t;
  double x_exact, y_exact;
  double x_runge, y_runge;
  double k1, k2, k3, k4, m1, m2, m3, m4;

  a = 0.0;
  b = 2*Pi;
  
  h = (b-a)/N;
  t = a; 
  
  // arxikes times
  x_runge = -(1.0/3.0);
  y_runge = 2.0;
  
  //anoigma arxeiou
  output=fopen("my_arxeio.dat", "w");
  fprintf(output, "t\t\t x_exact\t y_exact\n");
  
  printf("t\t\t x_exact\t y_exact\t RK(x)\t\t Rk(y)\n\n");
  
  for(j=0; j<=N; j++)
  {
           //ypologismos akrivous ari8mhtikhs lyshs
           x_exact = (2.0*sin(3.0*t)/3.0) - (cos(3.0*t)/3.0);
           y_exact = (2.0*cos(3.0*t)) + (sin(3*t));
           //ektypwsh apotelesmatwn
           printf("%f\t%f\t%f\t%f\t%f\n", t, x_exact, y_exact, x_runge, y_runge);
           fprintf(output, "%f\t%f\t%f\n", t, x_exact, y_exact); 
           //ypologismos typou ringe-kutta
           k1 = y_runge;
           m1 = -9.0*x_runge;
           k2 = y_runge + ((h/2.0)*m1);
           m2 = -9.0*(x_runge + ((h/2.0)*k2));
           k3 = y_runge + ((h/2.0)*m2);
           m3 = -9.0*(x_runge + ((h/2.0)*m2));
           k4 = y_runge + (h*m3);
           m4 = -9.0*(x_runge + (h*m3));
           x_runge = x_runge + (h/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
           y_runge = y_runge + (h/6.0)*(m1 + 2.0*m2 + 2.0*m3 + m4);
           //epomeno t
           t = t + h;
  }
 
  fclose(output); //kleisimo arxeiou
  
  system("PAUSE");	
  return 0;
}
