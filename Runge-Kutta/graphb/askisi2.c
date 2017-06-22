#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[])
{
  int N, j;
  double a, b, h, t, k;
  double x, y;
  double k1, k2, k3, k4, m1, m2, m3, m4;

  N = 100;
  a = 0.0;
  b = 0.7;
  
  h = (b-a)/N;
  t = a;
  k = 0.3;
  
  // arxikes times
  x = 1.0;
  y = 1.0;
  
  printf("t\t\t x\t\t y\n\n");
  
  for(j=0; j<=N; j++)
  {
           //ektypwsh apotelesmatwn
           printf("%f\t%f\t%f\n", t, x, y);
           //ypologismos lyshs systhmatos
           k1 = x*(pow(y+k,2)+((1-pow(k,2))*(1-x)));
           m1 = (k*x*((k*y)-1)) + (y*(pow(y+k,2) - 1 - pow(k,2)));
           k2 = (x+((h/2))*k1)*(pow((y+((h/2)*m1))+k,2)+((1-pow(k,2))*(1-(x+((h/2)*k1)))));
           m2 = (k*(x+((h/2)*k1))*((k*(y+((h/2)*m1)))-1)) + ((y+((h/2)*m1))*(pow((y+((h/2)*m1))+k,2) - 1 - pow(k,2)));
           k3 = (x+((h/2))*k2)*(pow((y+((h/2)*m2))+k,2)+((1-pow(k,2))*(1-(x+((h/2)*k2)))));
           m3 = (k*(x+((h/2)*k2))*((k*(y+((h/2)*m2)))-1)) + ((y+((h/2)*m2))*(pow((y+((h/2)*m2))+k,2) - 1 - pow(k,2)));
           k4 = (x+(h*k3))*(pow((y+(h*m3))+k,2)+((1-pow(k,2))*(1-(x+(h*k3)))));
           m4 = (k*(x+(h*k3))*((k*(y+(h*m3)))-1)) + (y*(pow((y+(h*m3))+k,2) - 1 - pow(k,2)));
           x = x + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
           y = y + (h/6)*(m1 + 2*m2 + 2*m3 + m4);
           //epomeno t
           t = t + h;
  }
  system("PAUSE");	
  return 0;
}
