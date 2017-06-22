#include <stdio.h>
#include <stdlib.h>
#include <math.h>

FILE *outputa;
FILE *outputb;
FILE *outputaerror, *outputberror;

typedef struct 
{
       double x;
       double x_;
       double error;
}values;

values rkf45(double h, double t);

int main(int argc, char *argv[])
{
  double a,b;
  double emin,emax;
  double h,t, f;
  int N,j;
  values value_final;
  
  a = 1.0;
  b = 2.56;
  N = 10;
  h = (b-a)/N;
  emax = pow(10,-6);
  emin = pow(10,-8);
  t = a;
  
  value_final.x =  2.0;
  value_final.error = 0.0;
  f = 1 + t + tan(t-1);
  
  outputa=fopen("aN.dat", "w");			 
  fprintf(outputa, "t(i)\t\t  x(i)\t\t e(i) \t\t f(i)\n");
  fprintf(outputa, "%f\t%f\t%f\t%f\n", t, value_final.x, value_final.error, f);

  outputaerror=fopen("aNerror.dat", "w");			 
  fprintf(outputaerror, "t(i)\t\terror\n");
  fprintf(outputaerror, "%f\t0.000000\n", t);
  
  for(j=0; j<N; j++)
  {
           if ((emin<=value_final.error) && (value_final.error<=emax))
           {
                t = t + h;
                value_final = rkf45(h,t);
                f = 1 + t + tan(t-1);
           }
     
           if(value_final.error<emin)
           {
                h = 2*h;
                t = t +h;
                value_final = rkf45(h,t);
                f = 1 + t + tan(t-1);
           }
       
           if(value_final.error>emax)
           {
                h = h/2;
                t = t + h;
                value_final = rkf45(h/2,t);
                f = 1 + t + tan(t-1);
           }
                    
           fprintf(outputa, "%f\t%f\t%f\t%f\n", t, value_final.x, value_final.error, f);
           fprintf(outputaerror, "%f\t%f\n", t, (100*fabs(value_final.x-f))/f);
  }
  
  fclose(outputa);
  fclose(outputaerror);
  
  emax = pow(10,-11);
  emin = pow(10,-13);
  h = (b-a)/N;
  t = a;
  
  value_final.x =  2.0;
  value_final.error = 0.0;
  f = 1 + t + tan(t-1);
  
  outputb=fopen("bN.dat", "w");			 
  fprintf(outputb, "t(i)\t\t  x(i)\t\t e(i) \t\t f(i)\n");
  fprintf(outputb, "%f\t%f\t%f\t%f\n", t, value_final.x, value_final.error, f);
  
  outputberror=fopen("bNerror.dat", "w");			 
  fprintf(outputberror, "t(i)\t\terror\n");
  fprintf(outputberror, "%f\t0.000000\n", t);

  for(j=0; j<N; j++)
  {
           if ((emin<=value_final.error) && (value_final.error<=emax))
           {
                t = t + h;
                value_final = rkf45(h,t);
                f = 1 + t + tan(t-1);
           }
     
           if(value_final.error<emin)
           {
                h = 2*h;
                t = t +h;
                value_final = rkf45(h,t);
                f = 1 + t + tan(t-1);
           }
       
           if(value_final.error>emax)
           {
                h = h/2;
                t = t + h;
                value_final = rkf45(h/2,t);
                f = 1 + t + tan(t-1);
           }
                    
           fprintf(outputb, "%f\t%f\t%f\t%f\n", t, value_final.x, value_final.error, f);
           fprintf(outputberror, "%f\t%f\n", t, (100*fabs(value_final.x-f))/f);
  }
  
  fclose(outputb);
  fclose(outputberror);
  
  system("PAUSE");	
  return 0;
}

values rkf45(double h, double t)
{
   values value;
   double k1,k2,k3,k4,k5,k6;
   
   value.x = 2.0;
   value.error = 0.0;
  
   k1 = 2 + pow(value.x-t-1,2);
   k2 = 2 + pow((value.x+((h/4.0)*k1))-(t+(h/4.0))-1.0,2);
   k3 = 2 + pow((value.x+((3.0/32.0)*h*k1)+((9.0/32.0)*h*k2))-(t+(3.0/8.0)*h)-1.0,2);
   k4 = 2 + pow((value.x+((1932.0/2197.0)*h*k1)-((7200.0/2197.0)*h*k2)+((7296.0/2197.0)*h*k3))-(t+(12.0/13.0)*h)-1.0,2);
   k5 = 2 + pow((value.x+((439.0/216.0)*h*k1)-(8*h*k2)+((3680.0/513.0)*h*k3)-((845.0/4104.0)*h*k4))-(t+h)-1.0,2);
   k6 = 2 + pow((value.x-((8.0/27.0)*h*k1)+(2*h*k2)-((3544.0/2565.0)*h*k3)-((1859.0/4104.0)*h*k4)-((11.0/40.0)*h*k5))-(t+(h/2))-1.0,2); 
  
   value.x_ = value.x +((25.0/216.0)*h*k1)+((1408.0/2565.0)*h*k3)+((2197.0/4104.0)*h*k4)-((h/5)*k5);
   value.x = value.x +((16.0/135.0)*h*k1)+((6656.0/12825.0)*h*k3)+((28561.0/56430.0)*h*k4)-((9.0/50.0)*h*k5)+((2.0/55.0)*h*k6); 
   value.error = fabs(value.x_-value.x);
   
   return value;
}
