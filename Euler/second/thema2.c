#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX 2.56		/* max for t */
#define MIN 1.0		/* min for t */

FILE *output10a, *output20a, *output40a;		
FILE *output10aerror, *output20aerror, *output40aerror;

int main(int argc, char *argv[])
{
    int N=10;
    double t, h, xRKF10t[N], x_RKF_10t[N], f10t, xRK10t[N], error;
    double x_euler10[N], x_opteuler10[N];
    double k1, k2, k3, k4, k5, k6;
    int j;
    
    h = (MAX-MIN)/N;
    
    //arxikes times
    t = MIN;
    xRKF10t[0] = x_RKF_10t[0] = xRK10t[0] = x_euler10[0] = x_opteuler10[0] = 2.0;
    f10t = 1 + t + tan(t-1);
    
    output10a=fopen("aN10.dat", "w");
    output10aerror=fopen("aN10log.dat", "w");	
    			 
    fprintf(output10a, "t(i)\t\t  RK\t\t RKF \t\t RKF_\t\t error\t\t f(i)\n");
    fprintf(output10a, "%f\t%f\t%f\t%f\t%f\t%f\n", t, xRK10t[0], xRKF10t[0], x_RKF_10t[0], xRKF10t[0]-x_RKF_10t[0], f10t);
    
    fprintf(output10aerror, "t(i)\t\t  RK\t\t RKF \t\t RKF_\n");
    fprintf(output10aerror, "%f\t0.000000\t0.000000\t0.000000\n", t);
    
    printf("euler\t\topt.euler\n");
    printf("%f\t%f\n", x_euler10[0], x_opteuler10[0]);
    
    for (j=1; j<=N; j++)			
    {
        t = MIN + j*h;
        
        //oi times gia Euler
        k1 = 2 + pow(x_euler10[j-1]-t-1,2);
        x_euler10[j] = x_euler10[j-1] + h*k1;
        //oi times gia veltistopoihmeno Euler
        k1 = 2 + pow(x_opteuler10[j-1]-t-1,2);
        k2 = 2 + pow((x_opteuler10[j-1]+(h*k1))-(t+h)-1,2);
        x_opteuler10[j] = x_opteuler10[j-1] + (h/2)*(k1+k2);
        //oi times gia runge-kutta-fehleberg
        k1 = 2 + pow(xRKF10t[j-1]-t-1,2);
        k2 = 2 + pow((xRKF10t[j-1]+((h/4.0)*k1))-(t+(h/4.0))-1.0,2);
        k3 = 2 + pow((xRKF10t[j-1]+((3.0/32.0)*h*k1)+((9.0/32.0)*h*k2))-(t+(3.0/8.0)*h)-1.0,2);
        k4 = 2 + pow((xRKF10t[j-1]+((1932.0/2197.0)*h*k1)-((7200.0/2197.0)*h*k2)+((7296.0/2197.0)*h*k3))-(t+(12.0/13.0)*h)-1.0,2);
        k5 = 2 + pow((xRKF10t[j-1]+((439.0/216.0)*h*k1)-(8*h*k2)+((3680.0/513.0)*h*k3)-((845.0/4104.0)*h*k4))-(t+h)-1.0,2);
        k6 = 2 + pow((xRKF10t[j-1]-((8.0/27.0)*h*k1)+(2*h*k2)-((3544.0/2565.0)*h*k3)-((1859.0/4104.0)*h*k4)-((11.0/40.0)*h*k5))-(t+(h/2))-1.0,2); 
        //o typos 2
        x_RKF_10t[j] = xRKF10t[j-1] +((25.0/216.0)*h*k1)+((1408.0/2565.0)*h*k3)+((2197.0/4104.0)*h*k4)-((h/5)*k5);
        //o typos 3
        xRKF10t[j] = xRKF10t[j-1] +((16.0/135.0)*h*k1)+((6656.0/12825.0)*h*k3)+((28561.0/56430.0)*h*k4)-((9.0/50.0)*h*k5)+((2.0/55.0)*h*k6); 
        //o typos 4
        error = fabs(x_RKF_10t[j]-xRKF10t[j]);
        //oi pragmatikes times
        f10t = 1 + t + tan(t-1);
        //oi rimes gia runge kutta
        k1 = 2 + pow(xRK10t[j-1]-t-1,2);
        k2 = 2 + pow((xRK10t[j-1]+((h/2.0)*k1))-(t+(h/2.0))-1,2);
        k3 = 2 + pow((xRK10t[j-1]+((h/2.0)*k2))-(t+(h/2.0))-1,2);
        k4 = 2 + pow((xRK10t[j-1]+(h*k3))-(t+h)-1,2);
        xRK10t[j] = xRK10t[j-1] + (h/6.0)*(k1 + (2*k2) + (2*k3) + k4);
        fprintf(output10a, "%f\t%f\t%f\t%f\t%f\t%f\n", t,  xRK10t[j], xRKF10t[j], x_RKF_10t[j], error, f10t);
        //gia posostiaio la8os
        fprintf(output10aerror, "%f\t%f\t%f\t%f\n", t,  (100*fabs(xRK10t[j]-f10t))/f10t, (100*fabs(xRKF10t[j]-f10t))/f10t, (100*fabs(x_RKF_10t[j]-f10t))/f10t);            
        printf("%f\t%f\n", x_euler10[j], x_opteuler10[j]);
    }
  
    fclose(output10a);
    fclose(output10aerror);
    
    //h idia diadikasia gia N=20
    N=20;
    double xRKF20t[N], x_RKF_20t[N], f20t, xRK20t[N];  
    h = (MAX-MIN)/N;
    
    t = MIN;
    xRKF20t[0] = x_RKF_20t[0] = xRK20t[0] = 2.0;
    f20t = 1 + t + tan(t-1);
    
    output20a=fopen("aN20.dat", "w");			 
    fprintf(output20a, "t(i)\t\t  RK\t\t RKF \t\t RKF_\t\t error\t\t f(i)\n");
    fprintf(output20a, "%f\t%f\t%f\t%f\t%f\t%f\n", t, xRK20t[0], xRKF20t[0], x_RKF_20t[0], xRKF20t[0]-x_RKF_20t[0], f20t);
    
    output20aerror=fopen("aN20log.dat", "w");
    fprintf(output20aerror, "t(i)\t\t  RK\t\t RKF \t\t RKF_\n");
    fprintf(output20aerror, "%f\t0.000000\t0.000000\t0.000000\n", t);
    
    for (j=1; j<=N; j++)			
    {
        t = MIN + j*h;
        
        k1 = 2 + pow(xRKF20t[j-1]-t-1,2);
        k2 = 2 + pow((xRKF20t[j-1]+((h/4.0)*k1))-(t+(h/4.0))-1.0,2);
        k3 = 2 + pow((xRKF20t[j-1]+((3.0/32.0)*h*k1)+((9.0/32.0)*h*k2))-(t+(3.0/8.0)*h)-1.0,2);
        k4 = 2 + pow((xRKF20t[j-1]+((1932.0/2197.0)*h*k1)-((7200.0/2197.0)*h*k2)+((7296.0/2197.0)*h*k3))-(t+(12.0/13.0)*h)-1.0,2);
        k5 = 2 + pow((xRKF20t[j-1]+((439.0/216.0)*h*k1)-(8*h*k2)+((3680.0/513.0)*h*k3)-((845.0/4104.0)*h*k4))-(t+h)-1.0,2);
        k6 = 2 + pow((xRKF20t[j-1]-((8.0/27.0)*h*k1)+(2*h*k2)-((3544.0/2565.0)*h*k3)-((1859.0/4104.0)*h*k4)-((11.0/40.0)*h*k5))-(t+(h/2))-1.0,2); 
        x_RKF_20t[j] = xRKF20t[j-1] +((25.0/216.0)*h*k1)+((1408.0/2565.0)*h*k3)+((2197.0/4104.0)*h*k4)-((h/5)*k5);
        xRKF20t[j] = xRKF20t[j-1] +((16.0/135.0)*h*k1)+((6656.0/12825.0)*h*k3)+((28561.0/56430.0)*h*k4)-((9.0/50.0)*h*k5)+((2.0/55.0)*h*k6); 
        error = fabs(x_RKF_20t[j]-xRKF20t[j]);
        f20t = 1 + t + tan(t-1);
        k1 = 2 + pow(xRK20t[j-1]-t-1,2);
        k2 = 2 + pow((xRK20t[j-1]+((h/2.0)*k1))-(t+(h/2.0))-1,2);
        k3 = 2 + pow((xRK20t[j-1]+((h/2.0)*k2))-(t+(h/2.0))-1,2);
        k4 = 2 + pow((xRK20t[j-1]+(h*k3))-(t+h)-1,2);
        xRK20t[j] = xRK20t[j-1] + (h/6.0)*(k1 + (2*k2) + (2*k3) + k4);
        fprintf(output20a, "%f\t%f\t%f\t%f\t%f\t%f\n", t,  xRK20t[j], xRKF20t[j], x_RKF_20t[j], error, f20t);            
        //gia posostiaio la8os
        fprintf(output20aerror, "%f\t%f\t%f\t%f\n", t,  (100*fabs(xRK20t[j]-f20t))/f20t, (100*fabs(xRKF20t[j]-f20t))/f20t, (100*fabs(x_RKF_20t[j]-f20t))/f20t);
    }
  
    fclose(output20a);
    fclose(output20aerror);
    
    //h idia diadikasia gia N=40
    N=40;
    double xRKF40t[N], x_RKF_40t[N], f40t, xRK40t[N];  
    h = (MAX-MIN)/N;
    
    t = MIN;
    xRKF40t[0] = x_RKF_40t[0] = xRK40t[0] = 2.0;
    f40t = 1 + t + tan(t-1);
    
    output40a=fopen("aN40.dat", "w");			 
    fprintf(output40a, "t(i)\t\t  RK\t\t RKF \t\t RKF_\t\t error\t\t f(i)\n");
    fprintf(output40a, "%f\t%f\t%f\t%f\t%f\t%f\n", t, xRK40t[0], xRKF40t[0], x_RKF_40t[0], xRKF40t[0]-x_RKF_40t[0], f40t);
    
    output40aerror=fopen("aN40log.dat", "w");
    fprintf(output40aerror, "t(i)\t\t  RK\t\t RKF \t\t RKF_\n");
    fprintf(output40aerror, "%f\t0.000000\t0.000000\t0.000000\n", t);
    
    for (j=1; j<=N; j++)			
    {
        t = MIN + j*h;
        
        k1 = 2 + pow(xRKF40t[j-1]-t-1,2);
        k2 = 2 + pow((xRKF40t[j-1]+((h/4.0)*k1))-(t+(h/4.0))-1.0,2);
        k3 = 2 + pow((xRKF40t[j-1]+((3.0/32.0)*h*k1)+((9.0/32.0)*h*k2))-(t+(3.0/8.0)*h)-1.0,2);
        k4 = 2 + pow((xRKF40t[j-1]+((1932.0/2197.0)*h*k1)-((7200.0/2197.0)*h*k2)+((7296.0/2197.0)*h*k3))-(t+(12.0/13.0)*h)-1.0,2);
        k5 = 2 + pow((xRKF40t[j-1]+((439.0/216.0)*h*k1)-(8*h*k2)+((3680.0/513.0)*h*k3)-((845.0/4104.0)*h*k4))-(t+h)-1.0,2);
        k6 = 2 + pow((xRKF40t[j-1]-((8.0/27.0)*h*k1)+(2*h*k2)-((3544.0/2565.0)*h*k3)-((1859.0/4104.0)*h*k4)-((11.0/40.0)*h*k5))-(t+(h/2))-1.0,2); 
        x_RKF_40t[j] = xRKF40t[j-1] +((25.0/216.0)*h*k1)+((1408.0/2565.0)*h*k3)+((2197.0/4104.0)*h*k4)-((h/5)*k5);
        xRKF40t[j] = xRKF40t[j-1] +((16.0/135.0)*h*k1)+((6656.0/12825.0)*h*k3)+((28561.0/56430.0)*h*k4)-((9.0/50.0)*h*k5)+((2.0/55.0)*h*k6); 
        error = fabs(x_RKF_40t[j]-xRKF40t[j]);
        f40t = 1 + t + tan(t-1);
        k1 = 2 + pow(xRK40t[j-1]-t-1,2);
        k2 = 2 + pow((xRK40t[j-1]+((h/2.0)*k1))-(t+(h/2.0))-1,2);
        k3 = 2 + pow((xRK40t[j-1]+((h/2.0)*k2))-(t+(h/2.0))-1,2);
        k4 = 2 + pow((xRK40t[j-1]+(h*k3))-(t+h)-1,2);
        xRK40t[j] = xRK40t[j-1] + (h/6.0)*(k1 + (2*k2) + (2*k3) + k4);
        fprintf(output40a, "%f\t%f\t%f\t%f\t%f\t%f\n", t,  xRK40t[j], xRKF40t[j], x_RKF_40t[j], error, f40t);            
        //gia posostiaio la8os
        fprintf(output40aerror, "%f\t%f\t%f\t%f\n", t,  (100*fabs(xRK40t[j]-f40t))/f40t, (100*fabs(xRKF40t[j]-f40t))/f40t, (100*fabs(x_RKF_40t[j]-f40t))/f40t);
    }
  
    fclose(output40a);
    fclose(output40aerror);
    
    system("PAUSE");	
    return 0;
}


  


