#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX 1.0		/* max for t */
#define MIN 0.0		/* min for t */

FILE *output10a, *output20a, *output10b, *output20b;
FILE *output10c, *output20c, *output10d, *output20d;
FILE *output10alog, *output20alog, *output10blog, *output20blog;

int main(int argc, char *argv[])
{
    int N=10; //to plh8os tou t

    double t, xa10t[N], fa10t, xb10t[N], fb10t, xc10t[N], xd10t[N];
    double x_euler_a10[N], x_euler_b10[N], x_euler_c10[N], x_euler_d10[N];
    double x_opteuler_a10[N], x_opteuler_b10[N], x_opteuler_c10[N], x_opteuler_d10[N];
    double k1, k2, k3, k4;
    int j;
    double h, step;

    h = (MAX-MIN)/N;
    step = h/2;

    //arxikes times
    t = 0;
    xa10t[0] = x_euler_a10[0] = x_opteuler_a10[0] = 1.0;
    xb10t[0] = x_euler_b10[0] = x_opteuler_b10[0] = 1.0;
    xc10t[0] = x_euler_c10[0] = x_opteuler_c10[0] = 0.5;
    xd10t[0] = x_euler_d10[0] = x_opteuler_d10[0] = -1.0;

    fa10t = exp(2*t) + (t/2);
    fb10t = - ((pow(t,2)+3)/(t-3));

    output10a=fopen("aN10.dat", "w");
    output10b=fopen("bN10.dat", "w");
    output10c=fopen("cN10.dat", "w");
    output10d=fopen("dN10.dat", "w");
    output10alog=fopen("aN10log.dat", "w");
    output10blog=fopen("bN10log.dat", "w");

    //ektypwsh arxikwn timwn sta antistoixa arxeia
    fprintf(output10a, "t(i)\t\t  Euler\t\t opt.Euler\t Runge\t\t f(i)\n");
    fprintf(output10a, "%f\t%f\t%f\t%f\t%f\n", t, x_euler_a10[0], x_opteuler_a10[0], xa10t[0], fa10t);

    fprintf(output10b, "t(i)\t\t  Euler\t\t opt.Euler\t Runge\t\t f(i)\n");
    fprintf(output10b, "%f\t%f\t%f\t%f\t%f\n", t, x_euler_b10[0], x_opteuler_b10[0], xb10t[0], fb10t);

    fprintf(output10c, "t(i)\t\t  Euler\t\t opt.Euler\t Runge\n");
    fprintf(output10c, "%f\t%f\t%f\t%f\n", t, x_euler_c10[0], x_opteuler_c10[0], xc10t[0]);

    fprintf(output10d, "t(i)\t\t  Euler\t\t opt.Euler\t Runge\n");
    fprintf(output10d, "%f\t%f\t%f\t%f\n", t, x_euler_d10[0], x_opteuler_d10[0], xd10t[0]);

    fprintf(output10alog, "t(i)\t\t  Euler\t\t opt.Euler\t Runge\n");
    fprintf(output10alog, "%f\t0.000000\t0.000000\t0.000000\n", t);

    fprintf(output10blog, "t(i)\t\t  Euler\t\t opt.Euler\t Runge\n");
    fprintf(output10blog, "%f\t0.000000\t0.000000\t0.000000\n", t);

    //to loop ypologismou
    for (j=1; j<=N; j++)
    {
        t = MIN + j*h;
        //gia thn 1h synarthsh (erwthma II)
        //times euler
        k1 = 0.5 - t + (2*x_euler_a10[j-1]);
        x_euler_a10[j] = x_euler_a10[j-1] + (h*k1);
        //times veltiwmenou euler
        k1 = 0.5 - t + (2*x_opteuler_a10[j-1]);
        k2 = 0.5 - (t+h) + (2*(x_opteuler_a10[j-1]+(h*k1)));
        x_opteuler_a10[j] = x_opteuler_a10[j-1] + (h/2)*(k1+k2);
        //times runge-kutta
        k1 = 0.5 - t + 2*xa10t[j-1];
        k2 = 0.5 - (t+step) + 2*(xa10t[j-1]+step*k1);
        k3 = 0.5 - (t+step) + 2*(xa10t[j-1]+step*k2);
        k4 = 0.5 - (t+h) + 2*(xa10t[j-1]+h*k3);
        xa10t[j] = xa10t[j-1] + (h/6)*(k1+2*k2+2*k3+k4);
        //pragmatikes times
        fa10t = exp(2*t) + (t/2);
        //ektypwsh timwn sto arxeio
        fprintf(output10a, "%f\t%f\t%f\t%f\t%f\n",t, x_euler_a10[j], x_opteuler_a10[j], xa10t[j], fa10t);
        //gia to posostiaio la8os
        fprintf(output10alog, "%f\t%f\t%f\t%f\n",t, (100*fabs(x_euler_a10[j]-fa10t))/fa10t, (100*fabs(x_opteuler_a10[j]-fa10t))/fa10t, (100*fabs(xa10t[j]-fa10t))/fa10t);
        //omoia gia thn 2h synarthsh (erwthma V)
        k1 = (pow(x_euler_b10[j-1],2)+(2*t*x_euler_b10[j-1]))/(3 + pow(t,2));
        x_euler_b10[j] = x_euler_b10[j-1] + h*k1;
        k1 = (pow(x_opteuler_b10[j-1],2)+(2*t*x_opteuler_b10[j-1]))/(3 + pow(t,2));
        k2 = (pow(x_opteuler_b10[j-1]+(h*k1),2)+(2*(t+h)*(x_opteuler_b10[j-1]+(h*k1))))/(3 + pow(t+h,2));
        x_opteuler_b10[j] = x_opteuler_b10[j-1] + (h/2)*(k1+k2);
        k1 = (pow(xb10t[j-1],2)+(2*t*xb10t[j-1]))/(3+pow(t,2));
        k2 = (pow(xb10t[j-1]+(step*k1),2)+(2*(t+step)*(xb10t[j-1]+(step*k2))))/(3+pow(t+step,2));
        k3 = (pow(xb10t[j-1]+(step*k2),2)+(2*(t+step)*(xb10t[j-1]+(step*k2))))/(3+pow(t+step,2));
        k4 = (pow(xb10t[j-1]+(h*k3),2)+(2*(t+h)*(xb10t[j-1]+(h*k3))))/(3+pow(t+h,2));
        xb10t[j] = xb10t[j-1] + (h/6)*(k1+2*k2+2*k3+k4);
        fb10t = - ((pow(t,2)+3)/(t-3));
        fprintf(output10b, "%f\t%f\t%f\t%f\t%f\n",t, x_euler_b10[j], x_opteuler_b10[j], xb10t[j], fb10t);
        //gia to posostiaio la8os
        fprintf(output10blog, "%f\t%f\t%f\t%f\n",t, (100*fabs(x_euler_b10[j]-fb10t))/fb10t, (100*fabs(x_opteuler_b10[j]-fb10t))/fb10t, (100*fabs(xb10t[j]-fb10t))/fb10t);
        //omoia gia thn 3h synarthsh (erwthma VI)
        k1 = pow(t,2) + pow(x_euler_c10[j-1],2);
        x_euler_c10[j] = x_euler_c10[j-1] + h*k1;
        k1 = pow(t,2) + pow(x_opteuler_c10[j-1],2);
        k2 = pow(t+h,2) + pow(x_opteuler_c10[j-1]+(h*k1),2);
        x_opteuler_c10[j] = x_opteuler_c10[j-1] + (h/2)*(k1+k2);
        k1 = pow(t,2) + pow(xc10t[j-1],2);
        k2 = pow(t+step,2) + pow(xc10t[j-1]+(step*k1),2);
        k3 = pow(t+step,2) + pow(xc10t[j-1]+(step*k2),2);
        k4 = pow(t+h,2) + pow(xc10t[j-1]+(h*k3),2);
        xc10t[j] = xc10t[j-1] + (h/6)*(k1+2*k2+2*k3+k4);
        fprintf(output10c, "%f\t%f\t%f\t%f\n",t, x_euler_c10[j], x_opteuler_c10[j], xc10t[j]);
        //omoia gia thn 4h synarthsh (erwthma VI)
        k1 = (pow(t,2) - pow(x_euler_d10[j-1],2))*sin(x_euler_d10[j-1]);
        x_euler_d10[j] = x_euler_d10[j-1] + h*k1;
        k1 = (pow(t,2) - pow(x_opteuler_d10[j-1],2))*sin(x_opteuler_d10[j-1]);
        k2 = (pow(t+h,2) - pow(x_opteuler_d10[j-1]+(h*k1),2))*sin(x_opteuler_d10[j-1]+(h*k1));
        x_opteuler_d10[j] = x_opteuler_d10[j-1] + (h/2)*(k1+k2);
        k1 = (pow(t,2) - pow(xd10t[j-1],2))*sin(xd10t[j-1]);
        k2 = (pow(t+step,2) - pow(xd10t[j-1]+(step*k1),2))*sin(xd10t[j-1]);
        k3 = (pow(t+step,2) - pow(xd10t[j-1]+(step*k2),2))*sin(xd10t[j-1]);
        k4 = (pow(t+h,2) - pow(xd10t[j-1]+(h*k3),2))*sin(xd10t[j-1]);
        xd10t[j] = xd10t[j-1] + (h/6)*(k1+2*k2+2*k3+k4);
        fprintf(output10d, "%f\t%f\t%f\t%f\n",t, x_euler_d10[j], x_opteuler_d10[j], xd10t[j]);
    }

    fclose(output10a);
    fclose(output10b);
    fclose(output10c);
    fclose(output10d);
    fclose(output10alog);
    fclose(output10blog);

    //h idia diadikasia gia N=20
    N=20;
    double xa20t[N], fa20t, xb20t[N], fb20t, xc20t[N],xd20t[N];
    double x_euler_a20[N], x_euler_b20[N], x_euler_c20[N], x_euler_d20[N];
    double x_opteuler_a20[N], x_opteuler_b20[N], x_opteuler_c20[N], x_opteuler_d20[N];

    h = (MAX-MIN)/N;

    t = 0;
    xa20t[0] = x_euler_a20[0] = x_opteuler_a20[0] = 1.0;
    xb20t[0] = x_euler_b20[0] = x_opteuler_b20[0] = 1.0;
    xc20t[0] = x_euler_c20[0] = x_opteuler_c20[0] = 0.5;
    xd20t[0] = x_euler_d20[0] = x_opteuler_d20[0] = -1.0;

    fa20t = exp(2*t) + (t/2);
    fb20t = - ((pow(t,2)+3)/(t-3));

    output20a=fopen("aN20.dat", "w");
    output20b=fopen("bN20.dat", "w");
    output20c=fopen("cN20.dat", "w");
    output20d=fopen("dN20.dat", "w");
    output20alog=fopen("aN20log.dat", "w");
    output20blog=fopen("bN20log.dat", "w");

    fprintf(output20a, "t(i)\t\t  Euler\t\t opt.Euler\t Runge\t\t f(i)\n");
    fprintf(output20a, "%f\t%f\t%f\t%f\t%f\n", t, x_euler_a20[0], x_opteuler_a20[0], xa20t[0], fa20t);

    fprintf(output20b, "t(i)\t\t  Euler\t\t opt.Euler\t Runge\t\t f(i)\n");
    fprintf(output20b, "%f\t%f\t%f\t%f\t%f\n", t, x_euler_b20[0], x_opteuler_b20[0], xb20t[0], fb20t);

    fprintf(output20c, "t(i)\t\t  Euler\t\t opt.Euler\t Runge\n");
    fprintf(output20c, "%f\t%f\t%f\t%f\n", t, x_euler_c20[0], x_opteuler_c20[0], xc20t[0]);

    fprintf(output20d, "t(i)\t\t  Euler\t\t opt.Euler\t Runge\n");
    fprintf(output20d, "%f\t%f\t%f\t%f\n", t, x_euler_d20[0], x_opteuler_d20[0], xd20t[0]);

    fprintf(output20alog, "t(i)\t\t  Euler\t\t opt.Euler\t Runge\n");
    fprintf(output20alog, "%f\t0.000000\t0.000000\t0.000000\n", t);

    fprintf(output20blog, "t(i)\t\t  Euler\t\t opt.Euler\t Runge\n");
    fprintf(output20blog, "%f\t0.000000\t0.000000\t0.000000\n", t);

    for (j=1; j<=N; j++)
    {
        t = MIN + j*h;

        k1 = 0.5 - t + (2*x_euler_a20[j-1]);
        x_euler_a20[j] = x_euler_a20[j-1] + (h*k1);
        k1 = 0.5 - t + (2*x_opteuler_a20[j-1]);
        k2 = 0.5 - (t+h) + (2*(x_opteuler_a20[j-1]+(h*k1)));
        x_opteuler_a20[j] = x_opteuler_a20[j-1] + (h/2)*(k1+k2);
        k1 = 0.5 - t + 2*xa20t[j-1];
        k2 = 0.5 - (t+step) + 2*(xa20t[j-1]+step*k1);
        k3 = 0.5 - (t+step) + 2*(xa20t[j-1]+step*k2);
        k4 = 0.5 - (t+h) + 2*(xa20t[j-1]+h*k3);
        xa20t[j] = xa20t[j-1] + (h/6)*(k1+2*k2+2*k3+k4);
        fa20t = exp(2*t) + (t/2);
        fprintf(output20a, "%f\t%f\t%f\t%f\t%f\n",t, x_euler_a20[j], x_opteuler_a20[j], xa20t[j], fa20t);
        //gia to posostiaio la8os
        fprintf(output20alog, "%f\t%f\t%f\t%f\n",t, (100*fabs(x_euler_a20[j]-fa20t))/fa20t, (100*fabs(x_opteuler_a20[j]-fa20t))/fa20t, (100*fabs(xa20t[j]-fa20t))/fa20t);

        k1 = (pow(x_euler_b20[j-1],2)+(2*t*x_euler_b20[j-1]))/(3 + pow(t,2));
        x_euler_b20[j] = x_euler_b20[j-1] + h*k1;
        k1 = (pow(x_opteuler_b20[j-1],2)+(2*t*x_opteuler_b20[j-1]))/(3 + pow(t,2));
        k2 = (pow(x_opteuler_b20[j-1]+(h*k1),2)+(2*(t+h)*(x_opteuler_b20[j-1]+(h*k1))))/(3 + pow(t+h,2));
        x_opteuler_b20[j] = x_opteuler_b20[j-1] + (h/2)*(k1+k2);
        k1 = (pow(xb20t[j-1],2)+(2*t*xb20t[j-1]))/(3+pow(t,2));
        k2 = (pow(xb20t[j-1]+(step*k1),2)+(2*(t+step)*(xb20t[j-1]+(step*k2))))/(3+pow(t+step,2));
        k3 = (pow(xb20t[j-1]+(step*k2),2)+(2*(t+step)*(xb20t[j-1]+(step*k2))))/(3+pow(t+step,2));
        k4 = (pow(xb20t[j-1]+(h*k3),2)+(2*(t+h)*(xb20t[j-1]+(h*k3))))/(3+pow(t+h,2));
        xb20t[j] = xb20t[j-1] + (h/6)*(k1+2*k2+2*k3+k4);
        fb20t = - ((pow(t,2)+3)/(t-3));
        fprintf(output20b, "%f\t%f\t%f\t%f\t%f\n",t, x_euler_b20[j], x_opteuler_b20[j], xb20t[j], fb20t);
        //gia to posostiaio la8os
        fprintf(output20blog, "%f\t%f\t%f\t%f\n",t, (100*fabs(x_euler_b20[j]-fb20t))/fb20t, (100*fabs(x_opteuler_b20[j]-fa20t))/fb20t, (100*fabs(xb20t[j]-fb20t))/fb20t);

        k1 = pow(t,2) + pow(x_euler_c20[j-1],2);
        x_euler_c20[j] = x_euler_c20[j-1] + h*k1;
        k1 = pow(t,2) + pow(x_opteuler_c20[j-1],2);
        k2 = pow(t+h,2) + pow(x_opteuler_c20[j-1]+(h*k1),2);
        x_opteuler_c20[j] = x_opteuler_c20[j-1] + (h/2)*(k1+k2);
        k1 = pow(t,2) + pow(xc20t[j-1],2);
        k2 = pow(t+step,2) + pow(xc20t[j-1]+(step*k1),2);
        k3 = pow(t+step,2) + pow(xc20t[j-1]+(step*k2),2);
        k4 = pow(t+h,2) + pow(xc20t[j-1]+(h*k3),2);
        xc20t[j] = xc20t[j-1] + (h/6)*(k1+2*k2+2*k3+k4);
        fprintf(output20c, "%f\t%f\t%f\t%f\n",t, x_euler_c20[j], x_opteuler_c20[j], xc20t[j]);

        k1 = (pow(t,2) - pow(x_euler_d20[j-1],2))*sin(x_euler_d20[j-1]);
        x_euler_d20[j] = x_euler_d20[j-1] + h*k1;
        k1 = (pow(t,2) - pow(x_opteuler_d20[j-1],2))*sin(x_opteuler_d20[j-1]);
        k2 = (pow(t+h,2) - pow(x_opteuler_d20[j-1]+(h*k1),2))*sin(x_opteuler_d20[j-1]+(h*k1));
        x_opteuler_d20[j] = x_opteuler_d20[j-1] + (h/2)*(k1+k2);
        k1 = (pow(t,2) - pow(xd20t[j-1],2))*sin(xd20t[j-1]);
        k2 = (pow(t+step,2) - pow(xd20t[j-1]+(step*k1),2))*sin(xd20t[j-1]);
        k3 = (pow(t+step,2) - pow(xd20t[j-1]+(step*k2),2))*sin(xd20t[j-1]);
        k4 = (pow(t+h,2) - pow(xd20t[j-1]+(h*k3),2))*sin(xd20t[j-1]);
        xd20t[j] = xd20t[j-1] + (h/6)*(k1+2*k2+2*k3+k4);
        fprintf(output20d, "%f\t%f\t%f\t%f\n",t, x_euler_d20[j], x_opteuler_d20[j], xd20t[j]);
    }

    fclose(output20a);
    fclose(output20b);
    fclose(output20c);
    fclose(output20d);
    fclose(output20alog);
    fclose(output20blog);

    system("PAUSE");
    return 0;
}




