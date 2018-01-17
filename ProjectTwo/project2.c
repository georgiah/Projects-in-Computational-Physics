/* Created 2016-08-25 by Georgia Hansford
 * glh@student.unimelb.edu.au
 * Last update: 2016-09-05

 * This program estimates solutions to the mass and density v. radius relations
 * for white dwarf stars, using Euler's method and the Runge-Kutta method to
 * solve the coupled differential equations. It first plots the solutions for
 * each method using a fixed step size of 0.01 and an initial density of 10.
 * Then it plots the solutions for a range of step sizes, in order to compare
 * and contrast the two methods. Going forward using just the Runge-Kutta
 * method, it plots the solutions for a range of central densities for a step
 * size of 0.00001. Finally, it allows for user input to manually vary the
 * values of the electron:nucleon ratio and the central density in order to
 * try and match the solutions to actual observed values for white dwarfs.  */


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "cpgplot.h"


#define MAX_STEPS      50000000
#define MAX_STEP_ERROR "Beware: calculations ended early as max array length has been reached\n"
#define TOLERANCE      1.0e-3
#define M0             5.67e33
#define R0             7.72e8
#define SOLAR_MASS     1.98e33
#define SOLAR_RADIUS   6.95e10

int stage=0;

/* Procedure: calculates the value of the function gamma
 * Inputs:    a double, x
 * Outputs:   a double, the value of the function  */
double gamma(double x)
{
    return x*x/(3*sqrt(1+x*x));
}

/* Procedure: calculates the value of dmdr
 * Inputs:    two doubles, the radius r and the density p
 * Outputs:   a double, the value of the derivative  */
double dmdr(double r, double p)
{
    return r*r*p;
}

/* Procedure: calculates the value of drhodr
 * Inputs:    three doubles, the radius r, the density p and the mass m
 * Outputs:   a double, the value of the derivative  */
double drhodr(double r, double p, double m)
{
/* If the density is less than zero, the cube root will return -nan, so take
 * the negative of the cube root of the absolute value if it is. */
    if(p<=0) {
	       return -m*p/(-gamma(pow(-p, 1/3.0))*r*r);
    }

    return -m*p/(gamma(pow(p, 1/3.0))*r*r);
}

/* Procedure: prints two graphs to the screen
 * Inputs:    three float arrays, representing the set of x values and the two
 *            sets of y values, and an int i, which is the array length
 * Outputs:   none  */
void print_two_graphs(float x[], float y1[], float y2[], int i) {
/* Set the text size, scale and labels for the first graph*/
    cpgsch(1.6);
    cpgenv(x[0], x[i-1], y1[0], y1[i-1], 0, 1);
    cpglab("Scaled radius, r", "Scaled mass, m", "dm/dr = r*r*rho");

    cpgbbuf();
    cpgline(i, x, y1);

/* Set the scale and labels for the second graph */
    cpgenv(x[0], x[i-1], y2[i-1], y2[0], 0, 1);
    cpglab("Scaled radius, r", "Scaled density, rho", "drho/dr = -m*rho/(gamma(cuberoot(rho))*r*r)");

    cpgline(i, x, y2);
}

/* Procedure: solves a set of coupled differential equations using Euler's
 *            method
 * Inputs:    the electron:nucleon ratio Ye, a double, the step size h, a
 *            double, the initial density p_c, a double, and two function
 *            pointers to the equations
 * Outputs:   none  */
void euler(double Ye, double h, double init_p) {
    // Initialise arrays, variables
    static float r[MAX_STEPS], m[MAX_STEPS], p[MAX_STEPS];
    r[0] = 0;
    m[0] = 0;
    p[0] = init_p;
    r[1] = h;
    m[1] = (float)h*h*h*init_p/3;
    p[1] = (float)init_p*(1-h*h*init_p/(6*gamma(pow(init_p, 1/3.0))));

    int i=1, j;

/* Increase r, calculating new values of m and p, until p = 0 */
    while(p[i]>TOLERANCE)
    {
	       i++;

/* Euler's method update equations */
        r[i] = (float)(r[i-1] + h);
        m[i] = (float)(m[i-1] + h*dmdr(r[i-1], p[i-1]));
        p[i] = (float)(p[i-1] + h*drhodr(r[i-1], p[i-1], m[i-1]));

        if(i>=MAX_STEPS -1)
        {
	           printf(MAX_STEP_ERROR);
             break;
        }
    }

    printf("Euler scaled mass: %e\t Euler scaled radius: %e\n", m[i-1], r[i-1]);

    print_two_graphs(r, m, p, i);
}

/* Procedure: solves a set of coupled differential equations using the
 *            Runge-Kutta method
 * Inputs:    the electron:nucleon ratio Ye, a double, the step size h, a
 *            double, the initial density p_c, a double, and two function
 *            pointers to the equations
 * Ouputs:    none  */
void runge_kutta(double Ye, double h, double init_p) {
    // Initialise arrays and variables
    static float r[MAX_STEPS], m[MAX_STEPS], p[MAX_STEPS];
    r[0] = 0;
    m[0] = 0;
    p[0] = init_p;
    r[1] = h;
    m[1] = (float)h*h*h*init_p/3;
    p[1] = (float)init_p*(1-h*h*init_p/(6*gamma(pow(init_p, 1/3.0))));

    int i=1, j;
    double f1, f2, f3, f4, g1, g2, g3, g4;

/* Increase r, calculating new values of m and p, until p=0 */
    while(p[i] > TOLERANCE)
    {
	       i++;

/* Runge-Kutta method update equations */
	       f1 = dmdr(r[i-1], p[i-1]);
	       g1 = drhodr(r[i-1], p[i-1], m[i-1]);

         f2 = dmdr(r[i-1] + h/2, p[i-1] + h*g1/2);
         g2 = drhodr(r[i-1] + h/2, p[i-1] + h*g1/2, m[i-1] + h*f1/2);

	       f3 = dmdr(r[i-1] + h/2, p[i-1] + h*g2/2);
	       g3 = drhodr(r[i-1] + h/2, p[i-1] + h*g2/2, m[i-1] + h*f2/2);

	       f4 = dmdr(r[i-1] + h, p[i-1] + h*g3);
	       g4 = drhodr(r[i-1] + h, p[i-1] + h*g3, m[i-1] + h*f3);

	       r[i] = (float)(r[i-1] + h);
	       m[i] = (float)(m[i-1] + h*(f1 + 2*f2 + 2*f3 + f4)/6);
	       p[i] = (float)(p[i-1] + h*(g1 + 2*g2 + 2*g3 + g4)/6);

	       if(i>=(MAX_STEPS-1))
	       {
	           printf(MAX_STEP_ERROR);
	           break;
	       }
    }

    printf("Scaled mass: %e\t Scaled radius: %e\n", m[i-1], r[i-1]);

    if(stage>=3)
    {
        printf("Solar mass: %e\t Solar radius: %e\n", m[i-1]*M0*Ye*Ye/SOLAR_MASS, r[i-1]*R0*Ye/SOLAR_RADIUS);
    }

    print_two_graphs(r, m, p, i);
}

int main(int argc, char **argv)
{
    double h = 0.01, Ye = 1, init_p = 10;

    if(cpgbeg(0, "?", 2, 1) != 1)
	       exit(EXIT_FAILURE);
    cpgask(1);

/* Stage 5.3.2 */
    printf("Calculating solutions for Euler's method using Ye = %.2lf , h = %lf  and p_c = %.2lf\n", Ye, h, init_p);
    euler(Ye, h, init_p);

    printf("Calculating solutions for Runge-Kutta method using Ye = %.2lf, h = %lf and p_c = %.2lf\n", Ye, h, init_p);
    runge_kutta(Ye, h, init_p);

    printf("Comparing methods for varying step size\n");
    for(h = 0.00001; h <= 0.1; h *= 10)
    {
	       printf("Calculating solutions for Euler's method using Ye = %.2lf , h = %lf  and p_c = %.2lf\n", Ye, h, init_p);
	       euler(Ye, h, init_p);

	       printf("Calculating solutions for Runge-Kutta method using Ye = %.2lf, h = %lf and p_c = %.2lf\n", Ye, h, init_p);
	       runge_kutta(Ye, h, init_p);
    }

/* Stage 5.3.3 */
    printf("Comparing central densities using Runge-Kutta method\n");

    h = 0.00001;
    stage = 3;

    for(init_p = 0.1; init_p <= 1.0e10; init_p *= 10)
    {
        printf("Calculating solutions for Runge-Kutta method using Ye = %.2lf, h = %lf and p_c = %.2lf\n", Ye, h, init_p);
        runge_kutta(Ye, h, init_p);
    }

/* Stage 5.3.5 */
    printf("Entering manual mode. Ctrl + C to exit\n");

    for(;;)
    {
        printf("Enter values for Ye, h and init_p: ");
        scanf("%lf %lf %lf", &Ye, &h, &init_p);

        runge_kutta(Ye, h, init_p);
    }

    cpgend();
    return 0;
}
