/* Created 2016-08-11 by Georgia Hansford
 * glh@student.unimelb.edu.au
 * Last update: 2016-08-22

 * This program evaluates a number of functions relevant to Project 1. First it
 * compares the Trapezoidal and Simpson's methods in numerically evaluating an
 * integral. It then uses Simpson's method and the false position algorithm to
 * determine the root of an integral equation. Putting these functions into
 * practice, it then determines the molecular vibrational energy levels for
 * two different potentials. First, a quadratic potential that can be evaluated
 * analytically, to test the program and gauge its accuracy. Then, the Morse
 * potential, a far more realistic function. The program outputs the calculated
 * energy levels for both potentials, and the known associated energies. It
 * optimises the variable 'a' in the Morse potential to fit the first and the
 * second energy levels of molecular hydrogen, and contrasts how the other
 * energy levels drift as a is varied.  */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define MAX_ITER 10000
#define MAX_ENERGY_LEVELS 15
#define MAX_OPTIMISATION_LEVELS  2
#define BORDER "--------------------------------------------------------------------------------\n"

/* Accuracy parameters */
int nsteps = 100000;
double start = 5.0;
double stepsize = 0.5;
double tolerance = 1.0e-7;

int nonzero_root_flag = 0;
int n = 0;
double curr_epsilon;
double a = 1;

/* Morse potential parameters */
double r_min = 0.74166;
double V0 = 4.747;
double gamma_coeff = 33.6567;
double energy_levels[15] = {-4.477, -3.962, -3.475, -3.017, -2.587, -2.185, -1.811, -1.466, -1.151, -.867, -.615, -.4, -.225, -.094, -0.017};


/* Procedure: calculates the value of the function log(1+x)/x for x>-1
 * Inputs:    a double, x
 * Outputs:   a double, the value of the function  */
double integrand441(double x)
{
    if (x<tolerance)
    {
	     return 1;
    }

    return log(1+x)/x;
}

/* Procedure: computes an integral using the trapezoidal method
 * Inputs:    two doubles, low and high, the terminals of the integral; an
 *            integer, numsteps, the number of divisions; and func, a function
 *            pointer to the integrand
 * Outputs:   a double, the numerical estimation of the integral  */
double trap_int(double low, double high, int numsteps, double (* func)(double))
{
    int i=0;
    double h = (high-low)/((double) numsteps);
    double x = low;

/* Trapezoidal formula: area = (x_a + 2*x_1 + 2*x_2 + ... + 2*x_n + x_b)*h/2 */
    double sum = func(x) + func(high);
    for (i=1; i<numsteps; i++)
    {
        x+=h;
        sum += 2*func(x);
    }

    return sum*h/2;
}

/* Procedure: computes an integral using Simpson's method
 * Inputs:    two doubles, low and high, the terminal of the integral; an
 *            integer, numsteps, the number of divisions; and func, a function
 *            pointer to the integrand
 * Outputs:   a double, the numerical estimation of the integral  */
double simp_int(double low, double high, int numsteps, double (* func)(double))
{
/* For Simpson's method, the number of steps must be even - check and fix here
 * if required. */
    if (numsteps%2!=0)
    {
        numsteps++;
    }

    int i;
    double h = (high-low)/((double) numsteps);
    double x = low;
    double sum = func(x) + func(high);

/* Simpson's formula: area = (x_a + 4*x_1 + 2*x_2 + 4*x_3 + 2*x_4 + ... + 2*x_n
 *                                + x_b)*h/3  */
    for (i=1; i<numsteps; i++)
    {
        x+=h;
        if (i%2==0)
        {
            sum += 2*func(x);
        }
        else
        {
            sum += 4*func(x);
        }
    }

    return sum*h/3;
}

/* Procedure: computes the value of the function x^2
 * Input:     a double, x
 * Output:    a double, the value of the function  */
double integrand442(double x)
{
    return x*x;
}

/* Procedure: computes the value of the integral function integral(t^2 from 0 to
 *            x) = x, using Simpson's method to determine the integral value
 * Input:     a double, x
 * Output:    a double, the value of the function  */
double function442(double x)
{
    return (simp_int(0, x, nsteps, integrand442) - x);
}

/* Procedure: computes the root of a function using a false-position algorithm.
 *            It first finds two points that straddle the root, and inputs them
 *            into a loop updating the points as per the algorithm. If the
 *            global 'nonzero_root_flag' is on, the function will not locate a
 *            root at zero, moving further away from the original guess if it
 *            gets too close to zero.
 * Inputs:    a double, guess, an initial guess of the root position; a double,
 *            fstep, a step size for locating points straddling the root; a
 *            function pointer, func, the function to find the root of
 * Outputs:   a double, the estimation of the root  */
double false_position(double guess, double fstep, double (*func)(double))
{
    int i;
    double xprev = guess;
    double yprev = func(xprev);
    double cstep = fstep;
    double xn = xprev + cstep;
    double yn = func(xn);

/* First, find two points that straddle the root (i.e. the function is positive
 * for one point, and negative for the other) by updating the previous guess by
 * the step size, terminating for a maximum number of iterations. When the
 * non zero root flag is on, if the next guess reaches the root at x=0, update
 * instead to a new guess further away from zero than the original guess.  */
    for(i=0; (i<MAX_ITER) && (yprev*yn>0); i++)
    {
  	     if (fabs(yn) > fabs(yprev))
  	     {
  	        cstep /= -2.0;
  	     }

        xprev = xn;
        xn += cstep;

        if(nonzero_root_flag)
        {
            if (xn<tolerance)
            {
                xn = guess + 2*i*fstep;
            }
        }

        yprev = yn;
        yn = func(xn);
    }

    if(i==MAX_ITER)
    {
        printf("Cannot locate straddling points. Early termination.\n");
        exit(EXIT_FAILURE);
    }

/* Now the false position algorithm can be used. Determine a new estimation
 * using the updating equation; select which of the previous points still
 * straddle the root; and repeat until a root is located.  */
    double xnext, ynext;

    while(fabs(xprev-xn)>tolerance)
    {
/* Updating the next point using the false position formula:
 * x_(i+1) = x_i - (f(x_i)*(x_i-x_(i-1))) / (f(x_i) - f(x_(i-1)))  */
	   xnext = xn - (yn*(xn-xprev))/(yn - yprev);
	   ynext = func(xnext);

    if ((yprev*ynext)>=0)
    {
        xprev = xn;
        yprev = yn;
    }

    xn = xnext;
    yn = ynext;

    if (fabs(yn)<tolerance)
    {
        break;
    }
  }

  return xn;
}

/* Procedure: computes the value of the normalised quadratic potential
 * Inputs:    a double, x, the normalised position
 * Outputs:   a double, the value of the potential  */
double normalised_quadratic_potential(double x)
{
    return 4*(x-1)*(x-2);
}

/* Procedure: computes the value of the integrand sqrt(epsilon - v(x)) for a
 *            quadratic potential
 * Inputs:    a double, x, the normalised position. curr_epsilon is a global
 *            variable that has already been updated for this evaluation
 * Outputs:   a double, the value of the integrand  */
double quadratic_integrand(double x)
{
    double arg = curr_epsilon - normalised_quadratic_potential(x);
/* The argument will need to be positive to prevent the square root from
 * returning -nan. If the argument is close to zero, force a return of 0.  */
    if (fabs(arg)<tolerance)
    {
        return 0;
    }

    return sqrt(arg);
}

/* Procedure: computes the action for a quadratic potential at an energy, using
 *            Simpson's method to determine the integral value. It uses the
 *            analytically determined functions for the terminals to calculate
 *            their value for this particular epsilon. It also updates the
 *            global curr_epsilon variable to be used in the integral.
 * Inputs:    a double, e, the normalised energy
 * Outputs:   a double, the value of the action  */
double quadratic_action(double e)
{
    double lower_term = .5*(3 - sqrt(1+e));
    double upper_term = .5*(3 + sqrt(1+e));
    curr_epsilon = e;

    return (simp_int(lower_term, upper_term, nsteps, quadratic_integrand) - (n+.5)*M_PI);
}

/* Procedure: computes the value of the normalised Morse potential
 * Inputs:    a double, x, the normalised position
 * Outputs:   a double, the value of the potential  */
double normalised_morse_potential(double x)
{
    return (pow(1-exp(-x + (r_min/a)), 2) -1);
}

/* Procedure: computes the value of an integrand
 * Inputs:    a double, x. curr_epsilon is used as a global variable.
 * Outputs:   a double, the input squared  */
double morse_integrand(double x)
{
    double arg = curr_epsilon - normalised_morse_potential(x);
/* The argument will need to be positive to prevent the square root from
 * returning -nan. If the argument is close to zero, force a return of 0.  */
    if (fabs(arg)<tolerance)
    {
        return 0;
    }

    return sqrt(arg);
}

/* Procedure: computes the value of the action for the Morse potential at a
 *            given energy, using Simpson's method to determine the intergral
 *            value. It uses the analytically determined functions for the
 *            terminals to calculate their value for this particular epsilon.
 *            It also updates the global curr_epsilon variable for the integral
 * Inputs:    a double, x
 * Outputs:   a double, the estimated value of the action  */
double morse_action(double e)
{
    double lower_term = (r_min/a) - log(1+sqrt(e + 1));
    double upper_term = (r_min/a) - log(1-sqrt(e + 1));
    curr_epsilon = e;

    return (gamma_coeff*a*simp_int(lower_term, upper_term, nsteps, morse_integrand) - (n+.5)*M_PI);
}

/* Procedure: computes the optimal value of a to match a given energy level of
 *            molecular hydrogen using the Morse potential
 * Inputs:    an int, level, the energy level to match to
 * Outputs:   a double, the value of a  */
void optimal_a(int level)
{
    n = level;
/* As epsilon only has the valid range (-1, 0), we encode strict accuracy
   parameters to prevent invalid responses.  */
    start = -0.5;
    stepsize = 0.01;

    double x = false_position(start, stepsize, morse_action);
    double a0=0, a1=0;
    double adj = 0.1;

    printf("Optimising a for n = %d ..... (this may take a moment)\n", n);
    while (fabs(x*V0 - (energy_levels[level])) > tolerance)
    {
        a0 = a1;
        a1 = a;

        if (fabs(x*V0) > fabs(energy_levels[level]))
        {
            a -= adj;
        }
        else
        {
            a += adj;
        }

        if (a0 == a)
        {
            adj /= 2;
        }

        x = false_position(start, stepsize, morse_action);
    }

    printf("Optimal value of a to fit energy level %d: %f\n", level, a);
}


int main(int argc, char **argv)
{
/* If a starting point and step size have been specified on the command line,
 * use these; otherwise use the default values 5.0 and .5  */
    if (argc>1)
    {
        nsteps = atoi(argv[1]);
        if(argc>2)
        {
            start = atof(argv[2]);
            stepsize = atof(argv[3]);
            tolerance = atof(argv[4]);
        }
    }

    printf("Running program with following parameters:\n");
    printf("Number of steps = %d; Start = %.4f; \nStepsize = %.4f; Tolerance = %.1e\n", nsteps, start, stepsize, tolerance);
    if (argc==1)
    {
        printf("If you would like to adjust these from the default values, please pass them as\ncommand line arguments\n");
    }

/* Stage 4.4.1 */
    printf("%s", BORDER);
    printf("Calculating the integral ln(1+x)/x from 0 to 1\n");
    printf("Trapezoidal method: %.10f; Simpson's method: %.10f\n", trap_int(0, 1, nsteps, integrand441), simp_int(0, 1, nsteps, integrand441));
    printf("Known value: %.10f\n", M_PI*M_PI/12);

/* Stage 4.4.2 */
    nonzero_root_flag = 1;
    printf("%s", BORDER);
    printf("Calculating the root of integral(t^2, from 0 to x) = x\n");
    printf("False position using Simpson's method for integral: %.6f\n", false_position(start, stepsize, function442));
    printf("Known value: %.6f\n", sqrt(3));

/* Quadratic potential stage */
    nonzero_root_flag = 0;
    printf("%s", BORDER);
    printf("Calculating energy levels for a quadratic potential V(r) = 4*V_0*(r/a - 1)*(r/a - 2)\n");
    for (n=0; n<MAX_ENERGY_LEVELS; n++)
    {
        printf("n = %d\t E = %9.6feV\t", n, false_position(start, stepsize, quadratic_action));
        if(n%2!=0)
        {
            printf("\n");
        }
    }

/* Morse potential stage */
    printf("\n%s", BORDER);

    int x;
    for(x=0; x<MAX_OPTIMISATION_LEVELS; x++)
    {
        optimal_a(x);
   	    for (n=0; n<15; n++)
        {
            double y = V0*false_position(start, stepsize, morse_action);
            printf("n = %d\t calculated E = %.3feV\t known E = %.3feV\t error = %5.1f%% \n", n, y, energy_levels[n], -fabs(y-energy_levels[n])*100/energy_levels[n]);
        }
    }

    return 0;
}
