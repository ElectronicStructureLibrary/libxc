#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <xc.h>

static double xc_trial_points[][5] = {
  /* rhoa      rhob    sigmaaa   sigmaab   sigmabb */
  {0.17E+01, 0.17E+01, 0.81E-11, 0.81E-11, 0.81E-11},
  {0.17E+01, 0.17E+01, 0.17E+01, 0.17E+01, 0.17E+01},
  {0.15E+01, 0.15E+01, 0.36E+02, 0.36E+02, 0.36E+02},
  {0.88E-01, 0.88E-01, 0.87E-01, 0.87E-01, 0.87E-01},
  {0.18E+04, 0.18E+04, 0.55E+00, 0.55E+00, 0.55E+00},
  {0.18E+04, 0.18E+04, 0.86E+04, 0.86E+04, 0.86E+04},
  {0.16E+04, 0.16E+04, 0.37E+10, 0.37E+10, 0.37E+10},
  {0.26E+00, 0.26E+00, 0.28E+00, 0.28E+00, 0.28E+00},
  {0.53E+05, 0.53E+05, 0.96E+05, 0.96E+05, 0.96E+05},
  {0.47E+05, 0.47E+05, 0.29E+14, 0.29E+14, 0.29E+14},
  {0.15E+00, 0.15E+00, 0.16E+00, 0.16E+00, 0.16E+00},
  {0.35E+01, 0.00E+00, 0.46E-10, 0.00E+00, 0.00E+00},
  {0.35E+01, 0.00E+00, 0.34E+01, 0.00E+00, 0.00E+00},
  {0.30E+01, 0.00E+00, 0.20E+03, 0.00E+00, 0.00E+00},
  {0.58E-01, 0.00E+00, 0.47E-01, 0.00E+00, 0.00E+00},
  {0.82E+02, 0.81E+02, 0.49E+07, 0.49E+07, 0.49E+07},
  {0.39E+02, 0.38E+02, 0.81E+06, 0.82E+06, 0.82E+06},
  {0.13E+00, 0.95E-01, 0.15E+00, 0.18E+00, 0.22E+00},
  {0.78E-01, 0.31E-01, 0.41E-02, 0.38E-02, 0.36E-02},
  {0.50E+02, 0.49E+02, 0.11E+06, 0.11E+06, 0.11E+06},
  {0.40E+02, 0.40E+02, 0.99E+05, 0.98E+05, 0.98E+05},
  {0.12E+00, 0.10E+00, 0.12E+00, 0.13E+00, 0.14E+00},
  {0.48E-01, 0.25E-01, 0.46E-02, 0.44E-02, 0.41E-02},
  {0.0, 0.0, 0.0, 0.0, 0.0}
};

typedef struct {
  int family;

  lda_type lda_func;
  gga_type gga_func;
} functionals_type;


void get_point(functionals_type *func, double point[5], double *e, double der[5])
{
  switch(func->family)
    {
    case XC_FAMILY_LDA:
      lda(&(func->lda_func), &(point[0]), e, &(der[0]), NULL);
      break;
    case XC_FAMILY_GGA:
      gga(&(func->gga_func), &(point[0]), &(point[2]),
	  e, &(der[0]), &(der[2]));
      break;
    }  
}

void first_derivative(functionals_type *func, double point[5], double der[5])
{
  int i;

  for(i=0; i<5; i++){
    const double delta = 5e-10;

    double dd, p[5], v[5];
    int j;
    
    dd = point[i]*delta;
    if(dd < delta) dd = delta;

    for(j=0; j<5; j++) p[j] = point[j];

    if(point[i]>=3.0*dd){ /* centered difference */
      double em1, em2, ep1, ep2;

      p[i] = point[i] + dd;
      get_point(func, p, &ep1, v);
      ep1 *= (p[0] + p[1]);

      p[i] = point[i] + 2*dd;
      get_point(func, p, &ep2, v);
      ep2 *= (p[0] + p[1]);

      p[i] = point[i] - dd;  /* backward point */
      get_point(func, p, &em1, v);
      em1 *= (p[0] + p[1]);

      p[i] = point[i] - 2*dd;  /* backward point */
      get_point(func, p, &em2, v);
      em2 *= (p[0] + p[1]);

      der[i]  = 1.0/2.0*(ep1 - em1);
      der[i] += 1.0/12.0*(em2 - 2*em1 + 2*ep1 - ep2);

      der[i] /= dd;

    }else{                   /* we use a 5 point forward difference */
      double e1, e2, e3, e4, e5;

      p[i] = point[i];
      get_point(func, p, &e1, v);
      e1 *= (p[0] + p[1]);

      p[i] = point[i] + dd;
      get_point(func, p, &e2, v);
      e2 *= (p[0] + p[1]);     /* convert to energy per unit volume */

      p[i] = point[i] + 2.0*dd;
      get_point(func, p, &e3, v);
      e3 *= (p[0] + p[1]);

      p[i] = point[i] + 3.0*dd;
      get_point(func, p, &e4, v);
      e4 *= (p[0] + p[1]);

      p[i] = point[i] + 4.0*dd;
      get_point(func, p, &e5, v);
      e5 *= (p[0] + p[1]);

      der[i]  =          (-e1 + e2);
      der[i] -=  1.0/2.0*( e1 - 2*e2 + e3);
      der[i] +=  1.0/3.0*(-e1 + 3*e2 - 3*e3 + e4);
      der[i] -=  1.0/4.0*( e1 - 4*e2 + 6*e3 - 4*e4 + e5);

      der[i] /= dd;
    }
  }
}


void print_error(char *what, double diff, functionals_type *func, double *p)
{
  static char *red="\033[31;1m", *norm="\033[0m";
  char *color;

  color = (diff > 5e-4) ? red : norm;
  
  printf("%s: %s%g%s\n", what, color, diff, norm);

  if(func != NULL){
    int j;
    double e, v_an[5], v_fd[5];

    printf("   point (% 8.2e, % 8.2e, % 8.2e, % 8.2e, % 8.2e)\n", 
	   p[0], p[1], p[2], p[3], p[4]);

    for(j=0; j<5; j++){
      v_fd[j] = v_an[j] = 0.0;
    }

    get_point(func, p, &e, v_an);
    first_derivative(func, p, v_fd);

    printf("  analit (% 8.2e, % 8.2e, % 8.2e, % 8.2e, % 8.2e)\n", 
	   v_an[0], v_an[1], v_an[2], v_an[3], v_an[4]);
    printf("      fd (% 8.2e, % 8.2e, % 8.2e, % 8.2e, % 8.2e)\n", 
	   v_fd[0], v_fd[1], v_fd[2], v_fd[3], v_fd[4]);
  }
}


void test_functional(int functional, int nspin)
{
  functionals_type func;
  const func_type *info;
  int i, j, p_max[5];
  double max_diff[5], avg_diff[5];

  /* initialize functional */
  func.family = family_from_id(functional);
  switch(func.family)
    {
    case XC_FAMILY_LDA:
      if(functional == XC_LDA_X)
	lda_x_init(&(func.lda_func), nspin, 3, 0);
      else
	lda_init(&(func.lda_func), functional, nspin);

      info = func.lda_func.func;
      break;
    case XC_FAMILY_GGA:
      gga_init(&(func.gga_func), functional, nspin);

      info = func.gga_func.func;
      break;
    default:
      fprintf(stderr, "Functional '%d' not found\n", functional);
      exit(1);
    }
  
  for(j=0; j<5; j++){
    avg_diff[j] = 0.0;

    p_max[j]    = 0;
    max_diff[j] = -1.0;
  }

  for(i=0; xc_trial_points[i][0]!=0.0; i++){
    double e, v_fd[5], v_an[5];

    for(j=0; j<5; j++){
      v_fd[j] = v_an[j] = 0.0;
    }

    /* first, get the analitic gradients */
    get_point(&func, xc_trial_points[i], &e, v_an);

    /* now get the numerical gradients */
    first_derivative(&func, xc_trial_points[i], v_fd);

    /* make statistics */
    for(j=0; j<5; j++){
      double diff = fabs(v_an[j] - v_fd[j]);

      avg_diff[j] += diff;
      if(diff > max_diff[j]){
	max_diff[j] = diff;
	p_max[j] = i;
      }

    }
  }

  for(j=0; j<5; j++){
    avg_diff[j] /= i;
  }

  /* print statistics */
  printf("Functional: %s\n", info->name);
  print_error("Avg. error vrho", (avg_diff[0] + avg_diff[1])/2.0, NULL, NULL);
  j = (max_diff[0] > max_diff[1]) ? 0 : 1;
  print_error("Max. error vrho", max_diff[j], &func, xc_trial_points[p_max[j]]);

  print_error("Avg. error vsig", (avg_diff[2] + avg_diff[3] + avg_diff[4])/3.0, NULL, NULL);
  j = (max_diff[2] > max_diff[3]) ? 2 : 3;
  j = (max_diff[j] > max_diff[4]) ? j : 4;
  print_error("Max. error vsig", max_diff[j], &func, xc_trial_points[p_max[j]]);
}

/*----------------------------------------------------------*/
int main(int argc, char *argv[])
{
  if(argc != 2){
    printf("Usage:\n%s funct\n", argv[0]);
    return 1;
  }
  
  test_functional(atoi(argv[1]), 2);

  return 0;
}
