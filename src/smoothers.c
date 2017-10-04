#include <math.h>
#include <stdio.h>
#define TOLERANCE 10.0e-10
#define MAX_PATTERNS 100000
#define MAX_PATTERN_LENGTH 2
// thats 1 input variable and 1 output variable per pattern

/* the one-dimensional gaussian kernel : */
double kernel(double t);

/* following needed for locally linear type of smoothing; called by l() below */
double g(int j, double patterns[][MAX_PATTERN_LENGTH], int num_patterns, double h, double x);

/* This function returns the <j>th basis function for linear
   smoothers, at point <x> when the bandwidth is <h>. It builds this
   function from the <patterns> passed. The type of smoothing is
   understood to be "locally constant" if type='0' and "locally linear" if
   type='1'. */
double l(int j, double patterns[][MAX_PATTERN_LENGTH], int num_patterns, double h, int type, double x);

/* This function computes the generalized cross validation error for a
   linear smoother built from the <patterns> passed, using a bandwidth
   of <h>. The type of smoothing is understood to be "locally
   constant" if type='0' and "locally linear" if type='1'. */
double Cgcv_error(double patterns[][MAX_PATTERN_LENGTH], int num_patterns, double h, int type);

/* This function returns the value at x of a linear smoother built
   from the patterns passed. <h> is the bandwidth and <type> the type
   of smoothing performed (see above) */
double Cprediction(double patterns[][MAX_PATTERN_LENGTH], int num_patterns, double h, int type, double x);

/* This function returns the effective number of degrees of freedom
   for a linear smoother built from the patterns passed. <h> is the
   bandwidth and <type> the type of smoothing performed; see above. */
double Ceffective_degrees(double patterns[][MAX_PATTERN_LENGTH], int num_patterns, double h, int type);

double kernel(double t) {
  return exp(-pow(t,2)/2.0)/sqrt(2.0*M_PI); //gaussian
  /* if (fabs(t) < 1) { */
  /*   return 70.0*pow(1-pow(fabs(t), 3), 3)/81.0; /\* tricube *\/ */
  /* } else { */
  /*   return 0; */
  /* } */
}

double g(int j, double patterns[][MAX_PATTERN_LENGTH], int num_patterns, double h, double x){
  double ret = 0;
  int i;
  for (i = 0; i != num_patterns; ++i) {
    ret += (patterns[i][0]-patterns[j][0])*(patterns[i][0]-x)*kernel((x-patterns[i][0])/h);
  }
  return ret*kernel((x-patterns[j][0])/h);
}

double l(int j, double patterns[][MAX_PATTERN_LENGTH], int num_patterns, double h, int type, double x){
    double denominator = 0;
    double numerator = 0;
    int i;
    if (type == 0) {
      for (i = 0; i != num_patterns; ++i) {
	double xi = patterns[i][0];
	double k = kernel( (x-xi)/h );
	denominator += k;
	if (i == j) {
	  numerator = k;
	}
      }
    } else {
      for (i = 0; i != num_patterns; ++i) {
	double xi = patterns[i][0];
	double k = g(i, patterns, num_patterns, h, x);
	denominator += k;
	if (i == j) {
	  numerator = k;
	}
      }
    }
    return numerator/denominator;
}

double Cgcv_error(double patterns[][MAX_PATTERN_LENGTH], int num_patterns, double h, int type){
  int division_by_zero = 0;
  double tolerance = TOLERANCE;

  // compute the in-sample error and simultaneously the effective
  // number of degrees of freedom
  double error = 0;
  double degrees = 0;
  int i, j;
  for (i = 0; i != num_patterns; ++i) {
    double xi = patterns[i][0];
    double yi = patterns[i][1];
    double yi_hat = 0;
    for (j = 0; j != num_patterns; ++j) {
      double yj = patterns[j][1];
      double lji = l(j, patterns, num_patterns, h, type, xi);
      yi_hat += lji*yj;
      if (i == j) {
	degrees += lji;
      }
    }
    error += pow(yi_hat - yi, 2);
  }
  error = error/num_patterns;

  // protect against degrees to close to number of patterns (division
  // by zero problem in cross-validation error)

  if (fabs(degrees - num_patterns) < tolerance) {
    degrees = 1;
    division_by_zero = 1;
  }

  // return the generalized cross validation error
  // to caller
  if (division_by_zero == 0) {
    return error/( pow(1.0 - degrees/num_patterns, 2) );
  } else {
    return 1.234;
      }
}

double Cprediction(double patterns[][MAX_PATTERN_LENGTH], 
		  int num_patterns, double h, int type, double x){
  int j;
  double y_hat = 0;
  for (j = 0; j != num_patterns; ++j) {
    double yj = patterns[j][1];
    y_hat += l(j, patterns, num_patterns, h, type, x)*yj;
  }
  return y_hat;
}

double Ceffective_degrees(double patterns[][MAX_PATTERN_LENGTH], int num_patterns, double h, int type){
  double degrees = 0;
  int i;
  for (i = 0; i != num_patterns; ++i) {
    degrees += l(i, patterns, num_patterns, h, type, patterns[i][0]);
  }
  return degrees;
}

int main() {
  double patterns[3][MAX_PATTERN_LENGTH];
  int i;
  for (i = 0; i != 3; ++i) {
    patterns[i][0] = 1.0*(i+1);
    patterns[i][1] = pow(1.0*(i+1), 2);
    printf("%f, %f \n", patterns[i][0], patterns[i][1]);
  }
  double result = Cgcv_error(patterns, 3, 0.6, 1);
  double pred = Cprediction(patterns, 3, 0.6, 1, 1);
  printf("The error is %f \n", result);
  printf("The prediction is %f", pred);
  return 0;
}
