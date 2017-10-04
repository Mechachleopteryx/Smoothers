#include <math.h>
#include <stdio.h>
#define MAX_PATTERNS 250
#define MAX_PATTERN_LENGTH 11

const double TOLERANCE = 10.0e-10;

double kernel(double t);
double gcv_error(double patterns[][MAX_PATTERN_LENGTH], int num_inputs, int num_patterns, double lambda);

double kernel(double t) {
  return exp(-pow(t,2)/2.0)/sqrt(2.0*M_PI);
}

void print_vector(double *vec, int dim){
  int i;
  for (i = 0; i != dim; ++i) {
    printf("%f ", vec[i]);
  }
  printf("\n");
}

double gcv_error(double patterns[][MAX_PATTERN_LENGTH], int num_inputs, int num_patterns, double lambda){
  // Note that for each pattern, patterns[i], the first num_inputs are the
  // inputs and the last element of the vector is the output.

  int division_by_zero = 0;

  // Compute the matrix x_ij = x_i - x_j; x_i being the vector of
  // inputs in pattern i, as well as K_ij= kernel(|x_i - x_j|/lambda)

  double v[10];
  double K[MAX_PATTERNS][MAX_PATTERNS];
  int i, j, k;
  for (i = 0; i != num_patterns; ++i) {
    for (j = 0; j != num_patterns; ++j) {
      for (k = 0; k != num_inputs; ++k) {
	v[k] = patterns[i][k] - patterns[j][k];
      }
      double norm = 0;
      for (k = 0; k != num_inputs; ++k) {
	norm += pow(v[k], 2);
      }
      norm = sqrt(norm);
      K[i][j]=kernel(norm/lambda);
    }
  }

  // compute l_ij = l_i(x_j), where l_i() is the ith "basis" function
  // for the model
  double l[MAX_PATTERNS][MAX_PATTERNS];
  for (i = 0; i != num_patterns; ++i) {
    double sum = 0;
    for (k = 0; k != num_patterns; ++k) {
      sum += K[i][k];
    }
    for (j = 0; j != num_patterns; ++j) {
      l[j][i] = K[i][j]/sum;
    }
  }

  // compute the in-sample error and simultaneously the effective
  // number of degrees of freedom
  double error = 0;
  double degrees = 0;
  for (i = 0; i != num_patterns; ++i) {
    double yi_hat = 0;
    for (j = 0; j != num_patterns; ++j) {
      double yj = patterns[j][num_inputs];
      yi_hat += l[j][i]*yj;
    }
    double yi = patterns[i][num_inputs];
    error += pow(yi_hat - yi, 2);
    degrees += l[i][i];
  }
  error = error/num_patterns;

  // protect against degrees to close to number of patterns (division
  // by zero problem in cross-validation error)
  if (abs(degrees - num_patterns) < TOLERANCE) {
    degrees = 1;
    division_by_zero = 1;
  }

  // debugging:
  /* for (i = 0; i != num_patterns; ++i) { */
  /*   printf("\n"); */
  /*   for (j = 0; j != num_patterns; ++j) { */
  /*     printf("%f ", l[j][i]); */
  /*   } */
  /* } */
  /* printf("\n"); */
  /* printf("degrees = %f\n", degrees); */

  // return the generalized cross validation error
  // to caller
  if (division_by_zero == 0) {
    return error/( pow(1.0 - degrees/num_patterns, 2) );
  } else {
    return 1.0;
      }
}

double test(double patterns[][MAX_PATTERN_LENGTH], int num_inputs, int num_patterns){

  int i, j;
  double ret = 0;
  for (i = 0; i != num_patterns; ++i){
    for (j = 0; j != num_inputs; ++j){
      ret += patterns[i][j];
    }
  }
  return ret;
}

int main() {
  double patterns[3][MAX_PATTERN_LENGTH];
  int i;
  for (i = 0; i != 3; ++i) {
    patterns[i][0] = 1.0*(i+1);
    patterns[i][1] = pow(1.0*(i+1), 2);
    printf("%f, %f \n", patterns[i][0], patterns[i][1]);
  }
  double result = gcv_error(patterns, 1, 3, 0.6);
  printf("The error is %f \n", result);
  return 0;
}
