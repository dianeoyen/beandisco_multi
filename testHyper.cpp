#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_result.h>
#include <cmath>
#include <stdio.h>

int main() {
  int maxU = 10;
  double** results = new double*[maxU];
  gsl_sf_result *out = new gsl_sf_result;
  int status;

  for (int i = 1; i <= maxU; i++) {
    results[i] = new double[i + 1];
    for (int j = 0; j <= i; j++) {
      status = gsl_sf_hyperg_2F1_e(1.0, (double)i, j + 2.0, 0.25, out);
      if (status) {
	printf("ERROR\t");
      } else {
	double result = out->val / (pow(4, i) * (j+1.0));
	printf("%f\t", result);
	results[i][j] = result;
      }
    }
    printf("\n");
  }
}
