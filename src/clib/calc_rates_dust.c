#ifdef  SMBH_RAD
#include <stdlib.h> 
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"
#include "phys_constants.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#define mass_h 1.67262171e-24   
#define pi_val 3.141592653589793
#define tiny 1.0e-20
#define huge 1.0e+20
#define tevk 1.1605e+4

extern int grackle_verbose;

int calc_rates_dust(chemistry_data *my_chemistry, chemistry_data_storage *my_rates)
{
  double s_Silicate = 3.20185; // bulk density [g/cc]
  double s_Graphite = 2.27949; // bulk density [g/cc]
  double a_dust = 0.1e-4;      // grain radius [cm]

  // Cross-section per unit grain mass
  my_rates->S_Silicate = 3.0 / (4.0 * s_Silicate * a_dust);
  my_rates->S_Graphite = 3.0 / (4.0 * s_Graphite * a_dust);

  double kp_Silicate[] =
  {  -1.056,   -0.995,   -0.936,   -0.878,   -0.820,   -0.764,   -0.709,   -0.638,   -0.570,   -0.504,
     -0.435,   -0.363,   -0.294,   -0.216,   -0.141,   -0.058,    0.027,    0.114,    0.206,    0.301,
      0.398,    0.497,    0.596,    0.694,    0.793,    0.891,    0.989,    1.088,    1.187,    1.286,
      1.387,    1.489,    1.593,    1.699,    1.808,    1.919,    2.032,    2.146,    2.257,    2.362,
      2.459,    2.546,    2.623,    2.691,    2.752,    2.806,    2.854,    2.895,    2.926,    2.946,
      2.953,    2.946,    2.925,    2.890,    2.843,    2.784,    2.714,    2.635,    2.547,    2.452,
      2.351,    2.243,    2.131,    2.014,    1.894,    1.771,    1.646,    1.520,    1.396,    1.299  };

  double kp_Graphite[] =
  {  -0.568,   -0.506,   -0.445,   -0.387,   -0.329,   -0.272,   -0.216,   -0.140,   -0.068,    0.000,
      0.074,    0.154,    0.228,    0.313,    0.393,    0.481,    0.567,    0.655,    0.745,    0.837,
      0.930,    1.026,    1.124,    1.223,    1.326,    1.431,    1.540,    1.653,    1.770,    1.892,
      2.018,    2.146,    2.273,    2.396,    2.512,    2.616,    2.707,    2.783,    2.843,    2.888,
      2.917,    2.931,    2.932,    2.920,    2.896,    2.863,    2.820,    2.770,    2.713,    2.652,
      2.588,    2.523,    2.458,    2.397,    2.343,    2.298,    2.266,    2.250,    2.251,    2.272,
      2.312,    2.372,    2.452,    2.552,    2.673,    2.816,    2.978,    3.156,    3.340,    3.523  };

  int iTd;
  int    kp_gr_NTd;
  double kp_gr_dTd;
  double kp_gr_Td0;

  kp_gr_NTd  =  70;
  kp_gr_dTd  =    0.050;
  kp_gr_Td0  =    0.000;

  my_rates->kp_gr_N      = malloc(sizeof(int));
  my_rates->kp_gr_Td     = malloc(kp_gr_NTd * sizeof(double));
  my_rates->kp_Silicate  = malloc(kp_gr_NTd * sizeof(double));
  my_rates->kp_Graphite  = malloc(kp_gr_NTd * sizeof(double));

  my_rates->kp_gr_Size = kp_gr_NTd;
  my_rates->kp_gr_N[0] = kp_gr_NTd;
  my_rates->kp_gr_dTd  = kp_gr_dTd;
  for(iTd = 0; iTd < kp_gr_NTd; iTd++) {
    my_rates->kp_gr_Td[iTd] = kp_gr_Td0 + (double)iTd * kp_gr_dTd;
  }
  for(iTd = 0; iTd < kp_gr_NTd; iTd++) {
    my_rates->kp_Silicate[iTd] = kp_Silicate[iTd];
    my_rates->kp_Graphite[iTd] = kp_Graphite[iTd];
  }

  return SUCCESS;
}
#endif
