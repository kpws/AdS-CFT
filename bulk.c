#include "bulk.h"
#include <math.h>

double f0(double z, double d0_1699, double d1_1700, double d2_1701, double d3_1702, double d4_1703, double d5_1704) {

   return (-d0_1699*pow(d2_1701, 2)*pow(z, 2) + 2*d0_1699*pow(z, 3) - 2*d0_1699 - d1_1700*pow(z, 7) - d1_1700*pow(z, 4) + 2*d1_1700*z)/(pow(z, 2)*(pow(z, 6) - 2*pow(z, 3) + 1));

}

double f1(double z, double d0_1699, double d1_1700, double d2_1701, double d3_1702, double d4_1703, double d5_1704) {

   return -2*pow(d0_1699, 2)*d2_1701/(pow(z, 2)*(pow(z, 3) - 1));

}

double f2(double z, double d0_1699, double d1_1700, double d2_1701, double d3_1702, double d4_1703, double d5_1704) {

   return (-2*pow(d0_1699, 2)*d4_1703*pow(z, 3) + 2*pow(d0_1699, 2)*d4_1703 - d4_1703*pow(z, 2) - 3*d5_1704*pow(z, 7) + 3*d5_1704*pow(z, 4))/(pow(z, 2)*(pow(z, 6) - 2*pow(z, 3) + 1));

}
