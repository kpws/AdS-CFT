#include "bulkODE.h"
#include <math.h>

double f1(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return (-d0_1704*pow(d2_1706, 2)*pow(z, 2) + 2*d0_1704*pow(z, 3) - 2*d0_1704 - d1_1705*pow(z, 7) - d1_1705*pow(z, 4) + 2*d1_1705*z)/(pow(z, 2)*(pow(z, 6) - 2*pow(z, 3) + 1));

}

double f3(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return -2*pow(d0_1704, 2)*d2_1706/(pow(z, 2)*(pow(z, 3) - 1));

}

double f5(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return (-18*pow(d0_1704, 2)*d4_1708*pow(z, 3) + 18*pow(d0_1704, 2)*d4_1708 - 18*pow(d0_1704, 2)*pow(z, 3)*cos(w*log(-z + 1)/3) + 18*pow(d0_1704, 2)*cos(w*log(-z + 1)/3) - 9*d4_1708*pow(w, 2)*pow(z, 2) - 27*d5_1709*pow(z, 7) + 27*d5_1709*pow(z, 4) + pow(w, 2)*pow(z, 6)*cos(w*log(-z + 1)/3) + 2*pow(w, 2)*pow(z, 5)*cos(w*log(-z + 1)/3) + 3*pow(w, 2)*pow(z, 4)*cos(w*log(-z + 1)/3) + 2*pow(w, 2)*pow(z, 3)*cos(w*log(-z + 1)/3) - 8*pow(w, 2)*pow(z, 2)*cos(w*log(-z + 1)/3) + 6*w*pow(z, 6)*sin(w*log(-z + 1)/3) + 3*w*pow(z, 5)*sin(w*log(-z + 1)/3) - 6*w*pow(z, 3)*sin(w*log(-z + 1)/3) - 3*w*pow(z, 2)*sin(w*log(-z + 1)/3))/(9*pow(z, 2)*(pow(z, 6) - 2*pow(z, 3) + 1));

}

double f7(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return (-18*pow(d0_1704, 2)*d6_1710*pow(z, 3) + 18*pow(d0_1704, 2)*d6_1710 - 18*pow(d0_1704, 2)*pow(z, 3)*sin(w*log(-z + 1)/3) + 18*pow(d0_1704, 2)*sin(w*log(-z + 1)/3) - 9*d6_1710*pow(w, 2)*pow(z, 2) - 27*d7_1711*pow(z, 7) + 27*d7_1711*pow(z, 4) + pow(w, 2)*pow(z, 6)*sin(w*log(-z + 1)/3) + 2*pow(w, 2)*pow(z, 5)*sin(w*log(-z + 1)/3) + 3*pow(w, 2)*pow(z, 4)*sin(w*log(-z + 1)/3) + 2*pow(w, 2)*pow(z, 3)*sin(w*log(-z + 1)/3) - 8*pow(w, 2)*pow(z, 2)*sin(w*log(-z + 1)/3) - 6*w*pow(z, 6)*cos(w*log(-z + 1)/3) - 3*w*pow(z, 5)*cos(w*log(-z + 1)/3) + 6*w*pow(z, 3)*cos(w*log(-z + 1)/3) + 3*w*pow(z, 2)*cos(w*log(-z + 1)/3))/(9*pow(z, 2)*(pow(z, 6) - 2*pow(z, 3) + 1));

}

double dfdz1(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return (-6*pow(z, 5) + 6*pow(z, 2))*(-d0_1704*pow(d2_1706, 2)*pow(z, 2) + 2*d0_1704*pow(z, 3) - 2*d0_1704 - d1_1705*pow(z, 7) - d1_1705*pow(z, 4) + 2*d1_1705*z)/(pow(z, 2)*pow(pow(z, 6) - 2*pow(z, 3) + 1, 2)) + (-2*d0_1704*pow(d2_1706, 2)*z + 6*d0_1704*pow(z, 2) - 7*d1_1705*pow(z, 6) - 4*d1_1705*pow(z, 3) + 2*d1_1705)/(pow(z, 2)*(pow(z, 6) - 2*pow(z, 3) + 1)) - 2*(-d0_1704*pow(d2_1706, 2)*pow(z, 2) + 2*d0_1704*pow(z, 3) - 2*d0_1704 - d1_1705*pow(z, 7) - d1_1705*pow(z, 4) + 2*d1_1705*z)/(pow(z, 3)*(pow(z, 6) - 2*pow(z, 3) + 1));

}

double dfdz3(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return 6*pow(d0_1704, 2)*d2_1706/pow(pow(z, 3) - 1, 2) + 4*pow(d0_1704, 2)*d2_1706/(pow(z, 3)*(pow(z, 3) - 1));

}

double dfdz5(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return (-6*pow(z, 5) + 6*pow(z, 2))*(-18*pow(d0_1704, 2)*d4_1708*pow(z, 3) + 18*pow(d0_1704, 2)*d4_1708 - 18*pow(d0_1704, 2)*pow(z, 3)*cos(w*log(-z + 1)/3) + 18*pow(d0_1704, 2)*cos(w*log(-z + 1)/3) - 9*d4_1708*pow(w, 2)*pow(z, 2) - 27*d5_1709*pow(z, 7) + 27*d5_1709*pow(z, 4) + pow(w, 2)*pow(z, 6)*cos(w*log(-z + 1)/3) + 2*pow(w, 2)*pow(z, 5)*cos(w*log(-z + 1)/3) + 3*pow(w, 2)*pow(z, 4)*cos(w*log(-z + 1)/3) + 2*pow(w, 2)*pow(z, 3)*cos(w*log(-z + 1)/3) - 8*pow(w, 2)*pow(z, 2)*cos(w*log(-z + 1)/3) + 6*w*pow(z, 6)*sin(w*log(-z + 1)/3) + 3*w*pow(z, 5)*sin(w*log(-z + 1)/3) - 6*w*pow(z, 3)*sin(w*log(-z + 1)/3) - 3*w*pow(z, 2)*sin(w*log(-z + 1)/3))/(9*pow(z, 2)*pow(pow(z, 6) - 2*pow(z, 3) + 1, 2)) + (-54*pow(d0_1704, 2)*d4_1708*pow(z, 2) - 6*pow(d0_1704, 2)*w*pow(z, 3)*sin(w*log(-z + 1)/3)/(-z + 1) + 6*pow(d0_1704, 2)*w*sin(w*log(-z + 1)/3)/(-z + 1) - 54*pow(d0_1704, 2)*pow(z, 2)*cos(w*log(-z + 1)/3) - 18*d4_1708*pow(w, 2)*z - 189*d5_1709*pow(z, 6) + 108*d5_1709*pow(z, 3) + pow(w, 3)*pow(z, 6)*sin(w*log(-z + 1)/3)/(3*(-z + 1)) + 2*pow(w, 3)*pow(z, 5)*sin(w*log(-z + 1)/3)/(3*(-z + 1)) + pow(w, 3)*pow(z, 4)*sin(w*log(-z + 1)/3)/(-z + 1) + 2*pow(w, 3)*pow(z, 3)*sin(w*log(-z + 1)/3)/(3*(-z + 1)) - 8*pow(w, 3)*pow(z, 2)*sin(w*log(-z + 1)/3)/(3*(-z + 1)) - 2*pow(w, 2)*pow(z, 6)*cos(w*log(-z + 1)/3)/(-z + 1) + 6*pow(w, 2)*pow(z, 5)*cos(w*log(-z + 1)/3) - pow(w, 2)*pow(z, 5)*cos(w*log(-z + 1)/3)/(-z + 1) + 10*pow(w, 2)*pow(z, 4)*cos(w*log(-z + 1)/3) + 12*pow(w, 2)*pow(z, 3)*cos(w*log(-z + 1)/3) + 2*pow(w, 2)*pow(z, 3)*cos(w*log(-z + 1)/3)/(-z + 1) + 6*pow(w, 2)*pow(z, 2)*cos(w*log(-z + 1)/3) + pow(w, 2)*pow(z, 2)*cos(w*log(-z + 1)/3)/(-z + 1) - 16*pow(w, 2)*z*cos(w*log(-z + 1)/3) + 36*w*pow(z, 5)*sin(w*log(-z + 1)/3) + 15*w*pow(z, 4)*sin(w*log(-z + 1)/3) - 18*w*pow(z, 2)*sin(w*log(-z + 1)/3) - 6*w*z*sin(w*log(-z + 1)/3))/(9*pow(z, 2)*(pow(z, 6) - 2*pow(z, 3) + 1)) - 2*(-18*pow(d0_1704, 2)*d4_1708*pow(z, 3) + 18*pow(d0_1704, 2)*d4_1708 - 18*pow(d0_1704, 2)*pow(z, 3)*cos(w*log(-z + 1)/3) + 18*pow(d0_1704, 2)*cos(w*log(-z + 1)/3) - 9*d4_1708*pow(w, 2)*pow(z, 2) - 27*d5_1709*pow(z, 7) + 27*d5_1709*pow(z, 4) + pow(w, 2)*pow(z, 6)*cos(w*log(-z + 1)/3) + 2*pow(w, 2)*pow(z, 5)*cos(w*log(-z + 1)/3) + 3*pow(w, 2)*pow(z, 4)*cos(w*log(-z + 1)/3) + 2*pow(w, 2)*pow(z, 3)*cos(w*log(-z + 1)/3) - 8*pow(w, 2)*pow(z, 2)*cos(w*log(-z + 1)/3) + 6*w*pow(z, 6)*sin(w*log(-z + 1)/3) + 3*w*pow(z, 5)*sin(w*log(-z + 1)/3) - 6*w*pow(z, 3)*sin(w*log(-z + 1)/3) - 3*w*pow(z, 2)*sin(w*log(-z + 1)/3))/(9*pow(z, 3)*(pow(z, 6) - 2*pow(z, 3) + 1));

}

double dfdz7(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return (-6*pow(z, 5) + 6*pow(z, 2))*(-18*pow(d0_1704, 2)*d6_1710*pow(z, 3) + 18*pow(d0_1704, 2)*d6_1710 - 18*pow(d0_1704, 2)*pow(z, 3)*sin(w*log(-z + 1)/3) + 18*pow(d0_1704, 2)*sin(w*log(-z + 1)/3) - 9*d6_1710*pow(w, 2)*pow(z, 2) - 27*d7_1711*pow(z, 7) + 27*d7_1711*pow(z, 4) + pow(w, 2)*pow(z, 6)*sin(w*log(-z + 1)/3) + 2*pow(w, 2)*pow(z, 5)*sin(w*log(-z + 1)/3) + 3*pow(w, 2)*pow(z, 4)*sin(w*log(-z + 1)/3) + 2*pow(w, 2)*pow(z, 3)*sin(w*log(-z + 1)/3) - 8*pow(w, 2)*pow(z, 2)*sin(w*log(-z + 1)/3) - 6*w*pow(z, 6)*cos(w*log(-z + 1)/3) - 3*w*pow(z, 5)*cos(w*log(-z + 1)/3) + 6*w*pow(z, 3)*cos(w*log(-z + 1)/3) + 3*w*pow(z, 2)*cos(w*log(-z + 1)/3))/(9*pow(z, 2)*pow(pow(z, 6) - 2*pow(z, 3) + 1, 2)) + (-54*pow(d0_1704, 2)*d6_1710*pow(z, 2) + 6*pow(d0_1704, 2)*w*pow(z, 3)*cos(w*log(-z + 1)/3)/(-z + 1) - 6*pow(d0_1704, 2)*w*cos(w*log(-z + 1)/3)/(-z + 1) - 54*pow(d0_1704, 2)*pow(z, 2)*sin(w*log(-z + 1)/3) - 18*d6_1710*pow(w, 2)*z - 189*d7_1711*pow(z, 6) + 108*d7_1711*pow(z, 3) - pow(w, 3)*pow(z, 6)*cos(w*log(-z + 1)/3)/(3*(-z + 1)) - 2*pow(w, 3)*pow(z, 5)*cos(w*log(-z + 1)/3)/(3*(-z + 1)) - pow(w, 3)*pow(z, 4)*cos(w*log(-z + 1)/3)/(-z + 1) - 2*pow(w, 3)*pow(z, 3)*cos(w*log(-z + 1)/3)/(3*(-z + 1)) + 8*pow(w, 3)*pow(z, 2)*cos(w*log(-z + 1)/3)/(3*(-z + 1)) - 2*pow(w, 2)*pow(z, 6)*sin(w*log(-z + 1)/3)/(-z + 1) + 6*pow(w, 2)*pow(z, 5)*sin(w*log(-z + 1)/3) - pow(w, 2)*pow(z, 5)*sin(w*log(-z + 1)/3)/(-z + 1) + 10*pow(w, 2)*pow(z, 4)*sin(w*log(-z + 1)/3) + 12*pow(w, 2)*pow(z, 3)*sin(w*log(-z + 1)/3) + 2*pow(w, 2)*pow(z, 3)*sin(w*log(-z + 1)/3)/(-z + 1) + 6*pow(w, 2)*pow(z, 2)*sin(w*log(-z + 1)/3) + pow(w, 2)*pow(z, 2)*sin(w*log(-z + 1)/3)/(-z + 1) - 16*pow(w, 2)*z*sin(w*log(-z + 1)/3) - 36*w*pow(z, 5)*cos(w*log(-z + 1)/3) - 15*w*pow(z, 4)*cos(w*log(-z + 1)/3) + 18*w*pow(z, 2)*cos(w*log(-z + 1)/3) + 6*w*z*cos(w*log(-z + 1)/3))/(9*pow(z, 2)*(pow(z, 6) - 2*pow(z, 3) + 1)) - 2*(-18*pow(d0_1704, 2)*d6_1710*pow(z, 3) + 18*pow(d0_1704, 2)*d6_1710 - 18*pow(d0_1704, 2)*pow(z, 3)*sin(w*log(-z + 1)/3) + 18*pow(d0_1704, 2)*sin(w*log(-z + 1)/3) - 9*d6_1710*pow(w, 2)*pow(z, 2) - 27*d7_1711*pow(z, 7) + 27*d7_1711*pow(z, 4) + pow(w, 2)*pow(z, 6)*sin(w*log(-z + 1)/3) + 2*pow(w, 2)*pow(z, 5)*sin(w*log(-z + 1)/3) + 3*pow(w, 2)*pow(z, 4)*sin(w*log(-z + 1)/3) + 2*pow(w, 2)*pow(z, 3)*sin(w*log(-z + 1)/3) - 8*pow(w, 2)*pow(z, 2)*sin(w*log(-z + 1)/3) - 6*w*pow(z, 6)*cos(w*log(-z + 1)/3) - 3*w*pow(z, 5)*cos(w*log(-z + 1)/3) + 6*w*pow(z, 3)*cos(w*log(-z + 1)/3) + 3*w*pow(z, 2)*cos(w*log(-z + 1)/3))/(9*pow(z, 3)*(pow(z, 6) - 2*pow(z, 3) + 1));

}

double j1_0(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return (-pow(d2_1706, 2)*pow(z, 2) + 2*pow(z, 3) - 2)/(pow(z, 2)*(pow(z, 6) - 2*pow(z, 3) + 1));

}

double j1_1(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return (-pow(z, 7) - pow(z, 4) + 2*z)/(pow(z, 2)*(pow(z, 6) - 2*pow(z, 3) + 1));

}

double j1_2(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return -2*d0_1704*d2_1706/(pow(z, 6) - 2*pow(z, 3) + 1);

}

double j1_3(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return 0;

}

double j1_4(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return 0;

}

double j1_5(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return 0;

}

double j1_6(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return 0;

}

double j1_7(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return 0;

}

double j3_0(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return -4*d0_1704*d2_1706/(pow(z, 2)*(pow(z, 3) - 1));

}

double j3_1(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return 0;

}

double j3_2(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return -2*pow(d0_1704, 2)/(pow(z, 2)*(pow(z, 3) - 1));

}

double j3_3(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return 0;

}

double j3_4(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return 0;

}

double j3_5(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return 0;

}

double j3_6(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return 0;

}

double j3_7(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return 0;

}

double j5_0(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return (-36*d0_1704*d4_1708*pow(z, 3) + 36*d0_1704*d4_1708 - 36*d0_1704*pow(z, 3)*cos(w*log(-z + 1)/3) + 36*d0_1704*cos(w*log(-z + 1)/3))/(9*pow(z, 2)*(pow(z, 6) - 2*pow(z, 3) + 1));

}

double j5_1(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return 0;

}

double j5_2(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return 0;

}

double j5_3(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return 0;

}

double j5_4(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return (-18*pow(d0_1704, 2)*pow(z, 3) + 18*pow(d0_1704, 2) - 9*pow(w, 2)*pow(z, 2))/(9*pow(z, 2)*(pow(z, 6) - 2*pow(z, 3) + 1));

}

double j5_5(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return (-27*pow(z, 7) + 27*pow(z, 4))/(9*pow(z, 2)*(pow(z, 6) - 2*pow(z, 3) + 1));

}

double j5_6(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return 0;

}

double j5_7(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return 0;

}

double j7_0(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return (-36*d0_1704*d6_1710*pow(z, 3) + 36*d0_1704*d6_1710 - 36*d0_1704*pow(z, 3)*sin(w*log(-z + 1)/3) + 36*d0_1704*sin(w*log(-z + 1)/3))/(9*pow(z, 2)*(pow(z, 6) - 2*pow(z, 3) + 1));

}

double j7_1(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return 0;

}

double j7_2(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return 0;

}

double j7_3(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return 0;

}

double j7_4(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return 0;

}

double j7_5(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return 0;

}

double j7_6(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return (-18*pow(d0_1704, 2)*pow(z, 3) + 18*pow(d0_1704, 2) - 9*pow(w, 2)*pow(z, 2))/(9*pow(z, 2)*(pow(z, 6) - 2*pow(z, 3) + 1));

}

double j7_7(double z, double d0_1704, double d1_1705, double d2_1706, double d3_1707, double d4_1708, double d5_1709, double d6_1710, double d7_1711, double w) {

   return (-27*pow(z, 7) + 27*pow(z, 4))/(9*pow(z, 2)*(pow(z, 6) - 2*pow(z, 3) + 1));

}
