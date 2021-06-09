#include <assert.h>
#include <math.h>

#include "av1/common/common.h"
#include "av1/common/spherical_pred.h"

void av1_sphere_to_plane_erp(double phi, double theta, int width, int height, double* x, double* y){
    // The 0.01 is for issues when comparing floating point numbers
    // TODO(yaoyaogoogle): Find a better way to solve the issue
    assert(phi >= -PI * 0.5 - 0.01 && phi <= PI * 0.5 + 0.01 &&
        theta >= -PI - 0.01 && theta <= PI + 0.01);

    // This should actually be a range related to 1/cos(phi) since x is distorted
    // TODO(yaoyaogoogle): Adjust the width of x according to the interpolation mode
    *x = theta/PI * (double)width * 0.5 + width * 0.5;
    // No minus sign for y since we only use an imaginary upside-down globe
    *y = phi/PI * 0.5 * (double)height * 0.5 + height * 0.5;
}

void av1_plane_to_sphere_erp(double x, double y, int width, int height, double *phi, double *theta){
    assert(x < width && x >=0 && y < height && y >=0);

    *theta = (x - width * 0.5) / (double)width * 2 * PI;
    *phi = (y - height * 0.5) / (double)height * PI;
}