#ifndef GEOMETRY_IMAGES_H
#define GEOMETRY_IMAGES_H

#include <Eigen/Core>

// Inputs:
//   V  #V by 3 list of mesh vertex positions
//   F  #F by 3 list of triangle indices into V
// Outputs:
//   U  #U by 2 list of mesh UV parameterization coordinates
//
void geometry_image(Eigen::MatrixXd & V,
                    Eigen::MatrixXi & F,
                    Eigen::MatrixXd & U);

#endif
