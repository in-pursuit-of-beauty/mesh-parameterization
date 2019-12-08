#include "parameterize_mesh.h"
#include <igl/read_triangle_mesh.h>
#include <igl/writeOBJ.h>
#include <igl/per_vertex_normals.h>
#include <Eigen/Core>
#include <string>
#include <iostream>

int main(int argc, char *argv[])
{
    // Load input mesh
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    const char* mesh_filepath = (argc > 1) ? argv[1] : "../data/beetle.obj";
    igl::read_triangle_mesh(mesh_filepath, V, F);

    // Compute texcoords as parameterization
    Eigen::MatrixXd TC;
    parameterize_mesh(V, F, TC);

    // Normals
    Eigen::MatrixXd VN;
    igl::per_vertex_normals(V, F, VN);  // #V x 3

    // Write new mesh with texcoords
    const char* out_filepath = (argc > 2) ? argv[2] : "out.obj";
    if (igl::writeOBJ(std::string(out_filepath), V, F, VN, F, TC, F))
        printf("Wrote mesh with texcoords to %s.\n", out_filepath);
    else
        printf("[-] ERROR: Failed to write mesh.\n");

    return EXIT_SUCCESS;
}
