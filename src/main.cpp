#include "geometry_images.h"
#include "io_utils.h"
#include <igl/writeOBJ.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>
#include <igl/facet_components.h>
#include <igl/jet.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <string>
#include <iostream>

int main(int argc, char *argv[])
{
    // Load input mesh
    Eigen::MatrixXd V, UV, UV_geom;
    Eigen::MatrixXi F;
    const char* mesh_filepath = (argc > 1) ? argv[1] : "../data/beetle.obj";
    bool use_igl_reader = true;
    if (use_igl_reader)
        igl::read_triangle_mesh(mesh_filepath, V, F);
    else
        io_utils::read_triangle_mesh(mesh_filepath, V, F);

    // Compute parameterization
    Eigen::MatrixXd Vcut;
    Eigen::MatrixXi Fcut;
    geometry_image(V, F, UV_geom);

    // Fit parameterization in unit sphere
    bool do_normalize = false;
    if (do_normalize)
    {
        const auto normalize = [](Eigen::MatrixXd & UV)
        {
            UV.rowwise() -= UV.colwise().mean().eval();
            UV.array() /= (UV.colwise().maxCoeff() - UV.colwise().minCoeff()).maxCoeff() / 2.0;
        };
        normalize(V);
        normalize(UV_geom);
    }

    // Write parameterization
    std::vector<float> uv_param(UV_geom.size());
    int num_vertices = V.rows();
    for (int vi = 0; vi < num_vertices; vi++)
    {
        uv_param[0 * num_vertices + vi] = UV_geom(vi, 0);
        uv_param[1 * num_vertices + vi] = UV_geom(vi, 1);
    }
    std::string param_outpath = std::string("uv_param.pfm");
    io_utils::write_pfm(param_outpath.c_str(), uv_param.data(), num_vertices, 2, 1);
    printf("Wrote per-vertex parameterization to %s.\n", param_outpath.c_str());

    // Write geometry images
    Eigen::MatrixXd N;
    igl::per_vertex_normals(V, F, N);
    int output_res = 256;
    std::vector<float> vertex_image(output_res * output_res * 3, 0.f);
    std::vector<float> normal_image(output_res * output_res * 3, 0.f);
    for (int vi = 0; vi < num_vertices; vi++)
    {
        // Forward mapping
        double u = (UV_geom(vi, 0) + 1.0) / 2.0;  // [0, 1]
        double v = (UV_geom(vi, 1) + 1.0) / 2.0;  // [0, 1]
        int x = static_cast<int>(u * (output_res - 1));
        int y = static_cast<int>(v * (output_res - 1));
        if (x >= 0 && x < output_res && y >= 0 && y < output_res)
        {
            vertex_image[(y * output_res + x) * 3 + 0] = V(vi, 0);
            vertex_image[(y * output_res + x) * 3 + 1] = V(vi, 1);
            vertex_image[(y * output_res + x) * 3 + 2] = V(vi, 2);
            normal_image[(y * output_res + x) * 3 + 0] = N(vi, 0);
            normal_image[(y * output_res + x) * 3 + 1] = N(vi, 1);
            normal_image[(y * output_res + x) * 3 + 2] = N(vi, 2);
        }
    }
    std::string vertices_outpath = std::string("vertices.pfm");
    std::string normals_outpath  = std::string("normals.pfm");
    io_utils::write_pfm(vertices_outpath.c_str(), vertex_image.data(), output_res, output_res, 3);
    io_utils::write_pfm(normals_outpath.c_str(),  normal_image.data(), output_res, output_res, 3);
    printf("Wrote vertex geometry image to %s.\n", vertices_outpath.c_str());
    printf("Wrote normal geometry image to %s.\n", normals_outpath.c_str());

    // Set up viewer
    printf("\n");
    igl::opengl::glfw::Viewer viewer;
    std::cout << R"(
  [space] Toggle whether displaying 3D surface or 2D parameterization
  C,c     Toggle checkerboard
)";
    bool plot_parameterization = false;
    bool visualize_facet_components = true;
    const auto & update = [&]()
    {
        if (plot_parameterization)
        {
            // Viewer wants 3D coordinates, so pad UVs with column of zeros
            viewer.data().set_vertices(
                (Eigen::MatrixXd(V.rows(), 3) <<
                    UV.col(0), Eigen::VectorXd::Zero(V.rows()), UV.col(1)).finished());
        }
        else
        {
            viewer.data().set_vertices(V);
        }
        if (!visualize_facet_components)
        {
            viewer.data().compute_normals();
            viewer.data().set_uv(UV * 10);
        }
    };
    viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer &, unsigned int key, int)
    {
        switch (key)
        {
            case ' ':
                plot_parameterization ^= 1;
                break;
            case 'l':
                UV = UV_geom;
                break;
            case 'C':
            case 'c':
                viewer.data().show_texture ^= 1;
                break;
            default:
                return false;
        }
        update();
        return true;
    };

    UV = UV_geom;
    viewer.data().set_mesh(V, F);
    if (visualize_facet_components)
    {
        Eigen::VectorXd components;
        igl::facet_components(F, components);

        Eigen::MatrixXd component_colors;
        igl::jet(components, true, component_colors);
        viewer.data().set_colors(component_colors);

        // Trigger per-face shading
        igl::per_face_normals(V, F, N);
        viewer.data().set_normals(N);
    }
    else
    {
        viewer.data().set_colors(N.array() * 0.5 + 0.5);
    }
    update();
    viewer.data().show_texture = false;
    viewer.data().show_lines = false;
    viewer.launch();

    return EXIT_SUCCESS;
}
