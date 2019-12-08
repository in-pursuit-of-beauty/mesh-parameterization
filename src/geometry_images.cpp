#include "geometry_images.h"
#include <deque>
#include <igl/slice.h>
#include <igl/cut_to_disk.h>
#include <igl/is_edge_manifold.h>
#include <igl/extract_manifold_patches.h>
#include <igl/remove_unreferenced.h>
#include <igl/remove_duplicates.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/edge_topology.h>
#include <igl/is_boundary_edge.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/is_border_vertex.h>
#include <igl/harmonic.h>
#include <igl/cut_mesh.h>
#include <igl/boundary_loop.h>
#include <igl/facet_components.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/exact_geodesic.h>
// #include <igl/copyleft/cgal/remesh_self_intersections.h>
// #include <igl/copyleft/cgal/outer_hull.h>
#include <cmath>
#include <unordered_set>

int check_components(const Eigen::MatrixXd & V,
                     const Eigen::MatrixXi & F,
                     int & max_component_id)
{
    Eigen::VectorXi components;  // per-face component IDs
    igl::facet_components(F, components);

    std::unordered_set<int> unique_component_ids;
    for (int ci = 0; ci < components.size(); ci++)
    {
        unique_component_ids.insert(components[ci]);
        max_component_id = std::max(components[ci], max_component_id);
    }

    return unique_component_ids.size();
}

void clean_mesh(Eigen::MatrixXd & V,
                Eigen::MatrixXi & F)
{
    // // Remesh self-intersections
    // Eigen::MatrixXd oldV = V;
    // Eigen::MatrixXi oldF = F;
    // Eigen::MatrixXi IF;
    // Eigen::VectorXi J;
    // Eigen::VectorXi IM;
    // igl::copyleft::cgal::remesh_self_intersections(
    //     oldV, oldF, {false, false, false}, V, F, IF, J, IM);
    // std::for_each(F.data(), F.data() + F.size(), [&IM](int & a) { a = IM(a); });
    // Eigen::MatrixXd SV;
    // Eigen::MatrixXi SF;
    // Eigen::VectorXi UIM;
    // igl::remove_unreferenced(V, F, SV, SF, UIM);
    // Eigen::VectorXi flip;
    // igl::copyleft::cgal::outer_hull(SV, SF, V, F, J, flip);

    // Remove unreferenced
    Eigen::MatrixXd oldV = V;
    Eigen::MatrixXi oldF = F;
    Eigen::VectorXi I;
    igl::remove_unreferenced(oldV, oldF, V, F, I);

    // // Remove duplicate vertices
    // oldV = V;
    // Eigen::MatrixXd SV;
    // Eigen::MatrixXi SVI;
    // Eigen::MatrixXi SVJ;
    // igl::remove_duplicate_vertices(oldV, 0.001, V, SVI, SVJ);
    // for (int fi; fi < F.rows(); fi++)
    // {
    //     F(fi, 0) = SVJ(F(fi, 0), 0);
    //     F(fi, 1) = SVJ(F(fi, 1), 0);
    //     F(fi, 2) = SVJ(F(fi, 2), 0);
    // }
}

// Remove all but the largest connected component.
void remove_non_max_components(Eigen::MatrixXd & V,
                               Eigen::MatrixXi & F,
                               const int & max_component_id)
{
    Eigen::VectorXi components;  // per-face component IDs
    igl::facet_components(F, components);

    // Determine component sizes
    std::vector<int> component_sizes(max_component_id + 1, 0);
    for (int fi = 0; fi < components.size(); fi++)
        component_sizes[components[fi]]++;

    // Identify largest component
    int largest_component = 0;
    int largest_component_size = -1;
    for (int ci = 0; ci < component_sizes.size(); ci++)
    {
        if (component_sizes[ci] > largest_component_size)
        {
            largest_component = ci;
            largest_component_size = component_sizes[ci];
        }
    }

    // Delete non-max components
    for (int fi = F.rows() - 1; fi >= 0; fi--)
    {
        if (components[fi] != largest_component)
        {
            // Delete the offending face
            unsigned int num_faces = F.rows() - 1;
            F.block(fi, 0, num_faces - fi, 3) = \
                F.block(fi + 1, 0, num_faces - fi, 3);
            F.conservativeResize(num_faces, 3);
        }
    }

    // Clean
    clean_mesh(V, F);
}

void igl_disk_boundary(const Eigen::MatrixXd & V,
                       const Eigen::MatrixXi & F,
                       Eigen::MatrixXd & Vcut,
                       Eigen::MatrixXi & Fcut,
                       Eigen::VectorXi & boundary)
{
    std::vector<std::vector<int>> cuts;
    igl::cut_to_disk(F, cuts);

    // cuts -> fCuts: translate the vertex-style cuts into the F oriented cuts representation
    // https://github.com/libigl/libigl/issues/1328#issuecomment-546332878
    Eigen::MatrixXi fCuts;
    fCuts.resizeLike(F);
    fCuts.setZero();
    std::vector<std::vector<int>> VF, VFi;
    igl::vertex_triangle_adjacency(V, F, VF, VFi);
    for (auto i : cuts[0]) {
        for (size_t f = 0; f < VF[i].size(); f++) {
            // VF[i][f] is the face that the specific vertex is in
            // VFi[i][f] is the index of the vertex in that face
            fCuts.row(VF[i][f])[VFi[i][f]] = 1;
        }
    }

    // Edges
    Eigen::MatrixXi E;
    Eigen::MatrixXi FE;
    Eigen::MatrixXi EF;
    igl::edge_topology(V, F, E, FE, EF);

    // Face-face topology
    Eigen::MatrixXi TT;
    Eigen::MatrixXi TTi;
    igl::triangle_triangle_adjacency(F, TT, TTi);

    // Border vertex identification
    std::vector<bool> V_border = igl::is_border_vertex(V, F);

    // Open mesh
    igl::cut_mesh(V, F, VF, VFi, TT, TTi, V_border, fCuts, Vcut, Fcut);

    // Get boundary indices into Vcut
    igl::boundary_loop(Fcut, boundary);
}

// Given a 3D mesh of arbitrary genus,
// find a cut that opens the mesh into a topological disk.
void initial_cut(const Eigen::MatrixXd & V,
                 const Eigen::MatrixXi & F,
                 float seed,
                 Eigen::ArrayXi & cut_edges)
{
    // Edges
    Eigen::MatrixXi E;
    Eigen::MatrixXi FE;
    Eigen::MatrixXi EF;
    igl::edge_topology(V, F, E, FE, EF);
    printf("V: %d | F: %d | E: %d\n", V.rows(), F.rows(), E.rows());

    // Vertex-face topology
    std::vector<std::vector<int>> VF;
    std::vector<std::vector<int>> VFi;
    igl::vertex_triangle_adjacency(V, F, VF, VFi);

    // Detect the boundary edges
    Eigen::ArrayXi boundary_edges;
    igl::is_boundary_edge(E, F, boundary_edges);

    Eigen::ArrayXi remaining_edges = 1 - boundary_edges;
    Eigen::ArrayXi remaining_vertices = Eigen::ArrayXi::Ones(V.rows());
    Eigen::ArrayXi remaining_faces = Eigen::ArrayXi::Ones(F.rows());

    // Remove seed triangle
    std::deque<int> new_boundary_edges;
    int seed_idx = int(seed * ((float) F.rows()));
    remaining_faces[seed_idx] = 0;
    new_boundary_edges.push_back(FE(seed_idx, 0));
    new_boundary_edges.push_back(FE(seed_idx, 1));
    new_boundary_edges.push_back(FE(seed_idx, 2));

    std::deque<int> suspect_vertices;
    auto remove_face = [&](const int & face,
                           const int & old_boundary_edge)
    {
        remaining_faces[face] = 0;
        Eigen::Vector3i face_edges = FE.row(face);
        if (face_edges(0) != old_boundary_edge)
            new_boundary_edges.push_back(face_edges(0));
        if (face_edges(1) != old_boundary_edge)
            new_boundary_edges.push_back(face_edges(1));
        if (face_edges(2) != old_boundary_edge)
            new_boundary_edges.push_back(face_edges(2));

        for (int vi = 0; vi < 3; vi++)
        {
            int vertex = F(face, vi);
            suspect_vertices.push_back(vertex);
        }
    };

    // Iteratively remove newly-created boundary edges
    while (!new_boundary_edges.empty())
    {
        int edge = new_boundary_edges.front();
        new_boundary_edges.pop_front();  // breadth-first

        int face0 = EF(edge, 0);
        int face1 = EF(edge, 1);
        bool face0_remains = face0 != -1 && remaining_faces[face0];
        bool face1_remains = face1 != -1 && remaining_faces[face1];

        if (face0_remains && !face1_remains)
        {
            remaining_edges[edge] = 0;
            remove_face(face0, edge);
        }
        else if (!face0_remains && face1_remains)
        {
            remaining_edges[edge] = 0;
            remove_face(face1, edge);
        }
    }

    // Iteratively remove newly-dangling edges
    while (!suspect_vertices.empty())
    {
        int vertex = suspect_vertices.front();
        suspect_vertices.pop_front();  // breadth-first

        if (vertex != -1 && remaining_vertices[vertex])
        {
            // Check whether vertex is adjacent to only one edge
            int num_adjacent_edges = 0;
            int adjacent_edge = -1;
            int other_vertex = -1;
            for (auto & face : VF[vertex])
            {
                for (int ei = 0; ei < 3; ei++)
                {
                    int edge = FE(face, ei);
                    if (edge != adjacent_edge &&
                        (remaining_edges[edge] || boundary_edges[edge]))
                    {
                        if (E(edge, 0) == vertex)
                        {
                            num_adjacent_edges++;
                            adjacent_edge = edge;
                            other_vertex = E(edge, 1);
                        }
                        else if (E(edge, 1) == vertex)
                        {
                            num_adjacent_edges++;
                            adjacent_edge = edge;
                            other_vertex = E(edge, 0);
                        }
                    }
                }
            }
            if (num_adjacent_edges == 1)
            {
                remaining_vertices[vertex] = 0;
                remaining_edges[adjacent_edge] = 0;
                suspect_vertices.push_back(other_vertex);
            }
        }
    }

    if (remaining_vertices.sum() == 1)
    {
        // Add back two adjacent edges to the remaining vertex
        for (int vertex = 0; vertex < V.rows(); vertex++)
        {
            if (remaining_vertices[vertex])
            {
                bool done = false;
                int num_edges_added = 0;
                int prev_edge_added = -1;
                for (auto & face : VF[vertex])
                {
                    for (int ei = 0; ei < 3; ei++)
                    {
                        int edge = FE(face, ei);
                        if (edge != prev_edge_added &&
                            !remaining_edges[edge] &&
                            (E(edge, 0) == vertex || E(edge, 1) == vertex))
                        {
                            remaining_edges[edge] = 1;
                            num_edges_added++;
                            prev_edge_added = edge;

                            if (num_edges_added >= 2)
                                break;
                        }
                    }
                    if (done)
                        break;
                }
                break;
            }
        }
    }

    // Straighten cut-paths
    // Replace each with constrained shortest path between cut-nodes

    cut_edges = remaining_edges + boundary_edges;  // re-add boundary edges
}

void open_cut(const Eigen::MatrixXd & V,
              const Eigen::MatrixXi & F,
              const Eigen::ArrayXi & cut_edges,
              Eigen::MatrixXd & Vcut,
              Eigen::MatrixXi & Fcut,
              Eigen::VectorXi & boundary)
{
    // Edges
    Eigen::MatrixXi E;
    Eigen::MatrixXi FE;
    Eigen::MatrixXi EF;
    igl::edge_topology(V, F, E, FE, EF);

    // Vertex-face topology
    std::vector<std::vector<int>> VF;
    std::vector<std::vector<int>> VFi;
    igl::vertex_triangle_adjacency(V, F, VF, VFi);

    // Face-face topology
    Eigen::MatrixXi TT;
    Eigen::MatrixXi TTi;
    igl::triangle_triangle_adjacency(F, TT, TTi);

    // Border vertex identification
    std::vector<bool> V_border = igl::is_border_vertex(V, F);

    // Define edges to cut
    Eigen::MatrixXi cuts = Eigen::MatrixXi::Zero(F.rows(), 3);
    for (int edge = 0; edge < cut_edges.size(); edge++)
    {
        if (cut_edges[edge])
        {
            for (int fi = 0; fi < 2; fi++)
            {
                int face = EF(edge, fi);
                if (face != -1)
                {
                    int ei = 0;
                    while (FE(face, ei) != edge) ei++;
                    cuts(face, ei) = 1;
                }
            }
        }
    }

    // Open mesh
    igl::cut_mesh(V, F, VF, VFi, TT, TTi, V_border, cuts, Vcut, Fcut);

    // Get boundary indices into Vcut
    igl::boundary_loop(Fcut, boundary);
}

// Map opened cut nodes to grid points on the boundary of the unit square.
void boundary_parameterization(const Eigen::MatrixXd & V,
                               const Eigen::MatrixXi & F,
                               const Eigen::VectorXi & boundary,
                               Eigen::MatrixXd & boundary_uv)
{
    // Output N x 2 list of 2D positions on the unit square boundary
    int n = boundary.size();
    boundary_uv = Eigen::MatrixXd::Zero(n, 2);
    double increment = 5.0 / static_cast<double>(n);
    double i = 0;
    for (int vi = 0; vi < n; vi++)
    {
        switch (int(std::floor(i)))
        {
            case 0:
                boundary_uv(vi, 0) = 1.0;
                boundary_uv(vi, 1) = i;
                break;
            case 1:
                boundary_uv(vi, 0) = 1.0 - (i - 1.0) * 2.0;
                boundary_uv(vi, 1) = 1.0;
                break;
            case 2:
                boundary_uv(vi, 0) = -1.0;
                boundary_uv(vi, 1) = 1.0 - (i - 2.0) * 2.0;
                break;
            case 3:
                boundary_uv(vi, 0) = -1.0 + (i - 3.0) * 2.0;
                boundary_uv(vi, 1) = -1.0;
                break;
            case 4:
                boundary_uv(vi, 0) = 1.0;
                boundary_uv(vi, 1) = -1.0 + (i - 4.0);
                break;
        }
        i += increment;
    }
}

// Simple interior parameterization.
// ---------------------------------
// Computes the UV parameterization for each vertex as a convex
// combination of the boundary vertices' UV parameterizations,
// where weights are based on geodesic distances to boundary vertices.
// -------------------------------------------------------------------
void interior_parameterization(const Eigen::MatrixXd & V,
                               const Eigen::MatrixXi & F,
                               const Eigen::VectorXi & boundary,
                               const Eigen::MatrixXd & boundary_uv,
                               Eigen::MatrixXd & U)
{
    U.resize(V.rows(), 2);
    Eigen::VectorXd boundary_u = boundary_uv.col(0).array();
    Eigen::VectorXd boundary_v = boundary_uv.col(1).array();
    for (int vi = 0; vi < V.rows(); vi++)
    {
        Eigen::VectorXi VS, FS, FT;
        VS.resize(1);
        VS << vi;
        Eigen::VectorXd d;  // geodesic distances
        igl::exact_geodesic(V, F, VS, FS, boundary, FT, d);
        Eigen::ArrayXd weights = 1.0 / (d.array() + 0.00000001);
        d = weights.matrix().normalized();
        U(vi, 0) = d.dot(boundary_u);  // "u" in parameterization
        U(vi, 1) = d.dot(boundary_v);  // "v" in parameterization
    }
}

void geometry_image(Eigen::MatrixXd & V,
                    Eigen::MatrixXi & F,
                    Eigen::MatrixXd & U)
{
    int max_component_id = -1;
    int num_components = check_components(V, F, max_component_id);
    printf("This mesh has %d connected component(s)\n", num_components);
    if (num_components != 1)
        printf("[-] Warning: number of connected components should be 1\n");
        printf("Reducing to largest connected component...\n");
        remove_non_max_components(V, F, max_component_id);
        num_components = check_components(V, F, max_component_id);
        printf("Done. Now the mesh has %d connected component(s)\n", num_components);

    // Extract manifold patches
    Eigen::VectorXi P;
    int num_patches = igl::extract_manifold_patches(F, P);
    printf("Extracted %d manifold patches.\n", num_patches);

    // Determine manifold patch sizes
    std::vector<int> patch_sizes(num_patches, 0);
    for (int fi = 0; fi < P.rows(); fi++)
        patch_sizes[P[fi]]++;

    // Take largest manifold patch
    int largest_patch = 0;
    int largest_patch_size = -1;
    for (int pi = 0; pi < num_patches; pi++)
        if (patch_sizes[pi] > largest_patch_size)
        {
            largest_patch = pi;
            largest_patch_size = patch_sizes[pi];
        }

    // Only retain largest patch
    printf("=======================================\n");
    printf("==> Processing patch %d.\n", largest_patch);
    printf("=======================================\n");

    // Collect faces for patch
    std::vector<int> patch_face_idxs;
    for (int fi = 0; fi < P.rows(); fi++)
        if (P[fi] == largest_patch)
            patch_face_idxs.push_back(fi);
    Eigen::MatrixXd Vpatch = V.replicate(1, 1);

    Eigen::VectorXi row_idxs = \
        Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(
            patch_face_idxs.data(), patch_face_idxs.size());
    Eigen::Vector3i col_idxs(0, 1, 2);
    Eigen::MatrixXi Fpatch;
    igl::slice(F, row_idxs, col_idxs, Fpatch);

    if (!igl::is_edge_manifold(Fpatch))
        printf("[-] Warning: the patch should be edge manifold at this point\n");

    clean_mesh(Vpatch, Fpatch);

    // Convert to topological disk, define boundary
    bool single_boundary_loop = false;
    Eigen::VectorXi boundary;
    if (single_boundary_loop)
    {
        igl::boundary_loop(Fpatch, boundary);
    }
    else
    {
        bool use_igl_disk_boundary = true;
        Eigen::MatrixXd Vcut;
        Eigen::MatrixXi Fcut;
        if (use_igl_disk_boundary)
        {
            igl_disk_boundary(Vpatch, Fpatch, Vcut, Fcut, boundary);
        }
        else
        {
            Eigen::ArrayXi cut_edges;
            initial_cut(Vpatch, Fpatch, 0.59f, cut_edges);
            printf("Number of edges in initial cut: %d\n", cut_edges.sum());
            open_cut(Vpatch, Fpatch, cut_edges, Vcut, Fcut, boundary);
        }
        printf("Finished opening cut.\n");
        Vpatch = Vcut;
        Fpatch = Fcut;
    }

    // Parameterize boundary
    bool circle_boundary = false;
    Eigen::MatrixXd boundary_uv;
    if (circle_boundary)
        igl::map_vertices_to_circle(Vpatch, boundary, boundary_uv);
    else
        boundary_parameterization(Vpatch, Fpatch, boundary, boundary_uv);
    printf("Finished parameterizing boundary.\n");

    // Parameterize interior
    clean_mesh(Vpatch, Fpatch);
    interior_parameterization(Vpatch, Fpatch, boundary, boundary_uv, U);
    printf("Finished parameterizing interior.\n");
    V = Vpatch;
    F = Fpatch;
    printf("Done processing patch %d.\n", largest_patch);
}
