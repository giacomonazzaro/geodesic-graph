#include <deque>

#include "ext/yocto-gl/yocto/yocto_math.h"
using namespace yocto;

// Data structure used for geodesic computation
struct geodesic_solver {
  static const int min_arcs = 12;
  struct graph_edge {
    int   node   = -1;
    float length = flt_max;
  };
#ifdef YOCTO_ABSEIL
  vector<short_vector<adjancency_list, min_arcs>> graph = {};
#else
  vector<vector<graph_edge>> graph = {};
#endif
};

// Build adjacencies between faces (sorted counter-clockwise)
vector<vec3i> face_adjacencies(const vector<vec3i>& triangles);

// Construct a a graph to compute geodesic distances
geodesic_solver make_geodesic_solver(const vector<vec3i>& triangles,
    const vector<vec3i>& adjacencies, const vector<vec3f>& positions);

// Compute geodesic distances
void update_geodesic_distances(vector<float>& distances,
    const geodesic_solver& solver, const vector<int>& sources,
    float max_distance = flt_max);

vector<float> compute_geodesic_distances(const geodesic_solver& solver,
    const vector<int>& sources, float max_distance = flt_max);

// Compute all shortest paths from source vertices to any other vertex.
// Paths are implicitly represented: each node is assignes its previous node in
// the path. Graph search early exits when reching end_vertex.
vector<int> compute_geodesic_paths(const geodesic_solver& solver,
    const vector<int>& sources, int end_vertex = -1);

// Sample vertices with a Poisson distribution using geodesic distances.
// Sampling strategy is farthest point sampling (FPS): at every step
// take the farthers point from current sampled set until done.
vector<int> sample_vertices_poisson(
    const geodesic_solver& solver, int num_samples);

// Compute the distance field needed to compute a voronoi diagram
vector<vector<float>> compute_voronoi_fields(
    const geodesic_solver& solver, const vector<int>& generators);

// Convert distances to colors
vector<vec4f> colors_from_field(const vector<float>& field, float scale = 1,
    const vec4f& c0 = {1, 1, 1, 1}, const vec4f& c1 = {1, 0.1, 0.1, 1});
