#include "ext/yocto-gl/yocto/yocto_commonio.h"
#include "ext/yocto-gl/yocto/yocto_sceneio.h"
#include "geodesic_graph.h"
using namespace yocto;

inline float get_time() {
  auto t = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  return float(t) * 1e-9;
}

int main(int num_args, const char* args[]) {
  string filename;
  int    source = 0;
  auto   cli    = make_cli("geodesic-graph", "");
  add_cli_option(cli, "model", filename, "Model filename", true);
  add_cli_option(cli, "source", source, "Vertex source index", false);
  parse_cli(cli, num_args, args);

  auto shape = sceneio_shape{};
  load_shape(filename, shape);

  auto time        = get_time();
  auto adjacencies = face_adjacencies(shape.triangles);
  auto solver      = make_geodesic_solver(
      shape.triangles, adjacencies, shape.positions);
  float build_time = get_time() - time;

  time             = get_time();
  auto  field      = compute_geodesic_distances(solver, {source});
  float solve_time = get_time() - time;

  printf("\n");
  printf("  Num triangles: %ld\n", shape.triangles.size());
  printf("  Num vertices:  %ld\n", field.size());
  printf("    --- \n");
  printf("  Graph build:    %fs\n", build_time);
  printf("  Geodesic solve: %fs\n", solve_time);
}
