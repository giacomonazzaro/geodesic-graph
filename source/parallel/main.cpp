#include <stdio.h>

#include <deque>
#include <vector>

#include "../ext/yocto-gl/yocto/yocto_common.h"
#include "../ext/yocto-gl/yocto/yocto_sceneio.h"
#include "mesh_ripper.h"
using namespace yocto;

#define here() printf("here: %d\n", __LINE__);
#define LOGGING 0
#define SAVE_MESHES 0

struct timer {
  double time = 0;
  timer() { reset(); }
  void reset() { time = get_time() * 1e-9; }
  void stop() { time = get_time() * 1e-9 - time; }
  void log(const string& s) {
#if LOGGING
    float t = get_time() * 1e-9 - time;
    printf("%f: %s\n", t, s.c_str());
#endif
  }
  // ~timer() { time = get_time() * 1e-9 - time; }

  operator double() const { return time; }
};

vector<vector<vec4f>> make_colors(const vector<vector<float>>& fields) {
  auto colors = vector<vector<vec4f>>(fields.size());
  for (int i = 0; i < fields.size(); i++) {
    colors[i]  = vector<vec4f>(fields[i].size());
    auto color = vec4f{cosf(2 * i), 0.05, sinf(3 * i), 1.0f};
    for (int k = 0; k < colors[i].size(); k++) {
      colors[i][k] = color;
      if (fields[i][k] != flt_max) {
        colors[i][k] += vec4f(sinf(1500 * fields[i][k]) * 0.5 + 0.5);
        colors[i][k].w = 1;
      }
    }
  }
  return colors;
}

vector<vector<float>> geodesic_distance_field(const MeshRipper& mesh,
    const vector<geodesic_solver>&                              graphs,
    const vector<vector<int>>&                                  _sources) {
  auto fields = vector<vector<float>>(mesh.patches.size());
  for (int i = 0; i < fields.size(); i++) {
    auto& field = fields[i];
    field       = vector<float>(mesh.patches[i].positions.size(), flt_max);
  }

  auto sources = _sources;
  for (int i = 0; i < fields.size(); i++) {
    for (auto& s : sources[i]) {
      fields[i][s] = 0;
    }
  }

  // Temporary data used for synchronization after each round.
  auto sync_data = vector<vector<std::tuple<int, int, float>>>(
      mesh.patches.size());
  auto f = [&](int k) {
    auto& s = sources[k];
    if (s.empty()) return;

    auto  t     = timer();
    auto& patch = mesh.patches[k];
    update_geodesic_distances(fields[k], graphs[k], s);
    // t.log("solve ["+std::to_string(s.size())+"] "+std::to_string(k));
  };
  auto make_sync_data = [&](int k) {
    if (sources[k].empty()) return;

    for (const auto& v : mesh.patches[k].boundaries) {
      auto [vertex, other_patch, other_vertex] = v;
      if (other_patch == -1) continue;
      if (fields[k][vertex] < fields[other_patch][other_vertex]) {
        sync_data[k].push_back({other_patch, other_vertex, fields[k][vertex]});
      }
    }
  };

  int round = 0;
  while (true) {
    round += 1;
    auto round_t = timer();

    auto t = timer();
    parallel_for(0, mesh.patches.size(), f);
    t.log("solve");

#if SAVE_MESHES
    save_mesh_ripper(
        "meshes/" + std::to_string(round) + ".ply", mesh, make_colors(fields));
#endif

    // Synchronize data for next round
    {
      t.reset();
      parallel_for(0, mesh.patches.size(), make_sync_data);

      bool we_have_syncs = false;
      for (auto& s : sources) {
        s.clear();
      }
      for (auto& sync : sync_data) {
        we_have_syncs |= not sync.empty();

        for (auto& v : sync) {
          auto [patch_id, vertex, value] = v;
          fields[patch_id][vertex]       = value;
          sources[patch_id].push_back(vertex);
        }
      }
      if (not we_have_syncs) break;
      for (auto& s : sync_data) {
        s.clear();
      }
      t.log("syncing");
    }

    round_t.log("round " + std::to_string(round) + "\n");
  }
  return fields;
}

int main(int num_args, const char* args[]) {
  auto cli = make_cli("mesh-ripper", "usage");

  string filename;
  int    num_patches = 12;
  add_cli_option(cli, "model", filename, "usage", true);
  add_cli_option(cli, "--patches", num_patches, "usage");
  parse_cli(cli, num_args, args);

  auto shape = sceneio_shape{};
  load_shape(filename, shape);
  auto& positions = shape.positions;
  auto& triangles = shape.triangles;

  assert(triangles.size());
  //  std::tie(triangles, positions) = subdivide_triangles(triangles, positions,
  //  1);
  printf("num_triangles: %ld\n", triangles.size());

  {  // Normalize in [-1, 1]^3
    bbox3f bbox;
    for (auto& p : positions) expand(bbox, p);
    vec3f center = (bbox.max + bbox.min) * 0.5f;
    float scale  = 1.0f / max(bbox.max - bbox.min);
    for (auto& p : positions) p = (p - center) * scale;
  }

  auto t    = timer();
  auto mesh = make_mesh_ripper(triangles, positions, num_patches);
  t.log("make mesh ripper");

  {  // serialial geodesic
    auto adjacencies = face_adjacencies(triangles);
    auto graph       = make_geodesic_solver(triangles, adjacencies, positions);

    auto t     = timer();
    auto field = compute_geodesic_distances(graph, {0}, flt_max);
    t.stop();
    printf("  SERIAL SOLVE: %lf\n", double(t));
  }

  {  // parallel geodesics
    t.reset();
    auto graphs      = vector<geodesic_solver>(mesh.patches.size());
    auto make_graphs = [&](int i) {
      auto& patch       = mesh.patches[i];
      auto  adjacencies = face_adjacencies(patch.triangles);
      graphs[i]         = make_geodesic_solver(
          patch.triangles, adjacencies, patch.positions);
    };
    parallel_for(0, mesh.patches.size(), make_graphs);
    t.log("make graphs");

    t.reset();
    auto sources = vector<vector<int>>(mesh.patches.size());
    sources[0]   = {0};
    sources[1]   = {0};
    auto fields  = geodesic_distance_field(mesh, graphs, sources);
    t.stop();
    printf("PARALLEL SOLVE: %lf\n", double(t));

    // save_mesh_ripper("mesh.ply", mesh, make_colors(fields));
  }
}
