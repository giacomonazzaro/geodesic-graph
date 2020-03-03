#include "../ext/yocto-gl/yocto/yocto_commonio.h"
#include "../geodesic_graph.h"
// #include "yocto-gl/yocto/yocto_math.h"
// #include "yocto-gl/yocto/yocto_shape.h"
using namespace yocto;

struct Patch {
  vector<vec3f> positions  = {};
  vector<vec3i> triangles  = {};
  vector<vec3i> boundaries = {};
};

struct MeshRipper {
  vector<Patch> patches;

  const vec3f& positions(int patch, int index) const {
    return patches[patch].positions[index];
  }
  vec3f& positions(int patch, int index) {
    return patches[patch].positions[index];
  }

  const vec3i& triangles(int patch, int index) const {
    return patches[patch].triangles[index];
  }
  vec3i& triangles(int patch, int index) {
    return patches[patch].triangles[index];
  }
};

void save_mesh_ripper(const string& filename, const MeshRipper& mesh,
                      const vector<vector<vec4f>>& patch_colors = {});
// {
//  assert(mesh.patches.size() == patch_colors.size());
//  auto triangles = vector<vec3i>{};
//  auto positions = vector<vec3f>{};
//  auto colors    = vector<vec4f>{};
//  for (int i = 0; i < mesh.patches.size(); i++) {
//    auto& patch = mesh.patches[i];
//    for (auto& t : patch.triangles) {
//      triangles.push_back(t + positions.size());
//    }
//    positions += patch.positions;
//
//    if (patch_colors.size()) colors += patch_colors[i];
//  }
//  save_shape(filename, {}, {}, triangles, {}, positions, {}, {}, colors, {});
//}

template <typename Update, typename Skip>
void visit_triangles(vector<int>& field, const vector<vec3i>& adjacency,
    const int source, Update&& update, Skip&& skip) {
  auto in_queue = vector<bool>(adjacency.size(), false);

  // setup queue
  auto queue       = std::deque<int>();
  in_queue[source] = true;
  queue.push_back(source);
  field[source] = 0;

  // Cumulative weights of elements in queue. Used to keep track of the
  // average weight of the queue.
  double cumulative_weight = 0.0;

  while (!queue.empty()) {
    auto node           = queue.front();
    auto average_weight = (float)cumulative_weight / queue.size();

    // Large Label Last (see comment at the beginning)
    for (auto tries = 0; tries < queue.size() + 1; tries++) {
      if (field[node] <= average_weight) break;
      queue.pop_front();
      queue.push_back(node);
      node = queue.front();
    }

    // Remove node from queue.
    queue.pop_front();
    in_queue[node] = false;
    cumulative_weight -= field[node];

    // Check early exit condition.
    if (skip(node)) continue;

    for (int i = 0; i < 3; i++) {
      // Distance of neighbor through this node
      auto new_distance = field[node] + 1;
      auto neighbor     = adjacency[node][i];
      if (neighbor == -1) continue;

      auto old_distance = field[neighbor];
      if (new_distance >= old_distance) continue;

      if (in_queue[neighbor]) {
        // If neighbor already in queue, don't add it.
        // Just update cumulative weight.
        cumulative_weight += new_distance - old_distance;
      } else {
        // If neighbor not in queue, add node to queue using Small Label
        // First (see comment at the beginning).
        if (queue.empty() || (new_distance < field[queue.front()]))
          queue.push_front(neighbor);
        else
          queue.push_back(neighbor);

        // Update queue information.
        in_queue[neighbor] = true;
        cumulative_weight += new_distance;
      }

      // Update distance of neighbor.
      field[neighbor] = new_distance;
      update(node, neighbor);
    }
  }
}

// Sampling strategy is farthest point sampling (FPS): at every step
// take the farthers triangle from current sampled set until done.
vector<int> sample_uniform_triangles(
    const vector<vec3i>& adjacency, int num_samples) {
  auto faces = vector<int>{};
  faces.reserve(num_samples);
  auto distances = vector<int>(adjacency.size(), int_max);
  while (true) {
    auto max_index =
        (int)(std::max_element(distances.begin(), distances.end()) -
              distances.begin());
    faces.push_back(max_index);
    if (faces.size() >= num_samples) break;
    distances[max_index] = 0.0f;
    // update_geodesic_distances(distances, adjacency, {max_index}, flt_max);
    auto update = [](int node, int neighbor) {};
    auto skip   = [](int node) { return false; };
    visit_triangles(distances, adjacency, max_index, update, skip);
  }
  return faces;
}

MeshRipper make_mesh_ripper(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, int num_patches) {
  auto mesh = MeshRipper{};

  // Init empty patches.
  mesh.patches = vector<Patch>(num_patches);

  // Make triangle adjacency.
  auto adjacency = face_adjacencies(triangles);

  // Tag faces to create uniform patches
  auto sources = sample_uniform_triangles(adjacency, num_patches);
  auto tags    = vector<int>(triangles.size(), -1);
  auto field   = vector<int>(triangles.size(), int_max);

  for (int i = 0; i < sources.size(); i++) {
    int  tag    = i;
    auto update = [&](int node, int neighbor) { tags[neighbor] = tag; };
    auto skip   = [&](int node) { return false; };
    visit_triangles(field, adjacency, sources[i], update, skip);
    tags[sources[i]] = tag;
  }

  // Map for each vertex its index in the patches (boundary
  // vertices can be in more than one patch).
  auto indices = vector<vector<int>>(num_patches);
  for (auto& i : indices) i = vector<int>(positions.size(), -1);

  // Map for each face its index in the patch (we know the patch
  // thanks to tags).
  auto face_map = vector<int>(triangles.size());

  // Fill patches with positions and triangles.
  for (int i = 0; i < triangles.size(); ++i) {
    int tag = tags[i];
    if (tag == -1) {
      // If a face was not tagged, something went wrong before.
      printf("face %d has tag -1\n", i);
      continue;
    }
    assert(tag != -1);

    // Add face to patch
    face_map[i] = mesh.patches[tag].triangles.size();
    mesh.patches[tag].triangles.push_back(triangles[i]);

    auto& tr = mesh.patches[tag].triangles.back();

    // For each vertex in the face.
    for (int k = 0; k < 3; ++k) {
      // Add vertex to patch.
      auto  mesh_index  = tr[k];
      auto& patch_index = indices[tag][mesh_index];
      if (patch_index == -1) {
        patch_index = mesh.patches[tag].positions.size();
        mesh.patches[tag].positions.push_back(positions[mesh_index]);
      }
      // Update triangle with new index.
      tr[k] = patch_index;
    }
  }

  // Create boundary data.
  for (int i = 0; i < triangles.size(); i++) {
    int tag = tags[i];
    if (tag == -1) continue;

    for (int k = 0; k < 3; k++) {
      int neighbor = adjacency[i][k];

      if (neighbor == -1) {
        int face        = face_map[i];
        int this_vertex = mesh.patches[tag].triangles[face][k];
        mesh.patches[tag].boundaries.push_back({this_vertex, -1, -1});
        continue;
      }

      int other_tag = tags[neighbor];
      if (other_tag == tag) continue;

      int face         = face_map[i];
      int this_vertex  = mesh.triangles(tag, face)[k];
      int other_vertex = -1;
      int other_face   = face_map[neighbor];
      for (int kk = 0; kk < 3; kk++) {
        // if(adjacency[neighbor][kk] == i) {
        //     // other_vertex =
        //     mesh.patches[other_tag].triangles[other_face][kk]; other_vertex =
        //     indices[other_tag][triangles[neighbor][kk]]; break;
        // }
        int  o = mesh.triangles(other_tag, other_face)[kk];
        auto a = mesh.positions(tag, this_vertex);
        auto b = mesh.positions(other_tag, o);
        // auto a = mesh.patches[tag].positions[this_vertex];
        // auto b = mesh.patches[other_tag].positions[other_vertex];
        if (a == b) {
          other_vertex = o;
          break;
        }
      }
      if (other_vertex == -1) {
        printf("%d [%d %d] [%d %d]\n", i, tag, face, other_tag, other_face);
        continue;
      }
      assert(other_vertex != -1);
      {
        auto a = mesh.patches[tag].positions[this_vertex];
        auto b = mesh.patches[other_tag].positions[other_vertex];
        if (a != b) printf("a: %d, b: %d\n", this_vertex, other_vertex);
        assert(a == b);
      }

      mesh.patches[tags[i]].boundaries.push_back(
          {this_vertex, tags[neighbor], other_vertex});
      // printf("(%d) %d %d %d\n", tags[i],
      // mesh.patches[tag].triangles[face][k], tags[neighbor], other_vertex);
    }
  }

  return mesh;
}
