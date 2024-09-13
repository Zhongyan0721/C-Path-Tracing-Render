#include "spdlog/spdlog.h"
#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc
// Optional. define TINYOBJLOADER_USE_MAPBOX_EARCUT gives robust trinagulation. Requires C++11
// #define TINYOBJLOADER_USE_MAPBOX_EARCUT
#include "tiny_obj_loader.h"
#include "common.h"
#include "triangle.h"
#include <iostream>
namespace muni
{
  Vec3f cal_face_normal(Vec3f v0, Vec3f v1, Vec3f v2)
  {
    Vec3f e1 = v1 - v0;
    Vec3f e2 = v2 - v0;
    return normalize(cross(e1, e2));
  }
  std::vector<Triangle> load_obj(std::string inputfile, int material_id = 7)
  {

    tinyobj::ObjReaderConfig reader_config;
    std::vector<Triangle> ret;
    tinyobj::ObjReader reader;

    if (!reader.ParseFromFile(inputfile, reader_config))
    {
      if (!reader.Error().empty())
      {
       spdlog::error("TinyObjReader: {}", reader.Warning().c_str());
      }
      exit(1);
    }

    if (!reader.Warning().empty())
    {
      spdlog::warn("TinyObjReader: {}", reader.Warning().c_str());
    }

    auto &attrib = reader.GetAttrib();
    auto &shapes = reader.GetShapes();

    // Loop over shapes
    for (size_t s = 0; s < shapes.size(); s++)
    {
      // Loop over faces(polygon)
      size_t index_offset = 0;
      for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++)
      {
        size_t fv = size_t(shapes[s].mesh.num_face_vertices[f]);
        Triangle tri;
        tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + 0];
        tri.v0 = Vec3f(attrib.vertices[3 * size_t(idx.vertex_index) + 0],
                       attrib.vertices[3 * size_t(idx.vertex_index) + 1],
                       attrib.vertices[3 * size_t(idx.vertex_index) + 2]);
        idx = shapes[s].mesh.indices[index_offset + 1];
        tri.v1 = Vec3f(attrib.vertices[3 * size_t(idx.vertex_index) + 0],
                       attrib.vertices[3 * size_t(idx.vertex_index) + 1],
                       attrib.vertices[3 * size_t(idx.vertex_index) + 2]);

        idx = shapes[s].mesh.indices[index_offset + 2];
        tri.v2 = Vec3f(attrib.vertices[3 * size_t(idx.vertex_index) + 0],
                        attrib.vertices[3 * size_t(idx.vertex_index) + 1],
                        attrib.vertices[3 * size_t(idx.vertex_index) + 2]);
        tri.face_normal = cal_face_normal(tri.v0, tri.v1, tri.v2);
        tri.emission = Vec3f(0, 0, 0);
        tri.material_id = material_id;
        ret.push_back(tri);
        index_offset += fv;
        // per-face material
        shapes[s].mesh.material_ids[f];
      }
    }
    return ret;
  }
}