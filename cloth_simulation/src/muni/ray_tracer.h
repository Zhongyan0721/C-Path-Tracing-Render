#pragma once
#include "common.h"
#include "ray_tracer.h"
#include "triangle.h"
#include "math_helpers.h"
#include <cstdint>
#include <memory>
#include <numeric>
#include <tuple>

namespace muni { namespace RayTracer {

struct BoundingBox3f {
    BoundingBox3f &include(const Vec3f &point) {
        min_point = linalg::min(min_point, point);
        max_point = linalg::max(max_point, point);
        return *this;
    }

    Vec3f get_center() const { return 0.5f * (min_point + max_point); }

    Vec3f get_corner(int index) const {
        Vec3f result;
        for (int i = 0; i < 3; ++i)
            result[i] = (index & (1 << i)) ? max_point[i] : min_point[i];
        return result;
    }

    bool bounds_overlap_triangle(const Triangle &tri) const {
        Vec3f tri_min = linalg::min(tri.v0, linalg::min(tri.v1, tri.v2));
        Vec3f tri_max = linalg::max(tri.v0, linalg::max(tri.v1, tri.v2));
        return (tri_min.x <= max_point.x && tri_max.x >= min_point.x &&
                tri_min.y <= max_point.y && tri_max.y >= min_point.y &&
                tri_min.z <= max_point.z && tri_max.z >= min_point.z);
    }

    std::tuple<bool, float, float> ray_intersect(const Vec3f &ray_pos,
                                                 const Vec3f &ray_dir) const {
        const Vec3f inv_dir = 1.f / ray_dir;
        const Vec3f lo = (min_point - ray_pos) * inv_dir;
        const Vec3f hi = (max_point - ray_pos) * inv_dir;
        const Vec3f tmin = min(lo, hi), tmax = max(lo, hi);
        const float t_near =
            std::fmax(0.f, std::fmax(tmin[0], std::fmax(tmin[1], tmin[2])));
        const float t_far = std::fmin(tmax[0], std::fmin(tmax[1], tmax[2]));
        return {t_near <= t_far, t_near, t_far};
    }

    Vec3f min_point;
    Vec3f max_point;
};

struct OctreeNode {
    OctreeNode(const BoundingBox3f &bounds,
               const std::vector<uint32_t> &triangle_indices,
               const bool is_leaf = false)
        : bounds(bounds), triangle_indices(triangle_indices), is_leaf(is_leaf) {
        for (int i = 0; i < 8; i++) children[i] = nullptr;
    }

    BoundingBox3f bounds;
    std::vector<uint32_t> triangle_indices;
    std::array<std::unique_ptr<OctreeNode>, 8> children;
    bool is_leaf;
};

struct Octree {
    std::tuple<bool, float, Triangle>
    basic_octree_traversal(const std::vector<Triangle> &triangles,
                           const OctreeNode &node, Vec3f ray_pos, Vec3f ray_dir,
                           const float t_max, const bool shadow_ray) const {
        // Check if they ray intersects the bounding box of the node
        auto [hit, t_near, t_far] = node.bounds.ray_intersect(ray_pos, ray_dir);
        if (!hit || t_far < 0 || t_near > t_max) return {false, 0, Triangle()};

        // If the node is a leaf, check the triangles
        if (node.is_leaf) {
            float t_min = std::numeric_limits<float>::infinity();
            Triangle nearest_tri;
            bool hit_one_tri = false;
            for (uint32_t tri_idx : node.triangle_indices) {
                const Triangle &tri = triangles[tri_idx];
                auto [hit, t] = Triangle::ray_triangle_intersect(
                    tri, ray_pos, ray_dir, EPS, t_max - ANYHIT_EPS);
                if (hit && shadow_ray) return {true, t, tri};
                if (hit && t < t_min) {
                    t_min = t;
                    nearest_tri = tri;
                    hit_one_tri = true;
                }
            }
            return {hit_one_tri, t_min, nearest_tri};
        }

        // Otherwise, recursively check the children
        float t_min = std::numeric_limits<float>::infinity();
        Triangle nearest_tri;
        bool hit_one_tri = false;
        for (int i = 0; i < 8; i++) {
            if (node.children[i] == nullptr) continue;
            auto [hit, t, tri] =
                basic_octree_traversal(triangles, *node.children[i], ray_pos,
                                       ray_dir, t_max, shadow_ray);
            if (hit && shadow_ray && t < t_max - ANYHIT_EPS)
                return {true, t, tri};
            if (hit && t < t_min) {
                t_min = t;
                nearest_tri = tri;
                hit_one_tri = true;
            }
        }
        return {hit_one_tri, t_min, nearest_tri};
    }

    std::unique_ptr<OctreeNode>
    build(const BoundingBox3f &bounds, const std::vector<Triangle> &triangles,
          const std::vector<uint32_t> &triangle_indices, int depth) {
        // If there is no triangle
        if (triangle_indices.empty()) return nullptr;

        // If there are too few triangles, return a leaf node
        if (triangle_indices.size() <= 16 || depth > 4) {
            num_leaf_node++;
            num_total_leaf_triangles += triangle_indices.size();
            return std::make_unique<OctreeNode>(bounds, triangle_indices, true);
        }

        // Otherwise, split the bounding box into 8 sub-boxes
        std::array<BoundingBox3f, 8> sub_boxes;
        for (int i = 0; i < 8; i++) {
            Vec3f min_point = Vec3f{
                std::min(bounds.get_center()[0], bounds.get_corner(i)[0]),
                std::min(bounds.get_center()[1], bounds.get_corner(i)[1]),
                std::min(bounds.get_center()[2], bounds.get_corner(i)[2])};
            Vec3f max_point = Vec3f{
                std::max(bounds.get_center()[0], bounds.get_corner(i)[0]),
                std::max(bounds.get_center()[1], bounds.get_corner(i)[1]),
                std::max(bounds.get_center()[2], bounds.get_corner(i)[2])};
            sub_boxes[i] = BoundingBox3f{min_point, max_point};
        }

        // Assign triangles to sub-boxes
        std::array<std::vector<uint32_t>, 8> sub_triangle_indices;
        for (uint32_t tri_idx : triangle_indices) {
            const Triangle &tri = triangles[tri_idx];
            for (int i = 0; i < 8; i++) {
                if (sub_boxes[i].bounds_overlap_triangle(tri)) {
                    sub_triangle_indices[i].push_back(tri_idx);
                }
            }
        }

        // Recursively build the children
        std::unique_ptr<OctreeNode> node =
            std::make_unique<OctreeNode>(bounds, triangle_indices);
        for (int i = 0; i < 8; i++) {
            node->children[i] = build(sub_boxes[i], triangles,
                                      sub_triangle_indices[i], depth + 1);
        }

        num_interior_node++;
        return node;
    }
    Octree() {
        
    }
    Octree(const std::vector<Triangle> &triangles) {
        build_octree(triangles);
    }
    void build_octree(const std::vector<Triangle> &triangles) {
      // Compute the bounding box of the scene
        BoundingBox3f bbox = BoundingBox3f();
        for (const Triangle &tri : triangles) {
            bbox.include(tri.v0).include(tri.v1).include(tri.v2);
        }

        // Fill with triange index: 0, 1, 2, ..., n
        std::vector<uint32_t> triangle_indices(triangles.size());
        std::iota(triangle_indices.begin(), triangle_indices.end(), 0);

        // Start to build the octree
        head = build(bbox, triangles, triangle_indices, 1);
    }
    std::unique_ptr<OctreeNode> head;
    uint32_t num_leaf_node = 0;
    uint32_t num_interior_node = 0;
    uint32_t num_total_leaf_triangles = 0;
};

/** Find the closest intersection of a ray with a list of triangles.
    \param[in] tri The triangle to intersect with.
    \param[in] ray_origin The origin of the ray.
    \param[in] ray_direction The direction of the ray.
    \return A tuple containing a boolean indicating whether the ray intersects
    the triangle, the t value of the intersection point along the ray, and the
    triangle that was hit.
*/
static std::tuple<bool, float, Triangle>
closest_hit(Vec3f ray_pos, Vec3f ray_dir, const Octree &octree,
            const std::vector<Triangle> &triangles) {
    // float t_min = std::numeric_limits<float>::infinity();
    // Triangle nearest_tri;
    // bool hit_one_tri = false;
    // for (const Triangle &tri : triangles) {
    //     auto [hit, t] = Triangle::ray_triangle_intersect(
    //         tri, ray_pos, ray_dir, EPS, std::numeric_limits<float>::infinity());
    //     if (t < EPS) hit = false;
    //     if (hit && t < t_min) {
    //         t_min = t;
    //         nearest_tri = tri;
    //         hit_one_tri = true;
    //     }
    // }
    // return {hit_one_tri, t_min, nearest_tri};
    return octree.basic_octree_traversal(
        triangles, *octree.head, ray_pos, ray_dir,
        std::numeric_limits<float>::infinity(), false);
}

/** Check if a ray intersects any triangle in a list.
    \param[in] tri The triangle to intersect with.
    \param[in] ray_origin The origin of the ray.
    \param[in] ray_direction The direction of the ray.
    \param[in] t_max The maximum t value to consider.
    \return True if the ray intersects any triangle, false otherwise.
*/
static bool any_hit(Vec3f ray_pos, Vec3f ray_dir, float t_max,
                    const Octree &octree,
                    const std::vector<Triangle> &triangles) {
    // for (const Triangle &tri : triangles) {
    //     auto [hit, t] = Triangle::ray_triangle_intersect(
    //         tri, ray_pos, ray_dir, EPS, t_max - ANYHIT_EPS);
    //     if (t < EPS) hit = false;
    //     if (hit && t < t_max - ANYHIT_EPS) { return true; }
    // }
    // return false;
    auto [hit, t, tri] = octree.basic_octree_traversal(
        triangles, *octree.head, ray_pos, ray_dir, t_max, true);
    return hit;
}

}}  // namespace muni::RayTracer
