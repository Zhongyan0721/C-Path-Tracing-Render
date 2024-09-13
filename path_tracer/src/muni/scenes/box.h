#include "common.h"
#include "triangle.h"
#include "material.h"
#include <array>
#include <variant>

namespace muni { namespace BoxScene {

static const float light_x = 0.195f;
static const float light_y = -0.355f;
static const float light_z = 0.545f;
static const float light_len_x = 0.16f;
static const float light_len_y = 0.16f;
static const float inv_light_area = 1 / (light_len_x * light_len_y);
static const Vec3f light_color{50.0f, 50.0f, 50.0f};
static const Vec3f light_normal{0.0f, 0.0f, -1.0f};

// Microfacet materials
const Microfacet Gold{.roughness = 0.0005f, .n1 = Vec3f{1.0f}, .n2 =Vec3f{0.2177f, 0.42659f,1.2425f } };
const Microfacet Iron{.roughness = 0.02f, .n1 = Vec3f{1.0f}, .n2 =Vec3f{2.8851f, 2.95f,2.6f } };
static const std::array<std::variant<Lambertian, Microfacet>, 7> materials = {
    // Back
    Lambertian{.albedo = Vec3f{0.874000013f, 0.874000013f, 0.875000000f}},
    // Bottom
    Lambertian{.albedo = Vec3f{0.874000013f, 0.874000013f, 0.875000000f}},
    // Left
    Lambertian{.albedo = Vec3f{0.0f, 0.2117f, 0.3765f}},
    // Right
    Lambertian{.albedo = Vec3f{0.996f, 0.7373f, 0.0667f}},
    // Top
    Lambertian{.albedo = Vec3f{0.874000013f, 0.874000013f, 0.875000000f}},
    // Bunny
    Iron,
    Gold,
};

static std::vector<Triangle> triangles = {
    // Light
    Triangle{.v0 = Vec3f{light_x, light_y + light_len_y, light_z},
             .v1 = Vec3f{light_x + light_len_x, light_y, light_z},
             .v2 = Vec3f{light_x, light_y, light_z},
             .face_normal = Vec3f{0.0f, 0.0f, -1.0f},
             .emission = light_color,
             .material_id = 0},
    Triangle{.v0 = Vec3f{light_x, light_y + light_len_y, light_z},
             .v1 = Vec3f{light_x + light_len_x, light_y + light_len_y, light_z},
             .v2 = Vec3f{light_x + light_len_x, light_y, light_z},
             .face_normal = Vec3f{0.0f, 0.0f, -1.0f},
             .emission = light_color,
             .material_id = 0},
    // Back
    Triangle{.v0 = Vec3f{0.000000133f, -0.559199989f, 0.548799932f},
             .v1 = Vec3f{0.555999935f, -0.559199989f, 0.000000040f},
             .v2 = Vec3f{0.000000133f, -0.559199989f, 0.000000040f},
             .face_normal = Vec3f{0.0f, 1.0f, 0.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 0},
    Triangle{.v0 = Vec3f{0.000000133f, -0.559199989f, 0.548799932f},
             .v1 = Vec3f{0.555999935f, -0.559199989f, 0.548799932f},
             .v2 = Vec3f{0.555999935f, -0.559199989f, 0.000000040f},
             .face_normal = Vec3f{0.0f, 1.0f, 0.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 0},
    // Bottom
    Triangle{.v0 = Vec3f{0.000000133f, -0.559199989f, 0.000000040f},
             .v1 = Vec3f{0.555999935f, -0.559199989f, 0.000000040f},
             .v2 = Vec3f{0.555999935f, -0.000000119f, 0.000000040f},
             .face_normal = Vec3f{0.0f, 0.0f, 1.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 1},
    Triangle{.v0 = Vec3f{0.000000133f, -0.559199989f, 0.000000040f},
             .v1 = Vec3f{0.555999935f, -0.000000119f, 0.000000040f},
             .v2 = Vec3f{0.000000133f, -0.000000119f, 0.000000040f},
             .face_normal = Vec3f{0.0f, 0.0f, 1.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 1},
    // Left
    Triangle{.v0 = Vec3f{0.555999935f, -0.000000119f, 0.548799932f},
             .v1 = Vec3f{0.555999935f, -0.000000119f, 0.000000040f},
             .v2 = Vec3f{0.555999935f, -0.559199989f, 0.000000040f},
             .face_normal = Vec3f{-1.0f, 0.0f, 0.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 2},
    Triangle{.v0 = Vec3f{0.555999935f, -0.000000119f, 0.548799932f},
             .v1 = Vec3f{0.555999935f, -0.559199989f, 0.000000040f},
             .v2 = Vec3f{0.555999935f, -0.559199989f, 0.548799932f},
             .face_normal = Vec3f{-1.0f, 0.0f, 0.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 2},
    // Right
    Triangle{.v0 = Vec3f{0.000000133f, -0.559199989f, 0.000000040f},
             .v1 = Vec3f{0.000000133f, -0.000000119f, 0.000000040f},
             .v2 = Vec3f{0.000000133f, -0.000000119f, 0.548799932f},
             .face_normal = Vec3f{1.0f, 0.0f, 0.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 3},
    Triangle{.v0 = Vec3f{0.000000133f, -0.559199989f, 0.000000040f},
             .v1 = Vec3f{0.000000133f, -0.000000119f, 0.548799932f},
             .v2 = Vec3f{0.000000133f, -0.559199989f, 0.548799932f},
             .face_normal = Vec3f{1.0f, 0.0f, 0.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 3},
    // Top
    Triangle{.v0 = Vec3f{0.000000133f, -0.000000119f, 0.548799932f},
             .v1 = Vec3f{0.555999935f, -0.559199989f, 0.548799932f},
             .v2 = Vec3f{0.000000133f, -0.559199989f, 0.548799932f},
             .face_normal = Vec3f{0.0f, 0.0f, -1.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 4},
    Triangle{.v0 = Vec3f{0.000000133f, -0.000000119f, 0.548799932f},
             .v1 = Vec3f{0.555999935f, -0.000000119f, 0.548799932f},
             .v2 = Vec3f{0.555999935f, -0.559199989f, 0.548799932f},
             .face_normal = Vec3f{0.0f, 0.0f, -1.0f},
             .emission = Vec3f{0.0f, 0.0f, 0.0f},
             .material_id = 4},
    };
}}  // namespace muni::BoxScene
