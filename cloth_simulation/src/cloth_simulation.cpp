#include "material.h"
#include "muni/camera.h"
#include "muni/common.h"
#include "muni/image.h"
#include "muni/material.h"
#include "muni/math_helpers.h"
#include "muni/obj_loader.h"
#include "muni/ray_tracer.h"
#include "muni/sampler.h"
#include "muni/scenes/box.h"
#include "muni/triangle.h"
#include "ray_tracer.h"
#include "spdlog/spdlog.h"
#include "triangle.h"
#include <cmath>
#include <iostream>

using namespace muni;

RayTracer::Octree octree{};

Vec3f offset_ray_origin(Vec3f ray_pos, Vec3f normal) {
    return ray_pos + EPS * normal;
}

bool is_emitter(const Triangle &tri) { return tri.emission != Vec3f{0.0f}; }

Vec3f eval_area_light(const Vec3f light_dir) {
    if (dot(light_dir, BoxScene::light_normal) > 0.0f)
        return BoxScene::light_color;
    return Vec3f{0.0f};
}

std::tuple<Vec3f, Vec3f, float> sample_area_light(Vec2f samples) {
    Vec3f light_position = Vec3f(BoxScene::light_x,BoxScene::light_y,BoxScene::light_z);
    Vec3f light_normal = BoxScene::light_normal;
    Vec3f sampled_point = light_position + 
                          (samples.x) * BoxScene::light_len_x * Vec3f{1.0f,0.f,0.f}+ 
                          (samples.y) * BoxScene::light_len_y * Vec3f{0.f,1.0f,0.f};
    float pdf = BoxScene::inv_light_area;

    return std::make_tuple(sampled_point, light_normal, pdf);
}


Vec3f shade_with_light_sampling(Triangle tri, Vec3f p, Vec3f wo) {
    Vec3f L_dir{0.0f};
    // Uniformly sample the light at x
    auto [sample_point, light_normal, temp_pdf] = sample_area_light(UniformSampler::next2d());
    // Shoot a ray from p to x
    Vec3f p_ = offset_ray_origin(p, tri.face_normal);
    Vec3f light_dir = p - sample_point;
    float distance_squred = length(light_dir) * length(light_dir);
    auto [hit_, t_min_, nearest_tri_] = RayTracer::closest_hit(p_, -light_dir, octree, BoxScene::triangles);
    light_dir = normalize(light_dir);
    unsigned int material_id = tri.material_id;
    // If the ray is not blocked in the middle
    if(is_emitter(nearest_tri_)){
        float cosine_1 = std::max(0.0f, dot(tri.face_normal, -light_dir));
        float cosine_2 = std::max(0.0f, dot(BoxScene::light_normal, light_dir));
        Vec3f f_r = Vec3f{0.f};
        if(std::holds_alternative<Lambertian>(BoxScene::materials[material_id])){
            f_r = get<Lambertian>(BoxScene::materials[material_id]).eval();
        }
        else if(std::holds_alternative<Cloth>(BoxScene::materials[material_id])){
            f_r = get<Cloth>(BoxScene::materials[material_id]).eval(wo, -light_dir, tri.face_normal, tri.tangent, tri.bitangent);
        }
        Vec3f L_i = eval_area_light(light_dir);
        L_dir = L_i * f_r * cosine_1* cosine_2 / distance_squred / temp_pdf;
    }
    // Contribution from other reflectors
    Vec3f L_indir{0.0f};
    // Test Russian Roulette with probability p_rr = 0.8f
    const float p_rr = 0.8f;
    if(UniformSampler::next1d()>p_rr){
        return L_dir;
    }
    // Randomly choose one direction wi
    Vec3f wi = Vec3f{0.f};
    float pdf = 0.f;
    if(std::holds_alternative<Lambertian>(BoxScene::materials[material_id])){
        auto [wi_temp, pdf_temp] = get<Lambertian>(BoxScene::materials[material_id]).sample(tri.face_normal, UniformSampler::next2d());
        wi = wi_temp;
        pdf = pdf_temp;
    }
    else if(std::holds_alternative<Cloth>(BoxScene::materials[material_id])){
        auto [wi_temp, pdf_temp] = get<Cloth>(BoxScene::materials[material_id]).sample(tri.face_normal, UniformSampler::next2d());
        wi = wi_temp;
        pdf = pdf_temp;
    }

    // Trace the new ray
    auto [hit, t_min, nearest_tri] = RayTracer::closest_hit(p_, wi, octree, BoxScene::triangles);
    if(hit){
        float cosine = std::max(0.0f, dot(tri.face_normal, wi));
        // If the ray hit a non-emitting object at q
        Vec3f f_r = Vec3f{0.f};
        if(std::holds_alternative<Lambertian>(BoxScene::materials[material_id])){
            f_r = get<Lambertian>(BoxScene::materials[material_id]).eval();
        }
        else if(std::holds_alternative<Cloth>(BoxScene::materials[material_id])){
            f_r = get<Cloth>(BoxScene::materials[material_id]).eval(wo, wi, tri.face_normal, tri.tangent, tri.bitangent);
        }
        Vec3f L_i = eval_area_light(light_dir);
        if(!is_emitter(nearest_tri)){
            Vec3f hit_point = p + wi * t_min;
            L_indir = shade_with_light_sampling(nearest_tri, hit_point, -wi) * f_r * cosine / pdf / p_rr;
        }
    }
    
    return L_dir + L_indir;
    // =============================================================================================
}

Vec3f path_tracing_with_light_sampling(Vec3f ray_pos, Vec3f ray_dir) {
    const auto [is_ray_hit, t_min, nearest_tri] =
        RayTracer::closest_hit(ray_pos, ray_dir, octree, BoxScene::triangles);
    if (!is_ray_hit) return Vec3f{0.0f};
    const Vec3f hit_position = ray_pos + t_min * ray_dir;
    if (is_emitter(nearest_tri)) return eval_area_light(-ray_dir);

    return shade_with_light_sampling(nearest_tri, hit_position, -ray_dir);
}









// main

int main(int argc, char **argv) {


    spdlog::info("\n"
                 "----------------------------------------------\n"
                 "Welcome to CS 190I Assignment 4: Cloth Materials\n"
                 "----------------------------------------------");
    const unsigned int max_spp = 64;
    const unsigned int image_width = 1024;
    const unsigned int image_height = 1024;
    // Some prepereations
    Image image{.width = image_width,
                .height = image_height,
                .pixels = std::vector<Vec3f>(image_width * image_height)};
    Camera camera{.vertical_field_of_view = 38.6f,
                  .aspect = static_cast<float>(image_width) / image_height,
                  .focal_distance = 0.8f,
                  .position = Vec3f{0.278f, 0.8f, 0.2744f},
                  .view_direction = Vec3f{0.0f, -1.0f, 0.0f},
                  .up_direction = Vec3f{0.0f, 0.0f, 1.0f},
                  .right_direction = Vec3f{-1.0f, 0.0f, 0.0f}};
    camera.init();
    UniformSampler::init(190);


    // =============================================================================================
    // Change the material ID after you have implemented the Cloth BRDF
    // Diffuse
    const int bunny_material_id = 7;
    // Iron
    // const int bunny_material_id = 5;
    // Gold
    // const int bunny_material_id = 6;
    
    // Load the scene
    // If program can't find the bunny.obj file, use xmake run -w . or move the bunny.obj file to the 
    // same directory as the executable file.
    const std::string obj_path = "./cloth_.obj";
    std::vector<Triangle> obj_triangles = load_obj(obj_path, bunny_material_id);
    BoxScene::triangles.insert(BoxScene::triangles.end(),
                               std::make_move_iterator(obj_triangles.begin()),
                               std::make_move_iterator(obj_triangles.end()));
    octree.build_octree(BoxScene::triangles);
    // =============================================================================================
    // Path Tracing with light sampling
    spdlog::info("Path Tracing with light sampling: rendering started!");
    for (int y = 0; y < image.height; y++) {
        if (y % 50 == 0) {
            spdlog::info("Rendering row {} / {} \r", y, image.height);
        }
        for (int x = 0; x < image.width; x++) {
            image(x, y) = Vec3f{0.0f};
            for (int sample = 0; sample < max_spp; sample++) {
                const float u = (x + UniformSampler::next1d()) / image.width;
                const float v = (y + UniformSampler::next1d()) / image.height;
                Vec3f ray_direction = camera.generate_ray(u, (1.0f - v));
                image(x, y) += clamp(path_tracing_with_light_sampling(
                                         camera.position, ray_direction),
                                     Vec3f(0.0f), Vec3f(50.0f));
            }
            image(x, y) /= (float)max_spp;
        }
    }
    spdlog::info("Path Tracing with light sampling: Rendering finished!");
    image.save_with_tonemapping("./path_tracing_with_light_sampling.png");


    // =============================================================================================
    // Path Tracing with MIS
    // spdlog::info("Path Tracing with MIS: Rendering started!");
    // for (int y = 0; y < image.height; y++) {
    //     if (y % 50 == 0) {
    //         spdlog::info("Rendering row {} / {} \r", y, image.height);
    //     }
    //     for (int x = 0; x < image.width; x++) {
    //         image(x, y) = Vec3f{0.0f};
    //         for (int sample = 0; sample < max_spp; sample++) {
    //             const float u = (x + UniformSampler::next1d()) / image.width;
    //             const float v = (y + UniformSampler::next1d()) / image.height;
    //             Vec3f ray_direction = camera.generate_ray(u, (1.0f - v));
    //             image(x, y) +=
    //                 clamp(path_tracing_with_MIS(camera.position, ray_direction),
    //                       Vec3f(0.0f), Vec3f(50.0f));
    //         }
    //         image(x, y) /= (float)max_spp;
    //     }
    // }
    // spdlog::info("Path Tracing with MIS: Rendering finished!");
    // image.save_with_tonemapping("./path_tracing_with_MIS.png");

    // =============================================================================================
    return 0;
}
