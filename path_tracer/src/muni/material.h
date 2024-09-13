#pragma once
#include "common.h"
#include "math_helpers.h"
#include <algorithm>
#include <cmath>
#include <winuser.h>

namespace muni {
struct Lambertian {
    Vec3f albedo;

    /** Evaluates the BRDF for the Lambertian material.
      \return The BRDF (fr) value.
  */
    Vec3f eval() const { 
        return albedo / M_PI;
     }

    /** Samples the BRDF for the Lambertian material.
      \param[in] normal The normal of the surface.
      \param[in] u A random number in (0,1)^2.
      \return A tuple containing the sampled direction in world space and the PDF.
    */
    std::tuple<Vec3f, float> sample(Vec3f normal, Vec2f u) const {
        float z = u.x;
        float r = sqrt(1-z*z);
        float phi = 2 * M_PI * u.y;
        float pdf = 1.0f / (2.0f * M_PI);
        Vec3f direction = {r*cos(phi),r*sin(phi),z};
        direction = from_local(direction, normal);
        return std::make_tuple(direction, pdf);
    }
    /** Computes the PDF for the Lambertian material.
      \param[in] wo The outgoing direction in world space.
      \param[in] wi The light incident direction in world space.
      \return The PDF value.
    */
    float pdf(Vec3f wo_world, Vec3f wi_world, Vec3f normal) const {
        return M_1_2PI;
    }

};

struct Microfacet {
    float roughness;
    // refraction indices for RGB channels
    Vec3f n1;
    Vec3f n2;

    /** Computes the Fresnel term for the microfacet material.
      \param[in] wi The light incident direction in local space.
      \return The Fresnel term.
    */
    Vec3f F(Vec3f wi) const {
        // =============================================================================================
        // TODO: Implement this function
        Vec3f R0 = pow((n1-n2)/(n1+n2),2.0f);
        float cosTheta = wi.z;
        return R0 + (Vec3f{1.0f} - R0) * pow(1.0f - cosTheta, 5.0f);
        // =============================================================================================
    }
    /** Computes the Beckmann normal distribution function for the microfacet material.
      \param[in] h The half vector in local space.
      \return The normal distribution function.
    */
    float D(Vec3f h) const {
        // =============================================================================================
        // TODO: Implement this function
        float alpha2 = roughness * roughness;
        float cosThetaH = h.z;
        float tanThetaH2 = (1.0f - cosThetaH * cosThetaH) / (cosThetaH * cosThetaH);
        // spdlog::warn(tanThetaH2 / alpha2);
        return exp(-tanThetaH2 / alpha2) / (M_PI * alpha2 * cosThetaH * cosThetaH * cosThetaH * cosThetaH);
        // return 0.02f / (M_PI * alpha2 * cosThetaH * cosThetaH * cosThetaH * cosThetaH);
        // return 1.0f;
        // =============================================================================================
    }

    float G1(Vec3f v) const {
        float tanTheta = sqrt((1.0f - v.z*v.z)/v.z/v.z);
        if (tanTheta == 0.0f) return 1.0f;
        float a = 1.0f / (roughness * tanTheta);
        if (a >= 1.6f) return 1.0f;
        return 1.0f / (1.0f + (1.0f - 1.259f * a + 0.396f * a * a) / (3.535f * a + 2.181f * a * a));
    }
    /** Computes the shadowing-masking function for the microfacet material.
      \param[in] wo The outgoing direction in local space.
      \param[in] wi The light incident direction in local space.
      \return The shadowing-masking value.
    */
    float G(Vec3f wo, Vec3f wi) const {
        // =============================================================================================
        // TODO: Implement this function
        return G1(wo) * G1(wi);
        // =============================================================================================
    }
    /** Evaluates the BRDF for the microfacet material.
      \param[in] wo_world The outgoing direction in world space.
      \param[in] wi_world The light incident direction in world space.
      \param[in] normal The normal of the surface.
      \return The BRDF (fr) value.
    */
    Vec3f eval(Vec3f wo_world, Vec3f wi_world, Vec3f normal) const {
        Vec3f wo = to_local(wo_world, normal);
        Vec3f wi = to_local(wi_world, normal);
        // =============================================================================================
        // TODO: Implement this function
        Vec3f h = normalize((wo + wi) / 2);
        
        Vec3f F = this->F(wi);
        float D = this->D(h);
        float G = this->G(wo, wi);
        // if(D!=0) spdlog::warn(D);
        // spdlog::warn(F);
        // spdlog::warn(D);
        // spdlog::warn(G);

        float denom = 4.0f * wi.z * wo.z;
        // spdlog::warn("r:{}",F * G / denom);
        return F * D * G / denom;
        // =============================================================================================
    }
    /** Computes the PDF for the microfacet material.
      \param[in] wo The outgoing direction in world space.
      \param[in] wi The light incident direction in world space.
      \param[in] normal The normal of the surface.
      \return The PDF value.
    */
    float pdf(Vec3f wo_world, Vec3f wi_world, Vec3f normal) const {
        // =============================================================================================
        // TODO: Implement this function
        Vec3f wo = to_local(wo_world, normal);
        Vec3f wi = to_local(wi_world, normal);
        Vec3f wh = normalize(wo + wi);

        float alpha = roughness;
        float theta_h = wh.z;
        float p_theta = 2.0f * sin(theta_h) / (pow(alpha, 2.0f) * pow(cos(theta_h), 3.0f)) * exp(-tan(theta_h) * tan(theta_h) / alpha / alpha);
        float p_phi = INV_TWOPI;
        float pdf_wh = p_theta * p_phi / sin(theta_h);
        float pdf_wi = pdf_wh / (4.0f * dot(wo, wh));
        
        return pdf_wi;
        // =============================================================================================
    }

    /** Samples the BRDF for the microfacet material.
      \param[in] wo_world The outgoing direction in world space.
      \param[in] normal The normal of the surface.
      \param[in] u A random number in (0,1)^2.
      \return A tuple containing the sampled direction in world space and the PDF.
    */
     std::tuple<Vec3f, float> sample(Vec3f wo_world, Vec3f normal, Vec2f u) const {
        // =============================================================================================
        // TODO: Implement this function
        Vec3f wo = to_local(wo_world, normal);
        float alpha = roughness;

        float phi_h = 2.0f * M_PI * u.x;
        float theta_h = atan(sqrt(-alpha * alpha * log(1.0f - u.y)));

        Vec3f wh = Vec3f(sin(theta_h) * cos(phi_h), sin(theta_h) * sin(phi_h), cos(theta_h));
        // spdlog::warn(wh);
        Vec3f wi = mirror_reflect(-wo, wh);
        // spdlog::warn(wi);
        // spdlog::warn(wo);
        Vec3f wi_world = from_local(wi, normal);

        float p_theta = 2.0f * sin(theta_h) / (pow(alpha, 2.0f) * pow(cos(theta_h), 3.0f)) * exp(-tan(theta_h) * tan(theta_h) / alpha / alpha);
        float p_phi = INV_TWOPI;
        float pdf_wh = p_theta * p_phi / sin(theta_h);
        float pdf_wi = pdf_wh / (4.0f * dot(wo, wh));
        // spdlog::warn(wi_world);
        // spdlog::warn(pdf_wi);
        // wi_world = from_local(Vec3f{0.f,0.f,1.0f},normal);
        return std::make_tuple(wi_world, pdf_wi);
        // =============================================================================================
      }

};
}  // namespace muni
