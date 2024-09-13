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

struct Cloth {
    float eta_0;
    float kd_0;
    float gamma_s_0;
    float gamma_v_0;
    Vec3f A_0;
    float alpha_0;
    std::vector<float> tangent_offset_u;

    float eta_1;
    float kd_1;
    float gamma_s_1;
    float gamma_v_1;
    Vec3f A_1;
    float alpha_1;
    std::vector<float> tangent_offset_v;

    Vec3f eval(Vec3f wo_world, Vec3f wi_world, Vec3f normal_, Vec3f tangent_, Vec3f bitangent_) const {
        // Vec3f wo = to_local(wo_world, normal_);
        // Vec3f wi = to_local(wi_world, normal_);
        Vec3f wo = wo_world;
        Vec3f wi = wi_world;

        Vec3f normal = normalize(normal_);
        Vec3f bitangent = normalize(cross(normal, normalize(tangent_)));
        Vec3f tangent = normalize(cross(bitangent, normal));

        // normal = to_local(normal, normal_);
        // bitangent = to_local(bitangent, normal_);
        // tangent = to_local(tangent, normal_);
        float Q = 0.0f;
        Vec3f fs = Vec3f{0.0f};
        Vec3f u_value = Vec3f{0.0f};
        Vec3f v_value = Vec3f{0.0f};
        for(auto tangent_offset:tangent_offset_u){
          
          float cos_temp = abs(dot(wi, tangent));
          float theta_i = asin(cos_temp);
          cos_temp = abs(dot(wo, tangent));
          float theta_o = asin(cos_temp);
          float cos_theta_i = cos(theta_i);
          float cos_theta_o = cos(theta_o);
          float theta_h = (theta_i + theta_o) / 2;
          float theta_d = (theta_i - theta_o) / 2;

          Vec3f wi_proj = wi - dot(wi, tangent) * tangent;
          Vec3f wo_proj = wo - dot(wo, tangent) * tangent;
          float cos_phi_i = dot(wi_proj, normal);
          float phi_i = acos(cos_phi_i);
          float cos_phi_o = dot(wo_proj, normal);
          float phi_o = acos(cos_phi_o);
          float phi_d = abs(phi_i - phi_o);

          Vec3f normal_temp = normalize(cross(normal, tangent));
          Vec3f wi_proj_ = wi - dot(wi, normal_temp) * normal_temp;
          Vec3f wo_proj_ = wo - dot(wo, normal_temp) * normal_temp;
          float cos_psi_i = dot(wi_proj_, normal);
          float cos_psi_o = dot(wo_proj_, normal);
          float psi_i = acos(cos_psi_i);
          float psi_o = acos(cos_psi_o);
          float psi_d = abs(psi_i - psi_o);

          float eta_i = 1.0f;
          // float Fr_cos_theta_i = cos(theta_d) * cos(phi_d * 0.5f);
          float Fr = fresnel(eta_0, eta_i, wi, normal);
          float Fr_2 = fresnel(eta_0, eta_i, wo, normal);
          float fr_s = Fr * cos(phi_d * 0.5f) * normalized_gaussian(gamma_s_0, theta_h);

          float Ft = 1.0f - Fr;
          float Ft_2 = 1.0f - Fr_2;
          float F = Ft * Ft_2;
          Vec3f fr_v = F * ((1.0f - kd_0) * normalized_gaussian(gamma_v_0, theta_h) + kd_0)
                        / (cos_theta_i + cos_theta_o) 
                        * A_0;

          fs = (Vec3f{fr_s, fr_s, fr_s} + fr_v) / pow(cos(theta_d), 2.0f);

          float m_i = std::max(cos_phi_i, 0.0f);
          float m_o = std::max(cos_phi_o, 0.0f);
          float u = u_gaussian(phi_d);
          float masking = (1.0f - u) * m_i * m_o + u * std::min(m_i, m_o);

          float p_i = std::max(cos_psi_i, 0.0f);
          float p_o = std::max(cos_psi_o, 0.0f);
          u = u_gaussian(psi_d);
          float p = (1.0f - u) * p_i * p_o + u * std::min(p_i, p_o);

          u_value += p * masking * fs * cos_theta_i;
          Q += alpha_0 * p / tangent_offset_u.size();
          // spdlog::warn(cos_phi_i);
        }
        u_value /= tangent_offset_u.size();


        tangent = normalize(cross(normal, tangent));

        for(auto tangent_offset:tangent_offset_v){
          
          Vec3f wi_proj = wi - dot(wi, tangent) * tangent;
          Vec3f wo_proj = wo - dot(wo, tangent) * tangent;
          float cos_temp = abs(dot(wi, tangent));
          float theta_i = asin(cos_temp);
          cos_temp = abs(dot(wo, tangent));
          float theta_o = asin(cos_temp);
          float cos_theta_i = cos(theta_i);
          float cos_theta_o = cos(theta_o);
          float theta_h = (theta_i + theta_o) / 2;
          float theta_d = (theta_i - theta_o) / 2;

          float cos_phi_i = abs(dot(wi_proj, normal));
          float phi_i = acos(cos_phi_i);
          float cos_phi_o = abs(dot(wo_proj, normal));
          float phi_o = acos(cos_phi_o);
          float phi_d = abs(phi_i - phi_o);

          Vec3f normal_temp = normalize(cross(normal, tangent));
          Vec3f wi_proj_ = wi - dot(wi, normal_temp) * normal_temp;
          Vec3f wo_proj_ = wo - dot(wo, normal_temp) * normal_temp;
          float cos_psi_i = dot(wi_proj_, normal);
          float cos_psi_o = dot(wo_proj_, normal);
          float psi_i = acos(cos_psi_i);
          float psi_o = acos(cos_psi_o);
          float psi_d = abs(psi_i - psi_o);

          float eta_i = 1.0f;
          // float Fr_cos_theta_i = cos(theta_d) * cos(phi_d * 0.5f);
          float Fr = fresnel(eta_1, eta_i, wi, normal);
          float Fr_2 = fresnel(eta_1, eta_i, wo, normal);
          float fr_s = Fr * cos(phi_d * 0.5f) * normalized_gaussian(gamma_s_1, theta_h);

          float Ft = 1.0f - Fr;
          float Ft_2 = 1.0f - Fr_2;
          float F = Ft * Ft_2;
          Vec3f fr_v = F * ((1.0f - kd_1) * normalized_gaussian(gamma_v_1, theta_h) + kd_1)
                        / (cos_theta_i + cos_theta_o) 
                        * A_1;

          fs = (Vec3f{fr_s, fr_s, fr_s} + fr_v) / pow(cos(theta_d), 2.0f);

          float m_i = std::max(cos_phi_i, 0.0f);
          float m_o = std::max(cos_phi_o, 0.0f);
          float u = u_gaussian(phi_d);
          float masking = (1.0f - u) * m_i * m_o + u * std::min(m_i, m_o);

          float p_i = std::max(cos_psi_i, 0.0f);
          float p_o = std::max(cos_psi_o, 0.0f);
          u = u_gaussian(psi_d);
          float p = (1.0f - u) * p_i * p_o + u * std::min(p_i, p_o);

          v_value += p * masking * fs * cos_theta_i;
          Q += alpha_1 * p / tangent_offset_v.size();
        }
        v_value /= tangent_offset_v.size();

        fs = u_value * alpha_0 + v_value * alpha_1;
        // spdlog::warn(u_value);
        Q += (1.0f - alpha_0 - alpha_1) * dot(wi, normal);
        if(Q > 0.0f){
            fs /= Q;
        }
        return fs;
    }

     std::tuple<Vec3f, float> sample(Vec3f normal, Vec2f u) const {
        float z = u.x;
        float r = sqrt(1-z*z);
        float phi = 2 * M_PI * u.y;
        float pdf = 1.0f / (2.0f * M_PI);
        Vec3f direction = {r*cos(phi),r*sin(phi),z};
        direction = from_local(direction, normal);
        return std::make_tuple(direction, pdf);
      }



    // Helper

    float normalized_gaussian(float beta, float theta) const {
        return exp(-theta * theta / (2.0f * beta * beta)) / (sqrt(2.0f * M_PI * beta * beta));
    }

    float fresnel(float eta_t, float eta_i, Vec3f wi, Vec3f n) const {
        bool is_total_internal = (1.0f - pow(eta_i / eta_t, 2) * (1.0f - pow(dot(wi, n), 2))) < 0.0f;
        if (is_total_internal) {
            return 1.0f;
        }
        Vec3f reflacted_dir = (wi - dot(wi, n) * n) * (eta_i / eta_t) - n * sqrt(1.0f - pow(eta_i / eta_t, 2.0f) * (1.0f - pow(dot(wi, n), 2.0f)));
        float cos_theta_i = dot(wi, n);
        float cos_theta_o = dot(reflacted_dir, -n);
        float rho_s = (eta_i * cos_theta_i - eta_t * cos_theta_o) / (eta_i * cos_theta_i + eta_t * cos_theta_o);
        float rho_p = (eta_i * cos_theta_o - eta_t * cos_theta_i) / (eta_i * cos_theta_o + eta_t * cos_theta_i);
        float fresnel = 0.5f * (pow(rho_s, 2) + pow(rho_p, 2));
        return fresnel;
    }

    float u_gaussian(float x) const {
        float sd = 20.0 * M_PI / 180.0; // Converting degrees to radians
        return exp(-x * x / (2.0 * sd * sd));
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
