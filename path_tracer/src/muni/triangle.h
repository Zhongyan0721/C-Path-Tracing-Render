#pragma once
#include "common.h"
#include "math_helpers.h"
#include <tuple>

namespace muni {
struct Triangle {
    Vec3f v0, v1, v2;
    Vec3f face_normal;
    Vec3f emission;
    unsigned int material_id;

    /** Ray-Triangle intersection based on "Watertight Ray/Triangle Intersection"
        Paper link: http://jcgt.org/published/0002/01/05/paper.pdf
        \param[in] tri The triangle to intersect with.
        \param[in] ray_origin The origin of the ray.
        \param[in] ray_direction The direction of the ray.
        \param[in] t_min The minimum t value of the intersection point along the ray.
        \param[in] t_max The maximum t value of the intersection point along the ray.
        \return A tuple containing a boolean indicating whether the ray intersects
        the triangle, and the t value of the intersection point along the ray.
    */
    static std::tuple<bool, float>
    ray_triangle_intersect(Triangle tri, Vec3f ray_origin, Vec3f ray_direction,
                           float t_min, float t_max) {
        const Vec3f abs_ray_direction = abs(ray_direction);
        unsigned int axis = 0;
        if (abs_ray_direction[1] > abs_ray_direction[0] &&
            abs_ray_direction[1] > abs_ray_direction[2])
            axis = 1;
        if (abs_ray_direction[2] > abs_ray_direction[0] &&
            abs_ray_direction[2] > abs_ray_direction[1])
            axis = 2;

        unsigned int kz = axis;
        unsigned int kx = (kz + 1) % 3;
        unsigned int ky = (kx + 1) % 3;
        if (ray_direction[kz] < 0.0f) {
            unsigned int swap = kx;
            kx = ky;
            ky = swap;
        }

        float Sx = ray_direction[kx] / ray_direction[kz];
        float Sy = ray_direction[ky] / ray_direction[kz];
        float Sz = 1.f / ray_direction[kz];

        const Vec3f A = tri.v0 - ray_origin;
        const Vec3f B = tri.v1 - ray_origin;
        const Vec3f C = tri.v2 - ray_origin;

        const float Ax = A[kx] - Sx * A[kz];
        const float Ay = A[ky] - Sy * A[kz];
        const float Bx = B[kx] - Sx * B[kz];
        const float By = B[ky] - Sy * B[kz];
        const float Cx = C[kx] - Sx * C[kz];
        const float Cy = C[ky] - Sy * C[kz];

        float U = Cx * By - Cy * Bx;
        float V = Ax * Cy - Ay * Cx;
        float W = Bx * Ay - By * Ax;

        if (U == 0.f || V == 0.f || W == 0.f) {
            double CxBy = static_cast<double>(Cx) * static_cast<double>(By);
            double CyBx = static_cast<double>(Cy) * static_cast<double>(Bx);
            U = (float)(CxBy - CyBx);
            double AxCy = static_cast<double>(Ax) * static_cast<double>(Cy);
            double AyCx = static_cast<double>(Ay) * static_cast<double>(Cx);
            V = (float)(AxCy - AyCx);
            double BxAy = static_cast<double>(Bx) * static_cast<double>(Ay);
            double ByAx = static_cast<double>(By) * static_cast<double>(Ax);
            W = (float)(BxAy - ByAx);
        }

        if ((U < 0.f || V < 0.f || W < 0.f) && (U > 0.f || V > 0.f || W > 0.f))
            return {false, 0.0f};

        float det = U + V + W;
        if (det == 0.f) return {false, 0.0f};

        const float Az = Sz * A[kz];
        const float Bz = Sz * B[kz];
        const float Cz = Sz * C[kz];
        const float T = U * Az + V * Bz + W * Cz;
        const float rcp_det = 1.f / det;

        const float t = T * rcp_det;
        const Vec3f barycentrics = Vec3f(U, V, W) * rcp_det;

        if (t < t_min || t > t_max) return {false, 0.0f};

        return {true, t};
    }
};

}  // namespace muni
