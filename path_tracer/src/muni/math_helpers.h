#pragma once
#include "common.h"

#undef M_PI
#define M_PI 3.14159265358979323846f
#define INV_PI 0.31830988618379067154f
#define INV_TWOPI 0.15915494309189533577f
#define INV_FOURPI 0.07957747154594766788f
#define SQRT_TWO 1.41421356237309504880f
#define INV_SQRT_TWO 0.70710678118654752440f

#define M_PI_2 1.57079632679489661923    // pi/2
#define M_PI_4 0.785398163397448309616   // pi/4
#define M_1_2PI 0.159154943091895335769  // 1/2pi

#define EPS 0.001f
#define ANYHIT_EPS 0.005f
namespace muni {

template<typename T> T length_squared(Vec3<T> v) {
    return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

/** Create a coordinate system from a single vector.
    \param[in] v1 The input vector (normal).
    \return A tuple containing the two vectors that form the coordinate system.
 */
std::tuple<Vec3f, Vec3f> coordinate_system(Vec3f v1) {
    float sign = std::copysign(1.0f, v1.z);
    float a = -1 / (sign + v1.z);
    float b = v1.x * v1.y * a;
    return {Vec3f(1 + sign * v1.x * v1.x * a, sign * b, -sign * v1.x),
            Vec3f(b, sign + v1.y * v1.y * a, -v1.y)};
}

/** Transform a vector from the local space to the world space.
    \param[in] v Vector in the local space.
    \param[in] n Normal of the surface.
    \return Vector in the world space.
 */
Vec3f from_local(Vec3f v, Vec3f n) {
    auto [x, y] = coordinate_system(n);
    return v.x * x + v.y * y + v.z * n;
}
Vec3f to_local(const Vec3f v, Vec3f n) { 
    auto [x, y] = coordinate_system(n);
    return Vec3f(dot(v, x), dot(v, y), dot(v, n)); 
}

/** Reflect a ray direction using the surface normal.
    \param[in] incident_dir The incident ray direction (heading to the surface). 
    \param[in] normal The normal of the surface at the hit point.
    \return The reflected ray direction.
*/
Vec3f mirror_reflect(const Vec3f incident_dir, const Vec3f normal) {
    return incident_dir - 2 * dot(incident_dir, normal) * normal;
}
}  // namespace muni
