#pragma once
#include "common.h"
#include "math_helpers.h"

namespace muni {
struct Camera {
    float vertical_field_of_view;
    float aspect;
    float focal_distance;
    Vec3f position;
    Vec3f view_direction;
    Vec3f up_direction;
    Vec3f right_direction;

    Vec3f bottom_left_corner;
    Vec3f up_vector;
    Vec3f right_vector;

    /** Initializes the camera.
     */
    void init() {
        Vec3f middle_of_image_plane =
            position + focal_distance * view_direction;
        float theta = vertical_field_of_view * M_PI / 180.0f;
        float image_plane_height =
            2.0f * focal_distance * std::tan(0.5f * theta);
        float image_plane_width = image_plane_height * aspect;
        up_vector = image_plane_height * up_direction;
        right_vector = image_plane_width * right_direction;
        bottom_left_corner =
            middle_of_image_plane - 0.5f * right_vector - 0.5f * up_vector;
    }

    /** Generates a ray direction from the camera.
        \param[in] u The horizontal coordinate on the image plane.
        \param[in] v The vertical coordinate on the image plane.
        \return The ray direction.
    */
    Vec3f generate_ray(float u, float v) const {
        Vec3f position_on_image_plane =
            bottom_left_corner + u * right_vector + v * up_vector;
        Vec3f ray_direction = normalize(position_on_image_plane - position);
        return ray_direction;
    }
};
}  // namespace muni
