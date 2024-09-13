#pragma once
#include <vector>
#include <string>
#include <cstdint>

#include "common.h"

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

namespace muni {
/** An image is just a 2D array of pixels.
*/
struct Image {
    int width, height;
    std::vector<Vec3f> pixels;

    /** Tone map the image using the ACES curve.
        See https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
        \param[in] value The value to tone map.
        \return The tone mapped value.
    */
    Vec3f tone_map_Aces(const Vec3f value) const {
        Vec3f color = 0.6f * value;
        float A = 2.51;
        float B = 0.03;
        float C = 2.43;
        float D = 0.59;
        float E = 0.14;
        color = (color * (A * color + B)) / (color * (C * color + D) + E);
        return color;
    }

    /** Save the image to a file.
        \param[in] filename The name of the file to save the image to.
        \return True if the image was saved successfully, false otherwise.
    */
    bool save(const std::string &filename) const {
        // Convert the floating-point data to 8-bit per channel.
        std::vector<uint8_t> data(3 * width * height);
        for (int i = 0; i < width * height; i++) {
            for (int j = 0; j < 3; j++) {
                data[3 * i + j] = static_cast<uint8_t>(
                    255.0f * std::max(0.0f, std::min(1.0f, pixels[i][j])));
            }
        }
        // Save the image to a png file.
        return stbi_write_png(filename.c_str(), width, height, 3, data.data(),
                              sizeof(uint8_t) * 3 * width) != 0;
    }
    /** Save the image to a file with tonemapping.
        \param[in] filename The name of the file to save the image to.
        \return True if the image was saved successfully, false otherwise.
    */
    bool save_with_tonemapping(const std::string &filename) const {
        // Convert the floating-point data to 8-bit per channel.
        std::vector<uint8_t> data(3 * width * height);
        for (int i = 0; i < width * height; i++) {
            for (int j = 0; j < 3; j++) {
                Vec3f pixel = tone_map_Aces(pixels[i]);
                data[3 * i + j] = static_cast<uint8_t>(
                    // 255.0f * std::max(0.0f, std::min(1.0f, pixels[i][j])));
                    255.0f * std::max(0.0f, std::min(1.0f, pixel[j])));
            }
        }
        // Save the image to a png file.
        return stbi_write_png(filename.c_str(), width, height, 3, data.data(),
                              sizeof(uint8_t) * 3 * width) != 0;
    }
    Vec3f &operator()(int x, int y) { return pixels[y * width + x]; }
    const Vec3f &operator()(int x, int y) const {
        return pixels[y * width + x];
    }
};
}  // namespace muni
