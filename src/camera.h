#pragma once
#include "hittable.h"
#include "material.h"
#include <span>
#include <string_view>
#include <fstream>
#include <vector>

double linear_to_gamma(double linear_component)
{
    return linear_component > 0 ? std::sqrt(linear_component) : 0;
}

bool write_bmp(std::string_view filename, int width, int height,std::span<const uint8_t> pixels) 
{
    unsigned char bmpHeader[] = {'B', 'M', 0, 0, 0, 0, 0, 0, 0, 0, 54, 0, 0, 0};

    unsigned int fileSize = 54 + width * height * 3;
    bmpHeader[2] = static_cast<unsigned char>(fileSize & 0xFF);
    bmpHeader[3] = static_cast<unsigned char>((fileSize >> 8) & 0xFF);
    bmpHeader[4] = static_cast<unsigned char>((fileSize >> 16) & 0xFF);
    bmpHeader[5] = static_cast<unsigned char>((fileSize >> 24) & 0xFF);

    unsigned char dibHeader[40] = {40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 24, 0};
    dibHeader[4] = static_cast<unsigned char>(width & 0xFF);
    dibHeader[5] = static_cast<unsigned char>((width >> 8) & 0xFF);
    dibHeader[8] = static_cast<unsigned char>(height & 0xFF);
    dibHeader[9] = static_cast<unsigned char>((height >> 8) & 0xFF);

    std::ofstream bmpFile(filename.data(), std::ios::binary);
    if (!bmpFile.is_open()) return false;
    
    bmpFile.write(reinterpret_cast<const char*>(bmpHeader), sizeof(bmpHeader));
    bmpFile.write(reinterpret_cast<const char*>(dibHeader), sizeof(dibHeader));
    bmpFile.write(reinterpret_cast<const char*>(pixels.data()), width * height * 3);
    return true;
}

struct camera 
{
    double aspect_ratio      = 1.0;
    int    image_width       = 100;
    int    samples_per_pixel = 10;
    int    max_depth         = 10;

    double vfov     = 90;
    vec3 lookfrom = vec3(0,0,0);
    vec3 lookat   = vec3(0,0,-1);
    vec3   vup      = vec3(0,1,0);

    double defocus_angle = 0;
    double focus_dist = 10; 

    void render(const hittable& world) 
    {
        initialize();

        std::vector<std::uint8_t> pixels(image_width * image_height * 3);
        for (int j = 0; j < image_height; j++) 
        {
           
            for (int i = 0; i < image_width; i++) 
            {
                vec3 pixel_vec3(0,0,0);
                for (int sample = 0; sample < samples_per_pixel; sample++) 
                {
                    ray r = get_ray(i, j);
                    pixel_vec3 += ray_vec3(r, max_depth, world);
                }
                vec3 clr = pixel_samples_scale * pixel_vec3;
                auto r = clr.x();
                auto g = clr.y();
                auto b = clr.z();

                r = linear_to_gamma(r);
                g = linear_to_gamma(g);
                b = linear_to_gamma(b);
                static const interval intensity(0.000, 0.999);
                int index = ((image_height - 1 - j) * image_width + i) * 3;
                pixels[index + 0] = int(255.99 * intensity.clamp(b));
                pixels[index + 1] = int(255.99 * intensity.clamp(g));
                pixels[index + 2] = int(255.99 * intensity.clamp(r));
            }
        }
        write_bmp("test.bmp",image_width,image_height,pixels);
    }

  private:
    int    image_height;
    double pixel_samples_scale;
    vec3 center;
    vec3 pixel00_loc;
    vec3   pixel_delta_u; 
    vec3   pixel_delta_v;
    vec3   u, v, w;
    vec3   defocus_disk_u;
    vec3   defocus_disk_v;

    void initialize() {
        image_height = int(image_width / aspect_ratio);
        image_height = (image_height < 1) ? 1 : image_height;

        pixel_samples_scale = 1.0 / samples_per_pixel;

        center = lookfrom;

        // Determine viewport dimensions.
        auto theta = degrees_to_radians(vfov);
        auto h = std::tan(theta/2);
        auto viewport_height = 2 * h * focus_dist;
        auto viewport_width = viewport_height * (double(image_width)/image_height);
        w = unit_vector(lookfrom - lookat);
        u = unit_vector(cross(vup, w));
        v = cross(w, u);


        vec3 viewport_u = viewport_width * u;
        vec3 viewport_v = viewport_height * -v;

        pixel_delta_u = viewport_u / image_width;
        pixel_delta_v = viewport_v / image_height;
        auto viewport_upper_left = center - (focus_dist * w) - viewport_u/2 - viewport_v/2;
        pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v);
        auto defocus_radius = focus_dist * std::tan(degrees_to_radians(defocus_angle / 2));
        defocus_disk_u = u * defocus_radius;
        defocus_disk_v = v * defocus_radius;
    }

    ray get_ray(int i, int j) const 
    {

        auto offset = sample_square();
        auto pixel_sample = pixel00_loc
                          + ((i + offset.x()) * pixel_delta_u)
                          + ((j + offset.y()) * pixel_delta_v);

        auto ray_origin = (defocus_angle <= 0) ? center : defocus_disk_sample();
        auto ray_direction = pixel_sample - ray_origin;

        return ray(ray_origin, ray_direction);
    }

    vec3 sample_square() const {
        return vec3(random_double() - 0.5, random_double() - 0.5, 0);
    }

    vec3 sample_disk(double radius) const {
        return radius * random_in_unit_disk();
    }

    vec3 defocus_disk_sample() const {
        auto p = random_in_unit_disk();
        return center + (p[0] * defocus_disk_u) + (p[1] * defocus_disk_v);
    }

    vec3 ray_vec3(const ray& r, int depth, const hittable& world) const {
        if (depth <= 0)
            return vec3(0,0,0);

        hit_record rec;

        if (world.hit(r, interval(0.001, infinity), rec)) {
            ray scattered;
            vec3 attenuation;
            if (rec.mat->scatter(r, rec, attenuation, scattered))
                return attenuation * ray_vec3(scattered, depth-1, world);
            return vec3(0,0,0);
        }

        vec3 unit_direction = unit_vector(r.direction());
        auto a = 0.5*(unit_direction.y() + 1.0);
        return (1.0-a)*vec3(1.0, 1.0, 1.0) + a*vec3(0.5, 0.7, 1.0);
    }
};

