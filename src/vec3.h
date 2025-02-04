#pragma once
#include <cmath>
#include <array>
#include "utils.h"


struct vec3 
{ 
    vec3() = default;
    vec3(float e0, float el, float e2) 
    { 
        e[0] = e0; 
        e[1] = el; 
        e[2] = e2; 
    } 
    float x() const { return e[0]; }
    float y() const { return e[1]; } 
    float z() const { return e[2]; } 
    float r() const { return e[0]; } 
    float g() const { return e[1]; } 
    float b() const { return e[2]; }
    const vec3& operator+() const { return *this; }
    vec3 operator-() const { return vec3(-e[0], -e[1], -e[2]); } 
    float operator[] (int i) const { return e[i]; }
    float& operator[] (int i) { return e[i]; };

    vec3& operator+= (const vec3 &v2)
    {
        e[0] += v2.e[0];
        e[1] += v2.e[1];
        e[2] += v2.e[2];
        return *this;
    }
    vec3& operator-=(const vec3 &v2)
    {
        e[0] -= v2.e[0];
        e[1] -= v2.e[1];
        e[2] -= v2.e[2];
        return *this;
    }
    vec3& operator*= (const vec3 &v2)
    {
        e[0] *= v2.e[0];
        e[1] *= v2.e[1];
        e[2] *= v2.e[2];
        return *this;
    }
    vec3& operator /= (const vec3 &v2)
    {
        e[0] /= v2.e[0];
        e[1] /= v2.e[1];
        e[2] /= v2.e[2];
        return *this;
    }
    vec3& operator*= (const float t)
    {
        e[0] *= t;
        e[1] *= t;
        e[2] *= t;
        return *this;
    }
    vec3& operator/=(const float t)
    {
        e[0] /= t;
        e[1] /= t;
        e[2] /= t;
        return *this;
    }
    static vec3 random() {
        return vec3(random_double(), random_double(), random_double());
    }

    static vec3 random(double min, double max) {
        return vec3(random_double(min,max), random_double(min,max), random_double(min,max));
    }
    float length() const noexcept
    {
        return sqrt(length_squared()); 
    } 
    float length_squared() const noexcept
    {
        return e[0]*e[0] + e[1]*e[1] + e[2]*e[2]; 
    }
    void make_unit_vector()
    {
        float k = 1.0 / sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
        e[0] *= k; e[1] *= k; e[2] *= k;
    }

    bool near_zero() const 
    {
        auto s = 1e-8;
        return (std::fabs(e[0]) < s) && (std::fabs(e[1]) < s) && (std::fabs(e[2]) < s);
    }

    std::array<float,3> e{};
};



vec3 operator+(const vec3 &v1, const vec3 &v2) 
{
    return vec3(v1.e[0] + v2.e[0], v1.e[1] + v2.e[1], v1.e[2] + v2.e[2]);
}

vec3 operator-(const vec3 &v1, const vec3 &v2) 
{
    return vec3(v1.e[0] - v2.e[0], v1.e[1] - v2.e[1], v1.e[2] - v2.e[2]);
}

vec3 operator*(const vec3 &v1, const vec3 &v2) 
{
    return vec3(v1.e[0] * v2.e[0], v1.e[1] * v2.e[1], v1.e[2] * v2.e[2]);
}

vec3 operator/(const vec3 &v1, const vec3 &v2) 
{
    return vec3(v1.e[0] / v2.e[0], v1.e[1] / v2.e[1], v1.e[2] / v2.e[2]);
}

vec3 operator*(float t, const vec3 &v) 
{
    return v * vec3(t,t,t);
}
vec3 operator*(const vec3 &v, float t) 
{
    return t * v;
}

vec3 operator/(vec3 v, float t) 
{
    return vec3(v.e[0]/t, v.e[1]/t, v.e[2]/t);
}

float dot(const vec3 &v1, const vec3 &v2) 
{
    return v1.e[0] * v2.e[0]+ v1.e[1] * v2.e[1]+ v1.e[2] * v2.e[2];
}

vec3 cross(const vec3 &v1, const vec3 &v2) {
    return vec3(v1.e[1] * v2.e[2] - v1.e[2] * v2.e[1],
                v1.e[2] * v2.e[0] - v1.e[0] * v2.e[2],
                v1.e[0] * v2.e[1] - v1.e[1] * v2.e[0]);
}


vec3 unit_vector(vec3 v) {
    return v / v.length();
}

vec3 random_in_unit_disk() {
    while (true) {
        auto p = vec3(random_double(-1,1), random_double(-1,1), 0);
        if (p.length_squared() < 1)
            return p;
    }
}

vec3 random_unit_vector() {
    while (true) {
        auto p = vec3::random(-1,1);
        auto lensq = p.length_squared();
        if (1e-160 < lensq && lensq <= 1.0)
            return p / sqrt(lensq);
    }
}

vec3 random_on_hemisphere(const vec3& normal) {
    vec3 on_unit_sphere = random_unit_vector();
    if (dot(on_unit_sphere, normal) > 0.0) 
        return on_unit_sphere;
    else
        return -on_unit_sphere;
}

vec3 reflect(const vec3& v, const vec3& n) {
    return v - 2*dot(v,n)*n;
}

vec3 refract(const vec3& uv, const vec3& n, double etai_over_etat) {
    auto cos_theta = std::fmin(dot(-uv, n), 1.0);
    vec3 r_out_perp =  etai_over_etat * (uv + cos_theta*n);
    vec3 r_out_parallel = -std::sqrt(std::fabs(1.0 - r_out_perp.length_squared())) * n;
    return r_out_perp + r_out_parallel;
}
struct ray 
{
    ray() {}

    ray(const vec3& origin, const vec3& direction) : orig(origin), dir(direction) {}

    const vec3& origin() const  { return orig; }
    const vec3& direction() const { return dir; }

    vec3 at(double t) const {
        return orig + t*dir;
    }

  private:
    vec3 orig;
    vec3 dir;
};


