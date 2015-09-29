#pragma once


#include "renderer_math.h"


class Random {
	unsigned int seed_[4];
public:
	unsigned int next(void) {
		const unsigned int t = seed_[0] ^ (seed_[0] << 11);
		seed_[0] = seed_[1];
		seed_[1] = seed_[2];
		seed_[2] = seed_[3];
		return seed_[3] = (seed_[3] ^ (seed_[3] >> 19)) ^ (t ^ (t >> 8));
	}

	double next01(void) {
		return (double)next() / UINT_MAX;
	}

	Random(const unsigned int initial_seed) {
		unsigned int s = initial_seed;
		for (int i = 1; i <= 4; i++){
			seed_[i - 1] = s = 1812433253U * (s ^ (s >> 30)) + i;
		}
	}
};

struct Vec {
	double x = 0.0, y = 0.0, z = 0.0;                  // position, also color (r,g,b) 
	Vec(double x_ = 0, double y_ = 0, double z_ = 0){ x = x_; y = y_; z = z_; }
	Vec operator+(const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }
	Vec operator-(const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }
	Vec operator*(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
	Vec operator*(double b) const { return Vec(x*b, y*b, z*b); }
	Vec operator/(const Vec &b) const { return Vec(x / b.x, y / b.y, z / b.z); }
	Vec mult(const Vec &b) const { return Vec(x*b.x, y*b.y, z*b.z); }
	Vec& norm(){ return *this = *this * (1 / sqrt(x*x + y*y + z*z)); }
	double dot(const Vec &b) const { return x*b.x + y*b.y + z*b.z; }
	double length() const { return sqrt(dot(*this)); }
	Vec operator%(const Vec&b)const{ return Vec(y*b.z - z*b.y, z*b.x - x*b.z, x*b.y - y*b.x); }// cross
};
struct Ray { Vec o, d; Ray(){ ; } Ray(Vec o_, Vec d_) : o(o_), d(d_) {} };

struct Radiance{ 
	Vec emit; 
	Vec transmit;
	Radiance(){ emit = Vec(0, 0, 0); transmit = Vec(1.0, 1.0, 1.0); }
	Radiance(Vec c, Vec a) :emit(c), transmit(a){}
	Vec apply(const Vec s) const { return emit + s * transmit; }
};

