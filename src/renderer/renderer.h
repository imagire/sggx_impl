
#pragma once


#include <climits>
#include <cmath> // exp, sqrt



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
	Vec operator%(Vec&b){ return Vec(y*b.z - z*b.y, z*b.x - x*b.z, x*b.y - y*b.x); }// cross
};
struct Ray { Vec o, d; Ray(){ ; } Ray(Vec o_, Vec d_) : o(o_), d(d_) {} };

struct Radiance{ 
	Vec emit; 
	double dampling; 
	Radiance(){ emit = Vec(0, 0, 0); dampling = 1.0; }
	Radiance(Vec c, double a) :emit(c),dampling(a){}
	Vec apply(const Vec s) const { return emit + s * dampling; }
};

enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance() 
struct Sphere {
	double rad;       // radius 
	Vec p, e, c;      // position, emission, color 
	Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive) 
	Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_) :
		rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}
	double intersect(const Ray &r) const { // returns distance, 0 if nohit 
		Vec op = p - r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0 
		double t, eps = 1e-4, b = op.dot(r.d), det = b*b - op.dot(op) + rad*rad;
		if (det<0) return 0; else det = sqrt(det);
		return (t = b - det)>eps ? t : ((t = b + det)>eps ? t : 0);
	}
};
Sphere spheres[] = {//Scene: radius, position, emission, color, material 
	Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.75, .25, .25), DIFF),//Left 
	Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, .75), DIFF),//Rght 
	Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(.75, .75, .75), DIFF),//Back 
	Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(.75, .75, .75), DIFF),//Frnt 
	Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(.75, .75, .75), DIFF),//Botm 
	Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.75, .75, .75), DIFF),//Top 
	Sphere(16.5, Vec(27, 16.5, 47), Vec(), Vec(1, 1, 1)*.999, SPEC),//Mirr 
	Sphere(16.5, Vec(73, 16.5, 78), Vec(), Vec(1, 1, 1)*.999, REFR),//Glas 
	Sphere(600, Vec(50, 681.6 - .27, 81.6), Vec(10, 7.3,4.7), Vec(), DIFF) //Lite 
};;
inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; }
inline int toInt(double x){ return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }
inline int intersect(const Ray &r, double &t){
	int id = -1;
	double n = sizeof(spheres) / sizeof(Sphere), d, inf = t = 1e20;
	for (int i = int(n); i--;) if ((d = spheres[i].intersect(r)) && d<t){ t = d; id = i; }
	return id;
}

//#define VOLUME_SIZE2
//#define VOLUME_SIZE4
#define VOLUME_SIZE8

#ifdef VOLUME_SIZE2
	#define VOLUME_X 2
	#define VOLUME_Y 2
	#define VOLUME_Z 2
	#define MASTER_DENSITY 0.02

uint32_t volume_color[VOLUME_Z][VOLUME_Y][VOLUME_X] =
{
	{
//		{ 0xffffffff, 0xffffffff },
//		{ 0xffffffff, 0xffffffff },
		{ 0x00ffffff, 0xffffffff },
		{ 0xffffffff, 0x00ffffff },
	},
	{
//		{ 0xffffffff, 0xffffffff },
//		{ 0xffffffff, 0xffffffff },
		{ 0x00ffffff, 0xffffffff },
		{ 0xffffffff, 0x00ffffff },
	},
};

#endif // VOLUME_SIZE2
#ifdef VOLUME_SIZE4
	#define VOLUME_X 4
	#define VOLUME_Y 4
	#define VOLUME_Z 4
	#define MASTER_DENSITY 0.02
//#define MASTER_DENSITY 0.02

	uint32_t volume_color[VOLUME_Z][VOLUME_Y][VOLUME_X] =
	{
		{
			{ 0x000000ff, 0xff0000ff, 0x000000ff, 0xff0000ff },
			{ 0xff0000ff, 0x000000ff, 0xff0000ff, 0x000000ff },
			{ 0x000000ff, 0xff0000ff, 0x000000ff, 0xff0000ff },
			{ 0xff0000ff, 0x000000ff, 0xff0000ff, 0x000000ff },
		},
		{
			{ 0xff00ff00, 0x0000ff00, 0xff00ff00, 0x0000ff00 },
			{ 0x0000ff00, 0xff00ff00, 0x0000ff00, 0xff00ff00 },
			{ 0xff00ff00, 0x0000ff00, 0xff00ff00, 0x0000ff00 },
			{ 0x0000ff00, 0xff00ff00, 0x0000ff00, 0xff00ff00 },
		},
		{
			{ 0x00ff0000, 0xffff0000, 0x00ff0000, 0xffff0000 },
			{ 0xffff0000, 0x00ff0000, 0xffff0000, 0x00ff0000 },
			{ 0x00ff0000, 0xffff0000, 0x00ff0000, 0xffff0000 },
			{ 0xffff0000, 0x00ff0000, 0xffff0000, 0x00ff0000 },
		},
		{
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
		},
	};

#endif // VOLUME_SIZE4
#ifdef VOLUME_SIZE8
	#define VOLUME_X 8
	#define VOLUME_Y 8
	#define VOLUME_Z 8
	#define MASTER_DENSITY 0.01

	uint32_t volume_color[VOLUME_Z][VOLUME_Y][VOLUME_X] =
	{
		{
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
		},
		{
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
		},
		{
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
		},
		{
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
		},
		{
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
		},
		{
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
		},
		{
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
		},
		{
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
			{ 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff },
			{ 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff, 0x00ffffff, 0xffffffff },
		},
	};

#endif // VOLUME_SIZE8

Vec cube_pos = Vec(20, 10.8, 55.6);
Vec cube_sca = Vec(60, 60, 120);
//Vec cube_pos = Vec(20, 10.8, 55.6);
//Vec cube_sca = Vec(60, 60, 1);
//#define DAMPING_CONST 0.01
//#define DAMPING_CONST 0.01

static inline Vec world_2_cube(const Vec in)
{
	return (in - cube_pos) / cube_sca;
}

static inline Vec cube_2_world(const Vec in)
{
	return in * cube_sca + cube_pos;
}

static inline Vec cube_2_world_SR(const Vec in)
{
	return in * cube_sca;
}

static inline bool AABB_test(const Vec in, const Vec out)
{
	Vec v_min = Vec(in.x < out.x ? in.x : out.x, in.y < out.y ? in.y : out.y, in.z < out.z ? in.z : out.z);
	Vec v_max = Vec(in.x < out.x ? out.x : in.x, in.y < out.y ? out.y : in.y, in.z < out.z ? out.z : in.z);

	return 0.0 < v_max.x && 0.0 <= v_max.y && 0.0 < v_max.z
		&& v_min.x < 1.0 && v_min.y < 1.0 && v_min.z < 1.0;
}

// in から in + diffで(0, 0, 0) - (1, 1, 1)を突き抜けるか判定
static inline bool intersect(const Vec in, const Vec diff, Ray &o)
{
	if (!AABB_test(in, in + diff)) return false;// 簡易テスト

	// 突き抜けた場合は、交差領域の始点と終点をo0, o1に入れてtrueで返る
	int state = 0;// 0:外側、1:internal 2: external again

	if (0.0 <= in.x && in.x < 1.0 &&
		0.0 <= in.y && in.y < 1.0 &&
		0.0 <= in.z && in.z < 1.0)
	{
		state = 1;
		o.o = in;
	}

	double t_min = 0.0;
	for (;state < 2; state++){
		double t_max = 1000000.0;
		Vec x_max;
		// x == 0
		double t = -in.x / diff.x;
		if (t_min < t && t < t_max){
			Vec x = in + diff * t;
			if (0.0 <= x.y && x.y < 1.0 &&
				0.0 <= x.z && x.z < 1.0)
			{
				t_max = t;
				x_max = x;
			}
		}
		// x == 1
		t = (1.0 - in.x) / diff.x;
		if (t_min < t && t < t_max){
			Vec x = in + diff * t;
			if (0.0 <= x.y && x.y < 1.0 &&
				0.0 <= x.z && x.z < 1.0)
			{
				t_max = t;
				x_max = x;
			}
		}
		// y == 0
		t = -in.y / diff.y;
		if (t_min < t && t < t_max){
			Vec x = in + diff * t;
			if (0.0 <= x.x && x.x < 1.0 &&
				0.0 <= x.z && x.z < 1.0)
			{
				t_max = t;
				x_max = x;
			}
		}
		// y == 1
		t = (1.0 - in.y) / diff.y;
		if (t_min < t && t < t_max){
			Vec x = in + diff * t;
			if (0.0 <= x.x && x.x < 1.0 &&
				0.0 <= x.z && x.z < 1.0)
			{
				t_max = t;
				x_max = x;
			}
		}
		// z == 0
		t = -in.z / diff.z;
		if (t_min < t && t < t_max){
			Vec x = in + diff * t;
			if (0.0 <= x.x && x.x < 1.0 &&
				0.0 <= x.y && x.y < 1.0)
			{
				t_max = t;
				x_max = x;
			}
		}
		// z == 1
		t = (1.0 - in.z) / diff.z;
		if (t_min < t && t < t_max){
			Vec x = in + diff * t;
			if (0.0 <= x.x && x.x < 1.0 &&
				0.0 <= x.y && x.y < 1.0)
			{
				t_max = t;
				x_max = x;
			}
		}

		if (0 == state){
			if (1.0 < t_max){
				// 交点がなければ交差しない
				return false;
			}
			else{
				o.o = x_max;
				t_min = t_max;
			}
		}
		else{
			if (1.0 < t_max){
				// 入って出なかったら最後の点が終点
				o.d = in + diff - o.o;
			}
			else{
				o.d = x_max - o.o;
			}
		}
	}

	return true;
}

static inline Radiance computep_participating_radiance(const int *id, const Vec i, const Vec o, const Radiance incoming)
{
	// ワールド空間に変換
	const Vec iw = cube_2_world(i);
	const Vec ow = cube_2_world(o);

	uint32_t voxel_data = volume_color[id[2]][id[1]][id[0]];
voxel_data |= 0xffffff;
	Vec c = Vec((1.0 / 255.0)*(double)(voxel_data & 0xff), (1.0 / 255.0)*(double)(voxel_data >> 8 & 0xff), (1.0 / 255.0)*(double)(voxel_data >> 16 & 0xff));
	double density = MASTER_DENSITY * (1.0 / 255.0)*(double)(voxel_data >> 24);

	double l = (ow-iw).length();
	double damping = l * density;
	damping = (1.0 < damping) ? 1.0 : damping;
	double damping_inv = 1.0 - damping;

	return Radiance(incoming.emit + c *incoming.dampling * damping, incoming.dampling * damping_inv);
}

Radiance ParticipatingEffect(Vec out, Vec in, Radiance backword)
{
	Ray ray;

	Vec x0 = world_2_cube(in);
	Vec x1 = world_2_cube(out);

	// x0-x1が(0,0,0)-(1,1,1)を突き抜けるか判定
	Vec d = x1 - x0;
	bool sign_x = 0.0 < d.x;
	bool sign_y = 0.0 < d.y;
	bool sign_z = 0.0 < d.z;

	Ray o;
	if (!intersect(x0, d, o)) return backword;// ボリュームと交差しませんでした。

	int id[3] = { 
		(int)((double)VOLUME_X * o.o.x + 0.00000 * d.x),
		(int)((double)VOLUME_Y * o.o.y + 0.00000 * d.y),
		(int)((double)VOLUME_Z * o.o.z + 0.00000 * d.z) };// エッジで判定するときわどい可能性があるので、少し進めて判定
	// todo: assert でひっかけてはみたい
	if (id[0] < 0) id[0] = 0;
	if (id[1] < 0) id[1] = 0;
	if (id[2] < 0) id[2] = 0;
	if (VOLUME_X <= id[0]) id[0] = VOLUME_X - 1;
	if (VOLUME_Y <= id[1]) id[1] = VOLUME_Y - 1;
	if (VOLUME_Z <= id[2]) id[2] = VOLUME_Z - 1;

	double t = 0.0;
	Vec x = o.o;
	do{
		double tx, ty, tz;
		if (sign_x){
			tx = ((id[0] + 1) / (double)VOLUME_X - o.o.x) / o.d.x;
		}
		else{
			tx = ((id[0]    ) / (double)VOLUME_X - o.o.x) / o.d.x;
		}
		if (sign_y){
			ty = ((id[1] + 1) / (double)VOLUME_Y - o.o.y) / o.d.y;
		}
		else{
			ty = ((id[1]    ) / (double)VOLUME_Y - o.o.y) / o.d.y;
		}
		if (sign_z){
			tz = ((id[2] + 1) / (double)VOLUME_Z - o.o.z) / o.d.z;
		}
		else{
			tz = ((id[2]    ) / (double)VOLUME_Z - o.o.z) / o.d.z;
		}
		if (tx < ty && tx < tz){
			t = (1.0 < tx) ? 1.0 : tx;
			Vec next = o.o + o.d * t;
			backword = computep_participating_radiance(id, x, next, backword);// id[0],id[1],id[2]でラディアンスの計算
			x = next;
			id[0] += (sign_x) ? 1 : -1;
			if (id[0] < 0) break;
			if (VOLUME_X <= id[0]) break;
		}
		else if (ty < tz && ty < tx){
			t = (1.0 < ty) ? 1.0 : ty;
			Vec next = o.o + o.d * t;
			backword = computep_participating_radiance(id, x, next, backword);// id[0],id[1],id[2]でラディアンスの計算
			x = next;
			id[1] += (sign_y) ? 1 : -1;
			if (id[1] < 0) break;
			if (VOLUME_Y <= id[1]) break;
		}
		else{
			t = (1.0 < tz) ? 1.0 : tz;
			Vec next = o.o + o.d * t;
			backword = computep_participating_radiance(id, x, next, backword);// id[0],id[1],id[2]でラディアンスの計算
			x = next;
			id[2] += (sign_z) ? 1 : -1;
			if (id[2] < 0) break;
			if (VOLUME_Z <= id[2]) break;
		}
	} while (t < 1.0);

	return backword;
}

Vec radiance(const Ray &r, int depth, const Radiance backword, Random &rnd){
	double t;                               // distance to intersection 

	// オブジェクトに衝突しなかった。
	int id = intersect(r, t); // id of intersected object 
	if (id < 0) return backword.emit; // if miss, return black 

	const Sphere &obj = spheres[id];        // the hit object 
	Vec x = r.o + r.d*t, // 衝突点
		n = (x - obj.p).norm(), // 衝突点の法線
		nl = n.dot(r.d)<0 ? n : n*-1, // レイの入射面からみた法線
		f = obj.c;// 物体の色

	Radiance pe = ParticipatingEffect(x, r.o, backword);

	double p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // 最大の強さの成分の色

	if (5<++depth){ // 多数回反射したとき
		if (rnd.next01() < p) 
			f = f * (1 / p);// 確率的に打ち切る色を取り込む
		else 
			return pe.apply(obj.e); // 打ち切る
	}

	if (100<depth) return pe.apply(obj.e); // 多すぎるとそこで打ち切り

	// 拡散反射
	if (obj.refl == DIFF){
		double r1 = 2 * M_PI*rnd.next01(), r2 = rnd.next01(), r2s = sqrt(r2);
		Vec w = nl, 
			u = ((fabs(w.x)>.1 ? Vec(0, 1) : Vec(1)) % w).norm(), // nl に関する従法線
			v = w % u;// nlに関する接ベクトル
		Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1 - r2)).norm();// 新しい拡散方向を定義
		return pe.apply(obj.e + f.mult(radiance(Ray(x, d), depth, pe, rnd)));
	}
	// 鏡面反射
	else if (obj.refl == SPEC){
		return pe.apply(obj.e + f.mult(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth, pe, rnd)));
	}
	// 屈折
	Ray reflRay(x, r.d - n * 2 * n.dot(r.d));     // 反射ベクトル
	bool into = n.dot(nl)>0;                // 入射していくならtrue
	double nc = 1, // 空気の屈折率
		nt = 1.5, // 透明物体の屈折率
		nnt = into ? nc / nt : nt / nc, // 入射前と入射後の屈折率の比（変異の向き、大きさに関与）
		ddn = r.d.dot(nl), // 屈折のための基本の変位
		cos2t;
	
	// 全反射の場合
	if ((cos2t = 1 - nnt*nnt*(1 - ddn*ddn)) < 0){
		if (into){
			return pe.apply(obj.e + f.mult(radiance(reflRay, depth, pe, rnd)));
		}
		else{
			return backword.apply(obj.e + f.mult(radiance(reflRay, depth, backword, rnd)));
		}
	}

	Vec tdir = (r.d*nnt - n*((into ? 1 : -1)*(ddn*nnt + sqrt(cos2t)))).norm(); // 屈折方向

	double a = nt - nc, 
		b = nt + nc, 
		R0 = a*a / (b*b), 
		c = 1 - (into ? -ddn : tdir.dot(n));
	double Re = R0 + (1 - R0)*c*c*c*c*c, // フレネル効果？
		Tr = 1 - Re, 
		P = .25 + .5*Re, 
		RP = Re / P, 
		TP = Tr / (1 - P);
	const Radiance &rr = into ? pe : backword;
	Vec out = obj.e + f.mult(depth>2 ? (rnd.next01() < P ?   // Russian roulette 
		radiance(reflRay, depth, rr, rnd)*RP : // 反射成分
		radiance(Ray(x, tdir), depth, rr, rnd)*TP) : // 屈折成分
		radiance(reflRay, depth, rr, rnd)*Re + radiance(Ray(x, tdir), depth, rr, rnd)*Tr);// 浅いところでは両方計算
	return pe.apply(out);
}

extern std::mutex mtx;

class renderer{
public:
	renderer(){}
	~renderer(){}
	inline void run(int w, int h, unsigned char *dst, bool *p_is_finished)
	{
		// global settings
		const int iteration_max = 256;
		Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm()); // cam pos, dir 
		Vec cx = Vec(w*.5135 / h), cy = (cx%cam.d).norm()*.5135;
		Vec *r = new Vec[4*w*h];// サブピクセルまで含めたフレームバッファ
		for (int i = 0; i < 4 * w*h; i++) r[i] = Vec(0, 0, 0);

		for (int l = 0; l < iteration_max; l++)
		{
			// ray tracing!!!
			const int samps = 8; // # samples iteration_max*samps がカメラからレイを飛ばす数になる
#pragma omp parallel for schedule(dynamic, 1)       // OpenMP 
			for (int y = 0; y < h; y++){                       // Loop over image rows 
				fprintf(stderr, "\r[iteration %d]: Rendering (%d spp) %5.2f%%", l, samps * 4, 100.*y / (h - 1));
				Random rnd(y + l * w * h + 1);
				for (unsigned short x = 0; x < w; x++){   // Loop cols 
					for (int sy = 0, i = (h - y - 1)*w + x; sy < 2; sy++){     // 2x2 subpixel rows 
						for (int sx = 0; sx < 2; sx++){        // 2x2 subpixel cols 
							Vec &ri = r[4 * i + 2 * sy + sx];
							for (int s = 0; s < samps; s++){
								double r1 = 2 * rnd.next01(), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
								double r2 = 2 * rnd.next01(), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
								Vec d = cx*(((sx + .5 + dx) / 2 + x) / w - .5) +
										cy*(((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
								ri = ri + radiance(Ray(cam.o + d * 140, d.norm()), 0, Radiance(), rnd)*(1. / samps);
							} // Camera rays are pushed ^^^^^ forward to start in interior 
						}
					}
				}
			}

			// 時間切れになっているようならおとなしく止める
			if (*p_is_finished) break;

			// 結果を画面イメージに書き出す (BMP用に上下反転)
			std::lock_guard<std::mutex> lock(mtx);
			const double tonemap_coeff = -0.25 / (double)(l+1);
			int addr = 0;
			for (int y = 0; y < h; y++){ 
				const Vec *v = &r[4 * (h - 1 - y) * w];
				for (unsigned short x = 0; x < w; x++){
					// サブピクセルの合成
					Vec c = v[0] + v[1] + v[2] + v[3];
					// tone mapping
					c.x = (c.x < 0.0) ? 0.0 : (1.0 - exp(c.x * tonemap_coeff));
					c.y = (c.y < 0.0) ? 0.0 : (1.0 - exp(c.y * tonemap_coeff));
					c.z = (c.z < 0.0) ? 0.0 : (1.0 - exp(c.z * tonemap_coeff));
					// gamma correct
					c.x = pow(c.x, 1.0 / 2.2);
					c.y = pow(c.y, 1.0 / 2.2);
					c.z = pow(c.z, 1.0 / 2.2);
					// 出力バッファに記録
					dst[addr + 2] = (unsigned char)(c.x * 255.99);
					dst[addr + 1] = (unsigned char)(c.y * 255.99);
					dst[addr + 0] = (unsigned char)(c.z * 255.99);
					addr += 3;
					v+=4;
				}
			}
		}
	}
};
