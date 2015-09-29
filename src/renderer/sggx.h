
// a test implementation of SGGX

// references:
// - The SGGX Microflake Distribution
// https://drive.google.com/file/d/0BzvWIdpUpRx_dXJIMk9rdEdrd00/view?usp=sharing
// - supplemental material
// https://drive.google.com/file/d/0BzvWIdpUpRx_djVyMG9jMnltdTg/view?usp=sharing
#pragma once

#include "renderer_math.h"


struct SGGX{
	double xx, yy, zz, xy, xz, yz;
	SGGX(double xx_, double yy_, double zz_, double xy_, double xz_, double yz_) :
		xx(xx_), yy(yy_), zz(zz_), xy(xy_), xz(xz_), yz(yz_) {}

	static SGGX interpolate(const SGGX& s0, const SGGX& s1){
		return SGGX(0.5*(s0.xx + s1.xx), 0.5*(s0.yy + s1.yy), 0.5*(s0.zz + s1.zz), 0.5*(s0.xy + s1.xy), 0.5*(s0.xz + s1.xz), 0.5*(s0.yz + s1.yz));
	}
};

struct SGGX_data{
	uint8_t sigma_x, sigma_y, sigma_z;
	int8_t rxy, rxz, ryz;
	SGGX_data(uint8_t sigx, uint8_t sigy, uint8_t sigz, int8_t xy, int8_t xz, int8_t yz) :
		sigma_x(sigx), sigma_y(sigy), sigma_z(sigz), rxy(xy), rxz(xz), ryz(yz) {}
	SGGX import()const{ 
		double x = (1.0 / 255.0) * (double)sigma_x;
		double y = (1.0 / 255.0) * (double)sigma_y;
		double z = (1.0 / 255.0) * (double)sigma_z;
		double xy = ((double)rxy) * ((rxy < 0.0) ? (1.0 / 128.0) : (1.0 / 127.0));
		double xz = ((double)rxz) * ((rxz < 0.0) ? (1.0 / 128.0) : (1.0 / 127.0));
		double yz = ((double)ryz) * ((ryz < 0.0) ? (1.0 / 128.0) : (1.0 / 127.0));
		return SGGX(x*x, y*y, z*z, xy*x*y, xz*x*z, yz*y*z);
	}
};
namespace sggx{
// build orthonormal basis (Building an Orthonormal Basis from a 3D Unit Vector Without Normalization, [Frisvad2012])
static inline void buildOrthonormalBasis(Vec& omega_1, Vec& omega_2, const Vec& omega_3)
{
	if (omega_3.z < -0.9999999)
	{
		omega_1 = Vec(0.0, -1.0, 0.0);
		omega_2 = Vec(-1.0, 0.0, 0.0);
	}
	else {
		const double a = 1.0 / (1.0 + omega_3.z);
		const double b = -omega_3.x * omega_3.y * a;
		omega_1 = Vec(1.0 - omega_3.x * omega_3.x * a, b, -omega_3.x);
		omega_2 = Vec(b, 1.0 - omega_3.y * omega_3.y * a, -omega_3.y);
	}
}

static inline Vec sample_VNDF(const Vec wi, const SGGX &S, Random &rnd)
{
	// generate sample (u, v, w)
	const double r = sqrt(rnd.next01());
	const double phi = 2.0*M_PI*rnd.next01();
	const double u = r*cos(phi);
	const double v = r*sin(phi);
	const double w = sqrt(1.0 - u*u - v*v);
	// build orthonormal basis
	Vec wk, wj;
	buildOrthonormalBasis(wk, wj, wi);
	// project S in this basis
	const double S_kk = wk.x*wk.x*S.xx + wk.y*wk.y*S.yy + wk.z*wk.z*S.zz
		+ 2.0 * (wk.x*wk.y*S.xy + wk.x*wk.z*S.xz + wk.y*wk.z*S.yz);
	const double S_jj = wj.x*wj.x*S.xx + wj.y*wj.y*S.yy + wj.z*wj.z*S.zz
		+ 2.0 * (wj.x*wj.y*S.xy + wj.x*wj.z*S.xz + wj.y*wj.z*S.yz);
	const double S_ii = wi.x*wi.x*S.xx + wi.y*wi.y*S.yy + wi.z*wi.z*S.zz
		+ 2.0 * (wi.x*wi.y*S.xy + wi.x*wi.z*S.xz + wi.y*wi.z*S.yz);
	const double S_kj = wk.x*wj.x*S.xx + wk.y*wj.y*S.yy + wk.z*wj.z*S.zz
		+ (wk.x*wj.y + wk.y*wj.x)*S.xy
		+ (wk.x*wj.z + wk.z*wj.x)*S.xz
		+ (wk.y*wj.z + wk.z*wj.y)*S.yz;
	const double S_ki = wk.x*wi.x*S.xx + wk.y*wi.y*S.yy + wk.z*wi.z*S.zz
		+ (wk.x*wi.y + wk.y*wi.x)*S.xy + (wk.x*wi.z + wk.z*wi.x)*S.xz + (wk.y*wi.z + wk.z*wi.y)*S.yz;
	const double S_ji = wj.x*wi.x*S.xx + wj.y*wi.y*S.yy + wj.z*wi.z*S.zz
		+ (wj.x*wi.y + wj.y*wi.x)*S.xy
		+ (wj.x*wi.z + wj.z*wi.x)*S.xz
		+ (wj.y*wi.z + wj.z*wi.y)*S.yz;
	// compute normal
	double sqrtDetSkji = sqrt(abs(S_kk*S_jj*S_ii - S_kj*S_kj*S_ii - S_ki*S_ki*S_jj - S_ji*S_ji*S_kk + 2.0*S_kj*S_ki*S_ji));
	double inv_sqrtS_ii = 1.0 / sqrt(S_ii);
	double tmp = sqrt(S_jj*S_ii - S_ji*S_ji);
	Vec Mk(sqrtDetSkji / tmp, 0.0, 0.0);
	Vec Mj(-inv_sqrtS_ii*(S_ki*S_ji - S_kj*S_ii) / tmp, inv_sqrtS_ii*tmp, 0);
	Vec Mi(inv_sqrtS_ii*S_ki, inv_sqrtS_ii*S_ji, inv_sqrtS_ii*S_ii);
	Vec wm_kji = (Mk * u + Mj * v + Mi * w).norm();
	// rotate back to world basis
	return wk * wm_kji.x + wj * wm_kji.y + wi * wm_kji.z;
}

inline double sigma(Vec wi, const SGGX &S)
{
	const double sigma_squared = wi.x*wi.x*S.xx + wi.y*wi.y*S.yy + wi.z*wi.z*S.zz
		+ 2.0 * (wi.x*wi.y*S.xy + wi.x*wi.z*S.xz + wi.y*wi.z*S.yz);
	return (0.0 < sigma_squared) ? sqrt(sigma_squared) : 0.0; // conditional to avoid numerical errors
}

static inline double D(Vec wm, const SGGX &S)
{
	const double detS =
		S.xx*S.yy*S.zz - S.xx*S.yz*S.yz - S.yy*S.xz*S.xz - S.zz*S.xy*S.xy + 2.0*S.xy*S.xz*S.yz;
	const double den = wm.x*wm.x*(S.yy*S.zz - S.yz*S.yz) + wm.y*wm.y*(S.xx*S.zz - S.xz*S.xz) + wm.z*wm.z*(S.xx*S.yy - S.xy*S.xy)
		+ 2.0*(wm.x*wm.y*(S.xz*S.yz - S.zz*S.xy) + wm.x*wm.z*(S.xy*S.yz - S.yy*S.xz) + wm.y*wm.z*(S.xy*S.xz - S.xx*S.yz));
	const double D = pow(abs(detS), 1.5) / (M_PI*den*den);
	return D;
}

static inline double eval_specular(const Vec in, const Vec out, const SGGX &S, Random &rnd)
{
	Vec wh = (in + out).norm();
	return 0.25 * D(wh, S) / sigma(in, S);
}

static inline Vec sample_specular(const Vec in, const SGGX &S, Random &rnd)
{
	// sample VNDF
	const Vec wm = sample_VNDF(in, S, rnd);

	// specular reflection
	const Vec out = in * (-1.0) + wm * 2.0 * in.dot(wm);
	return out;
}

static inline double eval_diffuse(const Vec in, const Vec out, const SGGX &S, Random &rnd)
{
	// sample VNDF
	const Vec wm = sample_VNDF(in, S, rnd);
	// eval diffuse
	double om = out.dot(wm);
	return (1.0 / M_PI) * (0.0 < om) ? om : 0.0;
}

static inline Vec sample_diffues(const Vec in, const SGGX &S, Random &rnd)
{
	// sample VNDF
	const Vec wm = sample_VNDF(in, S, rnd);

	Vec w1, w2;
	buildOrthonormalBasis(w1, w2, wm);

	double r1 = 2.0 * rnd.next01() - 1.0;
	double r2 = 2.0 * rnd.next01() - 1.0;

	// concentric map code from
	// http://psgraphics.blogspot.ch/2011/01/improved-code-for-concentric-map.html
	double phi, r;
	if (r1 * r1 <= DBL_EPSILON && r2 * r2 <= DBL_EPSILON) {
		r = phi = 0.0;
	}
	else if (r1*r1 > r2*r2) {
		r = r1;
		phi = (0.25 * M_PI) * (r2 / r1);
	}
	else {
		r = r2;
		phi = (0.5 * M_PI) - (r1 / r2) * (0.25 * M_PI);
	}
	double x = r * cos(phi);
	double y = r * sin(phi);
	double z = sqrt(1.0 - x * x - y * y);
	Vec out = w1 * x + w2 * y + wm * x;
	return out;
}
}
