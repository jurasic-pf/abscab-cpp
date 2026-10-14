
#include "abscab.hh"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

namespace abscab {

/** C99 does not define M_PI */
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

/** vacuum magnetic permeability, divided by pi */
const double MU_0_BY_PI = MU_0 / M_PI;

/** vacuum magnetic permeability, divided by 2 pi */
const double MU_0_BY_2_PI = MU_0 / (2.0 * M_PI);

/** vacuum magnetic permeability, divided by 4 pi */
const double MU_0_BY_4_PI = MU_0 / (4.0 * M_PI);

/** half of pi */
const double PI_2 = M_PI / 2.0;

/** sqrt of machine epsilon */
const double SQRT_EPS = sqrt(DBL_EPSILON);

/**
 * Compute the complete elliptic integral introduced in
 * "Numerical Calculation of Elliptic Integrals and Elliptic Functions. III"
 * by R. Bulirsch in "Numerische Mathematik" 13, 305-315 (1969):
 * cel(k_c, p, a, b) =
 * \int_0^{\pi/2} \frac{a \cos^2{\varphi} + b \sin^2{\varphi}}
 *                     {  \cos^2{\varphi} + p \sin^2{\varphi}}
 *                \frac{\mathrm{d}\varphi}
 *                     {\sqrt{\cos^2{\varphi} + k_c^2 \sin^2{\varphi}}}
 * @param k_c parameter k_c of cel(); absolute value must not be 0
 * @param p   parameter p of cel()
 * @param a   parameter a of cel()
 * @param b   parameter b of cel()
 * @return the value of cel(k_c, p, a, b)
 */
double cel(double k_c, double p, double a, double b) {
	if (k_c == 0.0) {
		if (b != 0.0) {
			// when k_c is zero and b != 0, cel diverges (?)
			return HUGE_VAL;
		} else {
			k_c = SQRT_EPS*SQRT_EPS;
		}
	} else {
		k_c = fabs(k_c);
	}

	double m = 1.0; // \mu
	double e = k_c; // \nu * \mu
	// In the iterations, \nu_i is stored in k_c.

	double f, g;

	// initialization
	if (p > 0.0) {
		p = sqrt(p);
		b = b / p;
	} else { // p <= 0
		double q;

		f = k_c*k_c;        // f = kc^2 (re-used here; later f = a_i)
		q = 1.0 - f;        // 1 - kc^2
		g = 1.0 - p;
		f -= p;             // kc^2 - p
		q *= b-a*p;         // (1 - kc^2)*(b-a*p)
		p = sqrt(f/g);      // sqrt((kc^2 - p)/(1-p))              --> p0 done here
		a = (a-b)/g;        // (a-b)/(1-p)                         --> a0 done here
		b = -q/(g*g*p)+a*p; // -(1 - kc^2)*(b-a*p)/( (1-p)^2 * p ) --> b0 done here
	}

	// iteration until convergence
	while (1) {
		f = a;    // f = a_i
		a += b/p; // a_{i+1} <-- a_i + b_i/p_i
		g = e/p;  // g = (\nu_i*\mu_i) / p_i
		b += f*g; // b_i + a_i * (\nu_i*\mu_i) /p_i
		b += b;   // b_{i+1} <-- 2*( b_i + a_i * (\nu_i*\mu_i) /p_i )
		p += g;   // p_{i+1} <-- p_i + (\nu_i*\mu_i) / p_i
		g = m;    // g = mu_i
		m += k_c; // mu_{i+1} <-- mu_i + \nu_i
		if (fabs(g - k_c) > g * SQRT_EPS) { // |\mu_i - \nu_i| > \mu_i*EPS == |1 - \nu_i / \mu_i| > EPS, but more robust!
			// not converged yet...
			k_c = sqrt(e); // k_c = sqrt(\nu_i * \mu_i)
			k_c += k_c;    // \nu_{i+1} <-- 2*sqrt(\nu_i * \mu_i)
			e = k_c*m;     // (\nu * \mu)_{i+1} <-- \nu_{i+1} * mu_{i+1} (also update product explicitly)
		} else {
			break;
		}
	}

	// final approximation:
	// \pi/2 * (a * \mu + b)/(\mu * (\mu + p))
	return PI_2 * (a*m+b) / (m*(m+p));
}

/**
 * Add a single contribution to the compensated summation.
 * The compensated sum is obtained by summing the final values of s, cs and ccs
 * after this method has been called for all contributions.
 *
 * @param contribution contribution to add to the sum
 * @param compSum[3]: {s, cs, ccs}: target for output
 */
void compAdd(double contribution, double *compSum) {
	double s   = compSum[0];
	double cs  = compSum[1];

	double t = s + contribution;
	double c;
	if (fabs(s) >= fabs(contribution)) {
		c = (s - t) + contribution;
	} else {
		c = (contribution - t) + s;
	}
	compSum[0] = t;

	double t2 = cs + c;
	double cc;
	if (fabs(cs) >= fabs(c)) {
		cc = (cs - t2) + c;
	} else {
		cc = (c - t2) + cs;
	}
	compSum[1] = t2;
	compSum[2] += cc;
}

/////// A_z of straight wire segment

/**
 * Compute the normalized axial component of magnetic vector potential of straight wire segment,
 * evaluated along axis of wire segment (rho = 0).
 * This is a special case for points away from the wire ("far-field") for zP < -1 or zP >= 2.
 *
 * @param zP normalized axial coordinate of evaluation location; must not be in [0, 1] (on wire segment)
 * @return normalized axial component of magnetic vector potential
 */
double sws_A_z_ax_f(double zP) {
	return atanh(1 / (fabs(zP) + fabs(1 - zP)));
}

/**
 * Compute the normalized axial component of magnetic vector potential of straight wire segment,
 * evaluated along axis of wire segment (rhoP = 0).
 * This is a special case for points close to the wire ("near-field") for -1 <= zP < 2.
 *
 * @param zP normalized axial coordinate of evaluation location; must not be in [0, 1] (on wire segment)
 * @return normalized axial component of magnetic vector potential
 */
double sws_A_z_ax_n(double zP) {
	// Two negative signs must be able to cancel each other here!
	return copysign(1.0, zP) * log(zP / (zP - 1)) / 2;
}

/**
 * Compute the normalized axial component of magnetic vector potential of straight wire segment,
 * evaluated along axis of wire segment (rho = 0).
 *
 * @param zP normalized axial coordinate of evaluation location; must not be in [0, 1] (on wire segment)
 * @return normalized axial component of magnetic vector potential
 */
double sws_A_z_ax(double zP) {
	if (zP < -1 || zP >= 2) {
		return sws_A_z_ax_f(zP);
	} else {
		return sws_A_z_ax_n(zP);
	}
}

/**
 * Compute the normalized axial component of the magnetic vector potential of a straight wire segment,
 * evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
 * This is a special case for points away from the wire ("far-field") for rhoP > 1.
 *
 * @param rhoP normalized radial coordinate of evaluation location; must not be zero (on wire segment)
 * @return normalized axial component of magnetic vector potential
 */
double sws_A_z_rad_f(double rhoP) {
	return atanh(1 / (rhoP + hypot(rhoP, 1)));
}

/**
 * Compute the normalized axial component of the magnetic vector potential of a straight wire segment,
 * evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
 * This is a special case for points close to the wire ("near-field") for rhoP <= 1.
 *
 * @param rhoP normalized radial coordinate of evaluation location; must not be zero (on wire segment)
 * @return normalized axial component of magnetic vector potential
 */
double sws_A_z_rad_n(double rhoP) {
	double cat = 1 / hypot(rhoP, 1);  // cos(atan(...)  )
	double sat = sin(atan(rhoP) / 2); // sin(atan(...)/2)
	double rc = rhoP * cat;
	double num = rc + 1 + cat;
	double den = rc + 2 * sat * sat;
	return log(num / den) / 2;
}

/**
 * Compute the normalized axial component of the magnetic vector potential of a straight wire segment,
 * evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
 *
 * @param rhoP normalized radial coordinate of evaluation location; must not be zero (on wire segment)
 * @return normalized axial component of magnetic vector potential
 */
double sws_A_z_rad(double rhoP) {
	if (rhoP > 1) {
		return sws_A_z_rad_f(rhoP);
	} else {
		return sws_A_z_rad_n(rhoP);
	}
}

/**
 * Compute the normalized axial component of the magnetic vector potential of a straight wire segment.
 * This formulation is useful for points away from the wire ("far-field")
 * at rhoP >= 1 or zP <= -1 or zP > 2.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized axial component of magnetic vector potential
 */
double sws_A_z_f(double rhoP, double zP) {
	double r_i = hypot(rhoP, zP);
	double r_f = hypot(rhoP, 1 - zP);
	return atanh(1 / (r_i + r_f));
}

/**
 * Compute the normalized axial component of the magnetic vector potential of a straight wire segment.
 * This formulation is useful for points close to the wire ("near-field")
 * at rhoP < 1 and -1 < zP <= 2.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized axial component of magnetic vector potential
 */
double sws_A_z_n(double rhoP, double zP) {
	double omz = 1 - zP;

	double r_i = hypot(rhoP, zP);
	double r_f = hypot(rhoP, omz);

	double alpha = atan2(rhoP, zP);
	double sinAlphaHalf = sin(alpha / 2);

	double beta = atan2(rhoP, omz);
	double sinBetaHalf = sin(beta / 2);

	double Ri_zP    = r_i * sinAlphaHalf * sinAlphaHalf; // 0.5 * (r_i - z')
	double Rf_p_zM1 = r_f * sinBetaHalf  * sinBetaHalf;  // 0.5 * (r_f - (1 - z'))

	double n = Ri_zP + Rf_p_zM1;

	return (log(1 + n) - log(n)) / 2;
}

/////// B_phi of straight wire segment

/**
 * Compute the normalized tangential component of the magnetic field of a straight wire segment,
 * evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @return normalized tangential component of magnetic field
 */
double sws_B_phi_rad(double rhoP) {
	return 1 / (rhoP * hypot(rhoP, 1));
}

/**
 * Compute the normalized tangential component of the magnetic field of a straight wire segment.
 * This formulation is useful for points away from the wire ("far-field")
 * at rhoP >= 1 or zP <= 0 or zP >= 1 or rhoP/(1-zP) >= 1 or rhoP/zP >= 1.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized tangential component of magnetic field
 */
double sws_B_phi_f(double rhoP, double zP) {
	double omz = 1 - zP;

	double r_i = hypot(rhoP, zP);
	double r_f = hypot(rhoP, omz);

	double num = rhoP * (1/r_i + 1/r_f);
	double den = rhoP * rhoP - zP * omz + r_i * r_f;

	return num / den;
}

/**
 * Compute the normalized tangential component of the magnetic field of a straight wire segment.
 * This formulation is useful for points close to the wire ("near-field")
 * at rhoP < 1 and 0 < zP < 1 and rhoP/(1-zP) < 1 and rhoP/zP < 1.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized tangential component of magnetic field
 */
double sws_B_phi_n(double rhoP, double zP) {
	double omz = 1 - zP;

	double r_i = hypot(rhoP, zP);
	double r_f = hypot(rhoP, omz);

	double num = rhoP * (1/r_i + 1/r_f);

	double alpha = atan2(rhoP, zP);
	double sinAlphaHalf = sin(alpha / 2);

	double beta = atan2(rhoP, omz);
	double sinBetaHalf = sin(beta / 2);

	// r_f * sin^2(beta/2) + (1 - z') * sin^2(alpha/2)
	double rfb_omza = r_f * sinBetaHalf * sinBetaHalf + omz * sinAlphaHalf * sinAlphaHalf;

	//     r_i * r_f - z' * (1 - z')
	// =   r_i * r_f - r_i * (1 - z') + r_i * (1 - z') - z' * (1 - z')
	// =   r_i * r_f - r_i * r_f * cos(beta)
	//   + r_i * (1 - z') + (1 - z') * r_i * cos(alpha)
	// =   r_i *    r_f   * (1 - cos(beta))
	//   + r_i * (1 - z') * (1 - cos(alpha))
	// = 2 * r_i * [ r_f * sin^2(beta/2) + (1 - z') * sin^2(alpha/2) ]
	double den = rhoP * rhoP + 2 * r_i * rfb_omza;

	return num / den;
}

///// A_phi of circular wire loop

/**
 * Compute the normalized tangential component of the magnetic vector potential of a circular wire loop.
 * This formulation is useful for points away from the wire ("far-field")
 * at rhoP < 1/2 or rhoP > 2 or |zP| >= 1.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized tangential component of magnetic vector potential
 */
double cwl_A_phi_f(double rhoP, double zP) {
	double sqrt_kCSqNum = hypot(zP, 1 - rhoP);
	double sqrt_kCSqDen = hypot(zP, 1 + rhoP);

	double kC = sqrt_kCSqNum / sqrt_kCSqDen;
	double kSq = 4 * rhoP / (sqrt_kCSqDen * sqrt_kCSqDen);

	double kCp1 = 1 + kC;
	double arg1 = 2 * sqrt(kC) / kCp1;
	double arg2 = 2 / (kCp1 * kCp1 * kCp1);
	double C = cel(arg1, 1, 0, arg2);

	return kSq/sqrt_kCSqDen * C;
}

/**
 * Compute the normalized tangential component of the magnetic vector potential of a circular wire loop.
 * This formulation is useful for points close to the wire ("near-field")
 * at 1/2 <= rhoP <= 2 and |zP| < 1.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized tangential component of magnetic vector potential
 */
double cwl_A_phi_n(double rhoP, double zP) {
	double rhoP_m_1 = rhoP - 1;

	double n = zP / rhoP_m_1;
	double m = 1 + 2 / rhoP_m_1;

	double num = n * n + 1;
	double den = n * n + m * m;

	double kCSq = num / den;

	double prefac = 1 / (fabs(rhoP - 1) * sqrt(den));
	double celPart = cel(sqrt(kCSq), 1, -1, 1);
	return prefac * celPart;
}

/**
 * Compute the normalized tangential component of the magnetic vector potential of a circular wire loop.
 * This formulation is useful for points along rhoP=1 with |zP| < 1.
 *
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized tangential component of magnetic vector potential
 */
double cwl_A_phi_v(double zP) {
	double absZp = fabs(zP);

	// 1/k_c
	double kCInv = sqrt(4 + zP * zP) / absZp;

	return cel(kCInv, 1, 1, -1) / absZp;
}

//////// B_rho of circular wire loop

/**
 * Compute the normalized radial component of the magnetic field of a circular wire loop.
 * This formulation is useful for points away from the wire ("far-field")
 * at rhoP < 1/2 or rhoP > 2 or |zP| >= 1.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized radial component of magnetic field
 */
double cwl_B_rho_f(double rhoP, double zP) {
	double sqrt_kCSqNum = hypot(zP, 1 - rhoP);
	double sqrt_kCSqDen = hypot(zP, 1 + rhoP);

	double kCSqNum = sqrt_kCSqNum * sqrt_kCSqNum;
	double kCSqDen = sqrt_kCSqDen * sqrt_kCSqDen;

	double kCSq = kCSqNum / kCSqDen;
	double kC = sqrt(kCSq);

	double D = cel(kC, 1, 0, 1);

	double kCp1 = 1 + kC;
	double arg1 = 2 * sqrt(kC) / kCp1;
	double arg2 = 2 / (kCp1 * kCp1 * kCp1);
	double C = cel(arg1, 1, 0, arg2);

	double prefac = 4 * rhoP / (kCSqDen * sqrt_kCSqDen * kCSqNum);

	return prefac * zP * (D - C);
}

/**
 * Compute the normalized radial component of the magnetic field of a circular wire loop.
 * This formulation is useful for points close to the wire ("near-field")
 * at 1/2 <= rhoP <= 2 and |zP| < 1.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized radial component of magnetic field
 */
double cwl_B_rho_n(double rhoP, double zP) {
	double rhoP_m_1 = rhoP - 1;
	double rd2 = rhoP_m_1 * rhoP_m_1;

	double n = zP / rhoP_m_1;
	double m = 1 + 2 / rhoP_m_1;

	double sqrt_kCSqNum = hypot(n, 1);
	double sqrt_kCSqDen = hypot(n, m);

	double kCSqNum = sqrt_kCSqNum * sqrt_kCSqNum;
	double kCSqDen = sqrt_kCSqDen * sqrt_kCSqDen;

	double kC = sqrt_kCSqNum / sqrt_kCSqDen;

	double D = cel(kC, 1, 0, 1);

	double kCp1 = 1 + kC;
	double arg1 = 2 * sqrt(kC) / kCp1;
	double arg2 = 2 / (kCp1 * kCp1 * kCp1);
	double C = arg2 * cel(arg1, 1, 0, 1);

	// z' / |rho' - 1|^5
	double zP_rd5 = zP / (fabs(rhoP_m_1) * rd2 * rd2);

	double prefac = 4 * rhoP / (kCSqDen * sqrt_kCSqDen * kCSqNum);

	return prefac * zP_rd5 * (D - C);
}

/**
 * Compute the normalized radial component of the magnetic field of a circular wire loop.
 * This formulation is useful for points along rhoP=1 with |zP| < 1.
 *
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized radial component of magnetic field
 */
double cwl_B_rho_v(double zP) {
	double zPSq = zP * zP;

	double kCSq = 1 / (1 + 4 / zPSq);
	double kC = sqrt(kCSq);

	double K = cel(kC, 1, 1, 1);
	double E = cel(kC, 1, 1, kCSq);

	return copysign(kC / 2 * ((2 / zPSq + 1) * E - K), zP);
}

////// B_z of circular wire loop

/**
 * Compute the normalized vertical component of the magnetic field of a circular wire loop.
 * This formulation is useful for certain points away from the wire ("far-field")
 * at rhoP < 1/2 or (rhoP <= 2 and |zP| >= 1).
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized vertical component of magnetic field
 */
double cwl_B_z_f1(double rhoP, double zP) {
	double sqrt_kCSqNum = hypot(zP, 1 - rhoP);
	double sqrt_kCSqDen = hypot(zP, 1 + rhoP);

	double kC = sqrt_kCSqNum / sqrt_kCSqDen;

	double K = cel(kC, 1, 1, 1);
	double E = cel(kC, 1, 1, kC * kC);
	double D = cel(kC, 1, 0, 1);

	double prefac = 1 / (sqrt_kCSqDen * sqrt_kCSqNum * sqrt_kCSqNum);
	double comb = (E - 2 * K + 2 * D);

	return prefac * (E + rhoP * comb);
}

/**
 * Compute the normalized vertical component of the magnetic field of a circular wire loop.
 * This formulation is useful for certain other points away from the wire ("far-field")
 * at rhoP > 2.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized vertical component of magnetic field
 */
double cwl_B_z_f2(double rhoP, double zP) {
	double sqrt_kCSqNum = hypot(zP, 1 - rhoP);
	double sqrt_kCSqDen = hypot(zP, 1 + rhoP);

	double kC = sqrt_kCSqNum / sqrt_kCSqDen;
	double kCSq = kC * kC;

	double zPSqP1 = zP * zP + 1;
	double rhoPSq = rhoP * rhoP;
	double t1 = zPSqP1 / rhoPSq + 1;
	double t2 = 2 / rhoP;

	// a is sqrt_kCSqDen normalized to rho'^2
	// b is sqrt_kCSqNum normalized to rho'^2
	// a == (z'^2 + (1 + rho')^2) / rho'^2 = (z'^2 + 1)/rho'^2 + 1  +  2/rho'
	// b == (z'^2 + (1 - rho')^2) / rho'^2 = (z'^2 + 1)/rho'^2 + 1  -  2/rho'
	double a = t1 + t2;
	double b = t1 - t2;

	// 1/prefac = sqrt( z'^2 + (1 + rho')^2)           * (z'^2 + (1 - rho')^2)
	//          = sqrt((z'^2 + (1 + rho')^2) / rho'^2) * (z'^2 + (1 - rho')^2) / rho'^2 * rho'^3
	//          = sqrt(a)                              * b                              * rho'^3
	double prefac = 1 / (sqrt(a) * b * rhoPSq * rhoP);

	double cdScale = 1 + (2 + zPSqP1 / rhoP) / rhoP;

	double E = cel(kC, 1, 1, kCSq);
	double D = cel(kC, 1, 0, 1);

	double kCP1 = 1 + kC;
	double arg1 = 2 * sqrt(kC) / kCP1;
	double arg2 = 2 / (kCP1 * kCP1 * kCP1);
	double C = arg2 * cel(arg1, 1, 0, 1);

	// use C - D for (2 * D - E)/kSq
	return prefac * (E + 4 * (C - D) / cdScale);
}

/**
 * Compute the normalized vertical component of the magnetic field of a circular wire loop.
 * This formulation is useful for points close to the wire ("near-field")
 * at 1/2 <= rhoP <= 2, but not rhoP=1, and |zP| <= 1.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized vertical component of magnetic field
 */
double cwl_B_z_n(double rhoP, double zP) {
	double rp1 = rhoP - 1;

	double n = zP / rp1;
	double m = 1 + 2 / rp1;

	double sqrt_kCSqNum = hypot(n, 1);
	double sqrt_kCSqDen = hypot(n, m);

	double kCSqDen = sqrt_kCSqDen * sqrt_kCSqDen;

	double kC = sqrt_kCSqNum / sqrt_kCSqDen;

	double prefac = 1 / (fabs(rp1) * rp1 * rp1 * kCSqDen * sqrt_kCSqDen);

	return prefac * cel(kC, kC * kC, 1 + rhoP, 1 - rhoP);
}

/**
 * Compute the normalized vertical component of the magnetic field of a circular wire loop.
 * This formulation is useful for points along rhoP=1 with |zP| <= 1.
 *
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized vertical component of magnetic field
 */
double cwl_B_z_v(double zP) {
	double kCSq = zP * zP / (4 + zP * zP);
	double kC = sqrt(kCSq);

	double f = zP * zP + 4;
	double prefac = 1 / (f * sqrt(f));

	return prefac * cel(kC, kCSq, 2, 0);
}

// --------------------------------------------------

/**
 * Compute the normalized axial component of the magnetic vector potential of a straight wire segment.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized axial component of magnetic vector potential
 */
double straightWireSegment_A_z(double rhoP, double zP) {
	if (rhoP == 0.0) {
		if (zP < 0 || zP > 1) {
			return sws_A_z_ax(zP);
		} else {
			fprintf(stderr, "evaluation locations on the wire segment (rho'=%g z'=%g) are not allowed\n", rhoP, zP);
			return 0.0;
		}
	} else if (zP == 0.0 || zP == 1.0) {
		return sws_A_z_rad(rhoP);
	} else if (rhoP >= 1.0 || zP <= -1.0 || zP > 2.0) {
		return sws_A_z_f(rhoP, zP);
	} else {
		return sws_A_z_n(rhoP, zP);
	}
}

/**
 * Compute the normalized tangential component of the magnetic field of a straight wire segment.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized tangential component of magnetic field
 */
double straightWireSegment_B_phi(double rhoP, double zP) {
	if (rhoP == 0.0) {
		if (zP < 0 || zP > 1) {
			return 0.0;
		} else {
			fprintf(stderr, "evaluation locations on the wire segment (rho'=%g z'=%g) are not allowed\n", rhoP, zP);
			return 0.0;
		}
	} else if (zP == 0.0 || zP == 1.0) {
		return sws_B_phi_rad(rhoP);
	} else if (rhoP >= zP || rhoP >= 1 - zP || zP < 0.0 || zP > 1.0) {
		return sws_B_phi_f(rhoP, zP);
	} else {
		return sws_B_phi_n(rhoP, zP);
	}
}

/**
 * Geometric part of magnetic vector potential computation for circular wire
 * loop at rho'=1, z'=0 (normalized coordinates). This routine selects special
 * case routines to get the most accurate formulation for given evaluation
 * coordinates.
 *
 * @param rhoP normalized radial evaluation position
 * @param zP   normalized vertical evaluation position
 * @return A_phi: toroidal component of magnetic vector potential: geometric
 *         part (no mu0*I/pi factor included)
 */
double circularWireLoop_A_phi(double rhoP, double zP) {
	if (rhoP == 0.0) {
		return 0.0;
	} else if (rhoP < 0.5 || rhoP > 2.0 || fabs(zP) >= 1.0) {
		return cwl_A_phi_f(rhoP, zP);
	} else if (rhoP != 1.0) {
		return cwl_A_phi_n(rhoP, zP);
	} else {
		if (zP != 0) {
			return cwl_A_phi_v(zP);
		} else {
			fprintf(stderr, "evaluation at location of wire loop (rho' = 1, z' = 0) is not defined\n");
			return 0.0;
		}
	}
}

/**
 * Geometric part of radial magnetic field computation for circular wire loop at
 * rho'=1, z'=0 (normalized coordinates). This routine selects special case
 * routines to get the most accurate formulation for given evaluation
 * coordinates.
 *
 * @param rhoP normalized radial evaluation position
 * @param zP   normalized vertical evaluation position
 * @return B_rho: radial component of magnetic field: geometric part (no
 *         mu0*I/(pi*a) factor included)
 */
double circularWireLoop_B_rho(double rhoP, double zP) {
	if (rhoP == 0.0 || zP == 0.0) {
		if (rhoP != 1.0) {
			return 0.0;
		} else {
			fprintf(stderr, "evaluation at location of wire loop (rho' = 1, z' = 0) is not defined\n");
			return 0.0;
		}
	} else if (rhoP < 0.5 || rhoP > 2.0 || fabs(zP) >= 1.0) {
		return cwl_B_rho_f(rhoP, zP);
	} else if (rhoP != 1.0) {
		return cwl_B_rho_n(rhoP, zP);
	} else {
		return cwl_B_rho_v(zP);
	}
}

/**
 * Geometric part of vertical magnetic field computation for circular wire loop
 * at rho'=1, z'=0 (normalized coordinates). This routine selects special case
 * routines to get the most accurate formulation for given evaluation
 * coordinates.
 *
 * @param rhoP normalized radial evaluation position
 * @param zP   normalized vertical evaluation position
 * @return B_z: vertical component of magnetic field: geometric part (no
 *         mu0*I/(pi*a) factor included)
 */
double circularWireLoop_B_z(double rhoP, double zP) {
	if (rhoP < 0.5 || (rhoP <= 2 && fabs(zP) > 1)) {
		return cwl_B_z_f1(rhoP, zP);
	} else if (rhoP > 2) {
		return cwl_B_z_f2(rhoP, zP);
	} else if (rhoP != 1.0) {
		return cwl_B_z_n(rhoP, zP);
	} else {
		if (zP != 0) {
			return cwl_B_z_v(zP);
		} else {
			fprintf(stderr, "evaluation at location of wire loop (rho' = 1, z' = 0) is not defined\n");
			return 0.0;
		}
	}
}

// --------------------------------------------------

/**
 * Compute the magnetic vector potential of a circular wire loop.
 *
 * @param center  [3: x, y, z] origin of loop (in meters)
 * @param normal  [3: x, y, z] normal vector of loop (in meters); will be
 *                normalized internally
 * @param radius  radius of the wire loop (in meters)
 * @param current loop current (in A)
 * @param evalPos [3: x, y, z][nEvalPos] evaluation locations (in meters)
 * @param vectorPotential [3: A_x, A_y, A_z][nEvalPos] Cartesian components of the magnetic
 *         vector potential evaluated at the given locations (in Tm); has to be allocated on entry
 */
void vectorPotentialCircularFilament(double *center, double *normal, double radius,
		double current, int nEvalPos, double *evalPos, double *vectorPotential) {

	if (!isfinite(radius) || radius <= 0.0) {
		printf("radius must be finite and positive, but is %g\n", radius);
		return;
	}

	double aPrefactor = MU_0_BY_PI * current;

	// squared length of normal vector
	double nLen2 = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];

	if (nLen2 == 0.0) {
		printf("length of normal vector must not be zero");
		return;
	}

	// length of normal vector
	double nLen = sqrt(nLen2);

	// unit normal vector of wire loop
	double eX = normal[0] / nLen;
	double eY = normal[1] / nLen;
	double eZ = normal[2] / nLen;

	for (int idxEval = 0; idxEval < nEvalPos; ++idxEval) {

		// vector from center of wire loop to eval pos
		double r0x = evalPos[3 * idxEval + 0] - center[0];
		double r0y = evalPos[3 * idxEval + 1] - center[1];
		double r0z = evalPos[3 * idxEval + 2] - center[2];

		// z position along normal of wire loop
		double alignedZ = eX * r0x + eY * r0y + eZ * r0z;

		// normalized z component of evaluation location in coordinate system of wire loop
		double zP = alignedZ / radius;

		// r0 projected onto axis of wire loop
		double rParallelX = alignedZ * eX;
		double rParallelY = alignedZ * eY;
		double rParallelZ = alignedZ * eZ;

		// vector perpendicular to axis of wire loop, pointing at evaluation pos
		double rPerpX = r0x - rParallelX;
		double rPerpY = r0y - rParallelY;
		double rPerpZ = r0z - rParallelZ;

		// perpendicular distance squared between evalPos and axis of wire loop
		double alignedRSq = rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ;

		// prevent division-by-zero when computing radial unit vector
		// A_phi is zero anyway on-axis --> no contribution expected
		if (alignedRSq > 0.0) {

			// perpendicular distance between evalPos and axis of wire loop
			double alignedR = sqrt(alignedRSq);

			// unit vector in radial direction
			double eRX = rPerpX / alignedR;
			double eRY = rPerpY / alignedR;
			double eRZ = rPerpZ / alignedR;

			// normalized rho component of evaluation location in coordinate system of wire loop
			double rhoP = alignedR / radius;

			// compute tangential component of magnetic vector potential, including current and mu_0
			double aPhi = aPrefactor * circularWireLoop_A_phi(rhoP, zP);

			// compute cross product between e_z and e_rho to get e_phi
			double ePhiX = eRY * eZ - eRZ * eY;
			double ePhiY = eRZ * eX - eRX * eZ;
			double ePhiZ = eRX * eY - eRY * eX;

			// add contribution from wire loop to result
			vectorPotential[3 * idxEval + 0] += aPhi * ePhiX;
			vectorPotential[3 * idxEval + 1] += aPhi * ePhiY;
			vectorPotential[3 * idxEval + 2] += aPhi * ePhiZ;
		}
	}
}

/**
 * Compute the magnetic field of a circular wire loop.
 *
 * @param center  [3: x, y, z] origin of loop (in meters)
 * @param normal  [3: x, y, z] normal vector of loop (in meters); will be
 *                normalized internally
 * @param radius  radius of the wire loop (in meters)
 * @param current loop current (in A)
 * @param evalPos [3: x, y, z][nEvalPos] evaluation locations (in meters)
 * @param magneticField [3: B_x, B_y, B_z][nEvalPos] Cartesian components of the magnetic
 *         field evaluated at the given locations (in T); has to be allocated on entry
 */
void magneticFieldCircularFilament(double *center, double *normal, double radius,
		double current, int nEvalPos, double *evalPos, double *magneticField) {

	if (!isfinite(radius) || radius <= 0.0) {
		printf("radius must be finite and positive, but is %g\n", radius);
		return;
	}

	double bPrefactor = MU_0_BY_PI * current / radius;

	// squared length of normal vector
	double nLen2 = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];

	if (nLen2 == 0.0) {
		printf("length of normal vector must not be zero");
		return;
	}

	// length of normal vector
	double nLen = sqrt(nLen2);

	// unit normal vector of wire loop
	double eX = normal[0] / nLen;
	double eY = normal[1] / nLen;
	double eZ = normal[2] / nLen;

	for (int idxEval = 0; idxEval < nEvalPos; ++idxEval) {

		// vector from center of wire loop to eval pos
		double r0x = evalPos[3 * idxEval + 0] - center[0];
		double r0y = evalPos[3 * idxEval + 1] - center[1];
		double r0z = evalPos[3 * idxEval + 2] - center[2];

		// z position along normal of wire loop
		double alignedZ = eX * r0x + eY * r0y + eZ * r0z;

		// normalized z component of evaluation location in coordinate system of wire loop
		double zP = alignedZ / radius;

		// r0 projected onto axis of wire loop
		double rParallelX = alignedZ * eX;
		double rParallelY = alignedZ * eY;
		double rParallelZ = alignedZ * eZ;

		// vector perpendicular to axis of wire loop, pointing at evaluation pos
		double rPerpX = r0x - rParallelX;
		double rPerpY = r0y - rParallelY;
		double rPerpZ = r0z - rParallelZ;

		// perpendicular distance squared between evalPos and axis of wire loop
		double alignedRSq = rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ;

		double rhoP;
		if (alignedRSq > 0.0) {
			// radial unit vector is only defined if evaluation pos is off-axis

			// perpendicular distance between evalPos and axis of wire loop
			double alignedR = sqrt(alignedRSq);

			// unit vector in radial direction
			double eRX = rPerpX / alignedR;
			double eRY = rPerpY / alignedR;
			double eRZ = rPerpZ / alignedR;

			// normalized rho component of evaluation location in coordinate system of wire loop
			rhoP = alignedR / radius;

			// compute radial component of normalized magnetic field
			// and scale by current and mu_0
			double bRho = bPrefactor * circularWireLoop_B_rho(rhoP, zP);

			// add contribution from B_rho of wire loop to result
			magneticField[3 * idxEval + 0] += bRho * eRX;
			magneticField[3 * idxEval + 1] += bRho * eRY;
			magneticField[3 * idxEval + 2] += bRho * eRZ;
		} else {
			rhoP = 0.0;
		}

		// compute vertical component of normalized magnetic field
		// and scale by current and mu_0
		double bZ = bPrefactor * circularWireLoop_B_z(rhoP, zP);

		// add contribution from B_z of wire loop to result
		magneticField[3 * idxEval + 0] += bZ * eX;
		magneticField[3 * idxEval + 1] += bZ * eY;
		magneticField[3 * idxEval + 2] += bZ * eZ;
	}
}

// --------------------------------------------------

/**
 * Compute the magnetic vector potential of a polygon filament
 * at a number of evaluation locations.
 *
 * @param vertices [3: x, y, z][numVertices] points along polygon; in m
 * @param current current along polygon; in A
 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
 * @param vectorPotential [3: x, y, z][numEvalPos] target array for magnetic vector potential at evaluation locations; in Tm
 * @param idxSourceStart first index in {@code vertices} to take into account
 * @param idxSourceEnd (last+1) index in {@code vertices} to take into account
 * @param idxEvalStart first index in {@code evalPos} to take into account
 * @param idxEvalEnd (last+1) index in {@code evalPos} to take into account
 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
 */
void kernelVectorPotentialPolygonFilament(
		double *vertices,
		double current,
		double *evalPos,
		double *vectorPotential,
		int idxSourceStart,
		int idxSourceEnd,
		int idxEvalStart,
		int idxEvalEnd,
		bool useCompensatedSummation) {

	double aPrefactor = MU_0_BY_2_PI * current;

	// setup compensated summation storage
	double *aXSum;
	double *aYSum;
	double *aZSum;
	if (useCompensatedSummation) {
		int numEvalPos = idxEvalEnd - idxEvalStart;

		// need three doubles (s, cs, ccs) per eval pos --> see compsum.h
		int numBytesToAllocate = 3 * numEvalPos * sizeof(double);

		aXSum = (double *) malloc (numBytesToAllocate);
		aYSum = (double *) malloc (numBytesToAllocate);
		aZSum = (double *) malloc (numBytesToAllocate);

		// initialize storage to zero
		memset(aXSum, 0, numBytesToAllocate);
		memset(aYSum, 0, numBytesToAllocate);
		memset(aZSum, 0, numBytesToAllocate);
	} else {
		aXSum = NULL;
		aYSum = NULL;
		aZSum = NULL;
	}

	double x_i = vertices[idxSourceStart * 3 + 0];
	double y_i = vertices[idxSourceStart * 3 + 1];
	double z_i = vertices[idxSourceStart * 3 + 2];

	for (int idxSource = idxSourceStart; idxSource < idxSourceEnd; ++idxSource) {

		double x_f = vertices[(idxSource + 1) * 3 + 0];
		double y_f = vertices[(idxSource + 1) * 3 + 1];
		double z_f = vertices[(idxSource + 1) * 3 + 2];

		// vector from start to end of i:th wire segment
		double dx = x_f - x_i;
		double dy = y_f - y_i;
		double dz = z_f - z_i;

		// squared length of wire segment
		double l2 = dx * dx + dy * dy + dz * dz;
		if (l2 == 0.0) {
			// skip zero-length segments: no contribution
			continue;
		}

		// length of wire segment
		double l = sqrt(l2);

		// unit vector parallel to wire segment
		double eX = dx / l;
		double eY = dy / l;
		double eZ = dz / l;

		for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {

			// vector from start of wire segment to eval pos
			double r0x = evalPos[3 * idxEval + 0] - x_i;
			double r0y = evalPos[3 * idxEval + 1] - y_i;
			double r0z = evalPos[3 * idxEval + 2] - z_i;

			// z position along axis of wire segment
			double alignedZ = eX * r0x + eY * r0y + eZ * r0z;

			// normalized z component of evaluation location in coordinate system of wire segment
			double zP = alignedZ / l;

			// vector perpendicular to axis of wire segment, pointing at evaluation pos
			double rPerpX = r0x - alignedZ * eX;
			double rPerpY = r0y - alignedZ * eY;
			double rPerpZ = r0z - alignedZ * eZ;

			// perpendicular distance between evalPos and axis of wire segment
			double alignedR = sqrt(rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ);

			// normalized rho component of evaluation location in coordinate system of wire segment
			double rhoP = alignedR / l;

			// compute parallel component of magnetic vector potential, including current and mu_0
			double aParallel = aPrefactor * straightWireSegment_A_z(rhoP, zP);

			// add contribution from wire segment to result
			if (useCompensatedSummation) {
				compAdd(aParallel * eX, &(aXSum[3 * (idxEval - idxEvalStart)]));
				compAdd(aParallel * eY, &(aYSum[3 * (idxEval - idxEvalStart)]));
				compAdd(aParallel * eZ, &(aZSum[3 * (idxEval - idxEvalStart)]));
			} else {
				vectorPotential[3 * idxEval + 0] += aParallel * eX;
				vectorPotential[3 * idxEval + 1] += aParallel * eY;
				vectorPotential[3 * idxEval + 2] += aParallel * eZ;
			}
		}

		// shift to next point
		x_i = x_f;
		y_i = y_f;
		z_i = z_f;
	}

	if (useCompensatedSummation) {
		// obtain compensated sums from summation objects
		for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {
			int relIdx = 3 * (idxEval - idxEvalStart);
			vectorPotential[3 * idxEval + 0] += aXSum[relIdx + 0] + aXSum[relIdx + 1] + aXSum[relIdx + 2];
			vectorPotential[3 * idxEval + 1] += aYSum[relIdx + 0] + aYSum[relIdx + 1] + aYSum[relIdx + 2];
			vectorPotential[3 * idxEval + 2] += aZSum[relIdx + 0] + aZSum[relIdx + 1] + aZSum[relIdx + 2];
		}

		free(aXSum);
		free(aYSum);
		free(aZSum);
	}
}

/**
 * Compute the magnetic vector potential of a polygon filament
 * at a number of evaluation locations.
 *
 * @param void (*vertexSupplier)(int i, double *point): callback to put i-th current carrier polygon vertex into point as [3: x, y, z]; in m
 * @param current current along polygon; in A
 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
 * @param vectorPotential [3: x, y, z][numEvalPos] target array for magnetic vector potential at evaluation locations; in Tm
 * @param idxSourceStart first index in {@code vertices} to take into account
 * @param idxSourceEnd (last+1) index in {@code vertices} to take into account
 * @param idxEvalStart first index in {@code evalPos} to take into account
 * @param idxEvalEnd (last+1) index in {@code evalPos} to take into account
 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
 */
void kernelVectorPotentialPolygonFilament(
		void (*vertexSupplier)(int i, double *point),
		double current,
		double *evalPos,
		double *vectorPotential,
		int idxSourceStart,
		int idxSourceEnd,
		int idxEvalStart,
		int idxEvalEnd,
		bool useCompensatedSummation) {

	double aPrefactor = MU_0_BY_2_PI * current;

	// setup compensated summation storage
	double *aXSum;
	double *aYSum;
	double *aZSum;
	if (useCompensatedSummation) {
		int numEvalPos = idxEvalEnd - idxEvalStart;

		// need three doubles (s, cs, ccs) per eval pos --> see compsum.h
		int numBytesToAllocate = 3 * numEvalPos * sizeof(double);

		aXSum = (double *) malloc (numBytesToAllocate);
		aYSum = (double *) malloc (numBytesToAllocate);
		aZSum = (double *) malloc (numBytesToAllocate);

		// initialize storage to zero
		memset(aXSum, 0, numBytesToAllocate);
		memset(aYSum, 0, numBytesToAllocate);
		memset(aZSum, 0, numBytesToAllocate);
	} else {
		aXSum = NULL;
		aYSum = NULL;
		aZSum = NULL;
	}

	// get first point from vertexSupplier
	double pointData[3];
	vertexSupplier(idxSourceStart, pointData);
	double x_i = pointData[0];
	double y_i = pointData[1];
	double z_i = pointData[2];

	for (int idxSource = idxSourceStart; idxSource < idxSourceEnd; ++idxSource) {

		// get next point from vertexSupplier
		vertexSupplier(idxSource+1, pointData);
		double x_f = pointData[0];
		double y_f = pointData[1];
		double z_f = pointData[2];

		// vector from start to end of i:th wire segment
		double dx = x_f - x_i;
		double dy = y_f - y_i;
		double dz = z_f - z_i;

		// squared length of wire segment
		double l2 = dx * dx + dy * dy + dz * dz;
		if (l2 == 0.0) {
			// skip zero-length segments: no contribution
			continue;
		}

		// length of wire segment
		double l = sqrt(l2);

		// unit vector parallel to wire segment
		double eX = dx / l;
		double eY = dy / l;
		double eZ = dz / l;

		for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {

			// vector from start of wire segment to eval pos
			double r0x = evalPos[3 * idxEval + 0] - x_i;
			double r0y = evalPos[3 * idxEval + 1] - y_i;
			double r0z = evalPos[3 * idxEval + 2] - z_i;

			// z position along axis of wire segment
			double alignedZ = eX * r0x + eY * r0y + eZ * r0z;

			// normalized z component of evaluation location in coordinate system of wire segment
			double zP = alignedZ / l;

			// vector perpendicular to axis of wire segment, pointing at evaluation pos
			double rPerpX = r0x - alignedZ * eX;
			double rPerpY = r0y - alignedZ * eY;
			double rPerpZ = r0z - alignedZ * eZ;

			// perpendicular distance between evalPos and axis of wire segment
			double alignedR = sqrt(rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ);

			// normalized rho component of evaluation location in coordinate system of wire segment
			double rhoP = alignedR / l;

			// compute parallel component of magnetic vector potential, including current and mu_0
			double aParallel = aPrefactor * straightWireSegment_A_z(rhoP, zP);

			// add contribution from wire segment to result
			if (useCompensatedSummation) {
				compAdd(aParallel * eX, &(aXSum[3 * (idxEval - idxEvalStart)]));
				compAdd(aParallel * eY, &(aYSum[3 * (idxEval - idxEvalStart)]));
				compAdd(aParallel * eZ, &(aZSum[3 * (idxEval - idxEvalStart)]));
			} else {
				vectorPotential[3 * idxEval + 0] += aParallel * eX;
				vectorPotential[3 * idxEval + 1] += aParallel * eY;
				vectorPotential[3 * idxEval + 2] += aParallel * eZ;
			}
		}

		// shift to next point
		x_i = x_f;
		y_i = y_f;
		z_i = z_f;
	}

	if (useCompensatedSummation) {
		// obtain compensated sums from summation objects
		for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {
			int relIdx = 3 * (idxEval - idxEvalStart);
			vectorPotential[3 * idxEval + 0] += aXSum[relIdx + 0] + aXSum[relIdx + 1] + aXSum[relIdx + 2];
			vectorPotential[3 * idxEval + 1] += aYSum[relIdx + 0] + aYSum[relIdx + 1] + aYSum[relIdx + 2];
			vectorPotential[3 * idxEval + 2] += aZSum[relIdx + 0] + aZSum[relIdx + 1] + aZSum[relIdx + 2];
		}

		free(aXSum);
		free(aYSum);
		free(aZSum);
	}
}

/**
 * Compute the magnetic field of a polygon filament
 * at a number of evaluation locations.
 *
 * @param vertices [3: x, y, z][numVertices] points along polygon; in m
 * @param current current along polygon; in A
 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
 * @param magneticField [3: x, y, z][numEvalPos] target array for magnetic field at evaluation locations; in T
 * @param idxSourceStart first index in {@code vertices} to take into account
 * @param idxSourceEnd (last+1) index in {@code vertices} to take into account
 * @param idxEvalStart first index in {@code evalPos} to take into account
 * @param idxEvalEnd (last+1) index in {@code evalPos} to take into account
 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
 */
void kernelMagneticFieldPolygonFilament(
		double *vertices,
		double current,
		double *evalPos,
		double *magneticField,
		int idxSourceStart,
		int idxSourceEnd,
		int idxEvalStart,
		int idxEvalEnd,
		bool useCompensatedSummation) {

	// needs additional division by length of wire segment!
	double bPrefactorL = MU_0_BY_4_PI * current;

	// setup compensated summation storage
	double *bXSum;
	double *bYSum;
	double *bZSum;
	if (useCompensatedSummation) {
		int numEvalPos = idxEvalEnd - idxEvalStart;

		// need three doubles (s, cs, ccs) per eval pos --> see compsum.h
		int numBytesToAllocate = 3 * numEvalPos * sizeof(double);

		bXSum = (double *) malloc (numBytesToAllocate);
		bYSum = (double *) malloc (numBytesToAllocate);
		bZSum = (double *) malloc (numBytesToAllocate);

		// initialize storage to zero
		memset(bXSum, 0, numBytesToAllocate);
		memset(bYSum, 0, numBytesToAllocate);
		memset(bZSum, 0, numBytesToAllocate);
	} else {
		bXSum = NULL;
		bYSum = NULL;
		bZSum = NULL;
	}

	double x_i = vertices[3 * idxSourceStart + 0];
	double y_i = vertices[3 * idxSourceStart + 1];
	double z_i = vertices[3 * idxSourceStart + 2];

	for (int idxSource = idxSourceStart; idxSource < idxSourceEnd; ++idxSource) {

		double x_f = vertices[3 * (idxSource + 1) + 0];
		double y_f = vertices[3 * (idxSource + 1) + 1];
		double z_f = vertices[3 * (idxSource + 1) + 2];

		// vector from start to end of i:th wire segment
		double dx = x_f - x_i;
		double dy = y_f - y_i;
		double dz = z_f - z_i;

		// squared length of wire segment
		double l2 = dx * dx + dy * dy + dz * dz;
		if (l2 == 0.0) {
			// skip zero-length segments: no contribution
			continue;
		}

		// length of wire segment
		double l = sqrt(l2);

		// assemble full prefactor for B_phi
		double bPrefactor = bPrefactorL / l;

		// unit vector parallel to wire segment
		double eX = dx / l;
		double eY = dy / l;
		double eZ = dz / l;

		for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {

			// vector from start of wire segment to eval pos
			double r0x = evalPos[3 * idxEval + 0] - x_i;
			double r0y = evalPos[3 * idxEval + 1] - y_i;
			double r0z = evalPos[3 * idxEval + 2] - z_i;

			// z position along axis of wire segment
			double alignedZ = eX * r0x + eY * r0y + eZ * r0z;

			// vector perpendicular to axis of wire segment, pointing at evaluation pos
			double rPerpX = r0x - alignedZ * eX;
			double rPerpY = r0y - alignedZ * eY;
			double rPerpZ = r0z - alignedZ * eZ;

			// perpendicular distance squared between evalPos and axis of wire segment
			double alignedRSq = rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ;

			// B_phi is zero along axis of filament
			if (alignedRSq > 0.0) {

				// perpendicular distance between evalPos and axis of wire segment
				double alignedR = sqrt(alignedRSq);

				// normalized rho component of evaluation location in coordinate system of wire segment
				double rhoP = alignedR / l;

				// normalized z component of evaluation location in coordinate system of wire segment
				double zP = alignedZ / l;

				// compute tangential component of magnetic vector potential, including current and mu_0
				double bPhi = bPrefactor * straightWireSegment_B_phi(rhoP, zP);

				// unit vector in radial direction
				double eRX = rPerpX / alignedR;
				double eRY = rPerpY / alignedR;
				double eRZ = rPerpZ / alignedR;

				// compute cross product between e_z and e_rho to get e_phi
				double ePhiX = eY * eRZ - eZ * eRY;
				double ePhiY = eZ * eRX - eX * eRZ;
				double ePhiZ = eX * eRY - eY * eRX;

				// add contribution from wire segment to result
				if (useCompensatedSummation) {
					compAdd(bPhi * ePhiX, &(bXSum[3 * (idxEval - idxEvalStart)]));
					compAdd(bPhi * ePhiY, &(bYSum[3 * (idxEval - idxEvalStart)]));
					compAdd(bPhi * ePhiZ, &(bZSum[3 * (idxEval - idxEvalStart)]));
				} else {
					magneticField[3 * idxEval + 0] += bPhi * ePhiX;
					magneticField[3 * idxEval + 1] += bPhi * ePhiY;
					magneticField[3 * idxEval + 2] += bPhi * ePhiZ;
				}
			}
		}

		// shift to next point
		x_i = x_f;
		y_i = y_f;
		z_i = z_f;
	}

	if (useCompensatedSummation) {
		// obtain compensated sums from summation objects
		for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {
			int relIdx = 3 * (idxEval - idxEvalStart);
			magneticField[3 * idxEval + 0] += bXSum[relIdx + 0] + bXSum[relIdx + 1] + bXSum[relIdx + 2];
			magneticField[3 * idxEval + 1] += bYSum[relIdx + 0] + bYSum[relIdx + 1] + bYSum[relIdx + 2];
			magneticField[3 * idxEval + 2] += bZSum[relIdx + 0] + bZSum[relIdx + 1] + bZSum[relIdx + 2];
		}

		free(bXSum);
		free(bYSum);
		free(bZSum);
	}
}

/**
 * Compute the magnetic field of a polygon filament
 * at a number of evaluation locations.
 *
 * @param void (*vertexSupplier)(int i, double *point): callback to put i-th current carrier polygon vertex into point as [3: x, y, z]; in m
 * @param current current along polygon; in A
 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
 * @param magneticField [3: x, y, z][numEvalPos] target array for magnetic field at evaluation locations; in T
 * @param idxSourceStart first index in {@code vertices} to take into account
 * @param idxSourceEnd (last+1) index in {@code vertices} to take into account
 * @param idxEvalStart first index in {@code evalPos} to take into account
 * @param idxEvalEnd (last+1) index in {@code evalPos} to take into account
 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
 */
void kernelMagneticFieldPolygonFilament(
		void (*vertexSupplier)(int i, double *point),
		double current,
		double *evalPos,
		double *magneticField,
		int idxSourceStart,
		int idxSourceEnd,
		int idxEvalStart,
		int idxEvalEnd,
		bool useCompensatedSummation) {

	// needs additional division by length of wire segment!
	double bPrefactorL = MU_0_BY_4_PI * current;

	// setup compensated summation storage
	double *bXSum;
	double *bYSum;
	double *bZSum;
	if (useCompensatedSummation) {
		int numEvalPos = idxEvalEnd - idxEvalStart;

		// need three doubles (s, cs, ccs) per eval pos --> see compsum.h
		int numBytesToAllocate = 3 * numEvalPos * sizeof(double);

		bXSum = (double *) malloc (numBytesToAllocate);
		bYSum = (double *) malloc (numBytesToAllocate);
		bZSum = (double *) malloc (numBytesToAllocate);

		// initialize storage to zero
		memset(bXSum, 0, numBytesToAllocate);
		memset(bYSum, 0, numBytesToAllocate);
		memset(bZSum, 0, numBytesToAllocate);
	} else {
		bXSum = NULL;
		bYSum = NULL;
		bZSum = NULL;
	}

	// get first point from vertexSupplier
	double pointData[3];
	vertexSupplier(idxSourceStart, pointData);
	double x_i = pointData[0];
	double y_i = pointData[1];
	double z_i = pointData[2];

	for (int idxSource = idxSourceStart; idxSource < idxSourceEnd; ++idxSource) {

		// get next point from vertexSupplier
		vertexSupplier(idxSource+1, pointData);
		double x_f = pointData[0];
		double y_f = pointData[1];
		double z_f = pointData[2];

		// vector from start to end of i:th wire segment
		double dx = x_f - x_i;
		double dy = y_f - y_i;
		double dz = z_f - z_i;

		// squared length of wire segment
		double l2 = dx * dx + dy * dy + dz * dz;
		if (l2 == 0.0) {
			// skip zero-length segments: no contribution
			continue;
		}

		// length of wire segment
		double l = sqrt(l2);

		// assemble full prefactor for B_phi
		double bPrefactor = bPrefactorL / l;

		// unit vector parallel to wire segment
		double eX = dx / l;
		double eY = dy / l;
		double eZ = dz / l;

		for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {

			// vector from start of wire segment to eval pos
			double r0x = evalPos[3 * idxEval + 0] - x_i;
			double r0y = evalPos[3 * idxEval + 1] - y_i;
			double r0z = evalPos[3 * idxEval + 2] - z_i;

			// z position along axis of wire segment
			double alignedZ = eX * r0x + eY * r0y + eZ * r0z;

			// vector perpendicular to axis of wire segment, pointing at evaluation pos
			double rPerpX = r0x - alignedZ * eX;
			double rPerpY = r0y - alignedZ * eY;
			double rPerpZ = r0z - alignedZ * eZ;

			// perpendicular distance squared between evalPos and axis of wire segment
			double alignedRSq = rPerpX * rPerpX + rPerpY * rPerpY + rPerpZ * rPerpZ;

			// B_phi is zero along axis of filament
			if (alignedRSq > 0.0) {

				// perpendicular distance between evalPos and axis of wire segment
				double alignedR = sqrt(alignedRSq);

				// normalized rho component of evaluation location in coordinate system of wire segment
				double rhoP = alignedR / l;

				// normalized z component of evaluation location in coordinate system of wire segment
				double zP = alignedZ / l;

				// compute tangential component of magnetic vector potential, including current and mu_0
				double bPhi = bPrefactor * straightWireSegment_B_phi(rhoP, zP);

				// unit vector in radial direction
				double eRX = rPerpX / alignedR;
				double eRY = rPerpY / alignedR;
				double eRZ = rPerpZ / alignedR;

				// compute cross product between e_z and e_rho to get e_phi
				double ePhiX = eY * eRZ - eZ * eRY;
				double ePhiY = eZ * eRX - eX * eRZ;
				double ePhiZ = eX * eRY - eY * eRX;

				// add contribution from wire segment to result
				if (useCompensatedSummation) {
					compAdd(bPhi * ePhiX, &(bXSum[3 * (idxEval - idxEvalStart)]));
					compAdd(bPhi * ePhiY, &(bYSum[3 * (idxEval - idxEvalStart)]));
					compAdd(bPhi * ePhiZ, &(bZSum[3 * (idxEval - idxEvalStart)]));
				} else {
					magneticField[3 * idxEval + 0] += bPhi * ePhiX;
					magneticField[3 * idxEval + 1] += bPhi * ePhiY;
					magneticField[3 * idxEval + 2] += bPhi * ePhiZ;
				}
			}
		}

		// shift to next point
		x_i = x_f;
		y_i = y_f;
		z_i = z_f;
	}

	if (useCompensatedSummation) {
		// obtain compensated sums from summation objects
		for (int idxEval = idxEvalStart; idxEval < idxEvalEnd; ++idxEval) {
			int relIdx = 3 * (idxEval - idxEvalStart);
			magneticField[3 * idxEval + 0] += bXSum[relIdx + 0] + bXSum[relIdx + 1] + bXSum[relIdx + 2];
			magneticField[3 * idxEval + 1] += bYSum[relIdx + 0] + bYSum[relIdx + 1] + bYSum[relIdx + 2];
			magneticField[3 * idxEval + 2] += bZSum[relIdx + 0] + bZSum[relIdx + 1] + bZSum[relIdx + 2];
		}

		free(bXSum);
		free(bYSum);
		free(bZSum);
	}
}

// --------------------------------------------------

static void polygonFilamentSoA(
		bool potential, int numVertices,
		const double *vx, const double *vy, const double *vz, double current,
		int numEvalPos, const double *ex, const double *ey, const double *ez,
		double *ox, double *oy, double *oz,
		int numProcessors, bool useCompensatedSummation);

// array-of-structs interface: convert to SoA and add the result to the target
static void polygonFilamentAoS(
		bool potential, int numVertices, const double *vertices, double current,
		int numEvalPos, const double *evalPos, double *target,
		int numProcessors, bool useCompensatedSummation) {
	if (numVertices < 2 || numProcessors < 1 || current == 0.0) {
		// let the SoA routine report the problem
		polygonFilamentSoA(potential, numVertices, NULL, NULL, NULL, current, 0, NULL, NULL, NULL,
				NULL, NULL, NULL, numProcessors, useCompensatedSummation);
		return;
	}
	const size_t nv = (size_t) numVertices, ne = (size_t) numEvalPos;
	double *buffer = (double *) calloc(3 * nv + 6 * ne, sizeof(double));
	if (buffer == NULL) {
		printf("failed to allocate temporary arrays\n");
		return;
	}
	double *v[3] = {buffer, buffer + nv, buffer + 2 * nv};
	double *e[3] = {buffer + 3 * nv, buffer + 3 * nv + ne, buffer + 3 * nv + 2 * ne};
	double *o[3] = {buffer + 3 * nv + 3 * ne, buffer + 3 * nv + 4 * ne, buffer + 3 * nv + 5 * ne};
	for (size_t i = 0; i < nv; ++i) {
		for (int c = 0; c < 3; ++c) {
			v[c][i] = vertices[3 * i + c];
		}
	}
	for (size_t i = 0; i < ne; ++i) {
		for (int c = 0; c < 3; ++c) {
			e[c][i] = evalPos[3 * i + c];
		}
	}
	polygonFilamentSoA(potential, numVertices, v[0], v[1], v[2], current, numEvalPos, e[0], e[1], e[2],
			o[0], o[1], o[2], numProcessors, useCompensatedSummation);
	for (size_t i = 0; i < ne; ++i) {
		for (int c = 0; c < 3; ++c) {
			target[3 * i + c] += o[c][i];
		}
	}
	free(buffer);
}

/**
 * Compute the magnetic vector potential of a polygon filament
 * at a number of evaluation locations.
 *
 * @param numVertices number of points along polygon
 * @param vertices [numVertices][3: x, y, z] points along polygon; in m
 * @param current current along polygon; in A
 * @param numEvalPos number of evaluation locations
 * @param evalPos [numEvalPos][3: x, y, z] evaluation locations; in m
 * @param vectorPotential [numEvalPos][3: x, y, z] target array for magnetic vector potential at evaluation locations; in Tm
 * @param numProcessors number of processors to use for parallelization
 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
 */
void vectorPotentialPolygonFilament(
		int numVertices,
		double *vertices,
		double current,
		int numEvalPos,
		double *evalPos,
		double *vectorPotential,
		int numProcessors,
		bool useCompensatedSummation) {
	polygonFilamentAoS(true, numVertices, vertices, current, numEvalPos, evalPos, vectorPotential,
			numProcessors, useCompensatedSummation);
}

/**
 * Compute the magnetic vector potential of a polygon filament
 * at a number of evaluation locations.
 *
 * @param numVertices number of vertices to obtain from vertexSupplier
 * @param vertexSupplier callback to put i-th current carrier polygon vertex into point as [3: x, y, z]; in m
 * @param current current along polygon; in A
 * @param numEvalPos number of evaluation locations
 * @param evalPos [numEvalPos][3: x, y, z] evaluation locations; in m
 * @param vectorPotential [numEvalPos][3: x, y, z] target array for magnetic vector potential at evaluation locations; in Tm
 * @param numProcessors number of processors to use for parallelization
 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
 */
void vectorPotentialPolygonFilament(
		int numVertices,
		void (*vertexSupplier)(int i, double *point),
		double current,
		int numEvalPos,
		double *evalPos,
		double *vectorPotential,
		int numProcessors,
		bool useCompensatedSummation) {

	if (numVertices < 2) {
		printf("need at least 2 vertices, but only got %d\n", numVertices);
		return;
	}

	int numSegments = numVertices - 1;

	if (numProcessors < 1) {
		printf("need at least 1 processor, but only got %d\n", numProcessors);
		return;
	}

	if (current == 0.0) {
		return;
	}

	if (numProcessors == 1) {
		// single-threaded call
		int idxSourceStart = 0;
		int idxSourceEnd   = numSegments;
		int idxEvalStart   = 0;
		int idxEvalEnd     = numEvalPos;
		kernelVectorPotentialPolygonFilament(
				vertexSupplier, current,
				evalPos,
				vectorPotential,
				idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
				useCompensatedSummation);
	} else {
		// use multithreading

		if (numSegments > numEvalPos) {
			// parallelize over numSegments

			// Note that each thread needs its own copy of the vectorPotential array,
			// so this approach might need quite some memory in case the number of
			// threads and the number of evaluation points is large.

			int nThreads;
			int nSourcePerThread;
			int nSourceRemainder;
			if (numSegments < numProcessors) {
				nThreads = numSegments;

				nSourcePerThread = 1;
				nSourceRemainder = 0;
			} else {
				nThreads = numProcessors;

				nSourcePerThread = numSegments / nThreads;
				nSourceRemainder = numSegments % nThreads;
			}
			
			int idxEvalStart   = 0;
			int idxEvalEnd     = numEvalPos;

			int numBytes = nThreads * numEvalPos * 3 * sizeof(double);
			double *vectorPotentialContributions = (double *) malloc(numBytes);
			if (vectorPotentialContributions == NULL) {
				printf("failed to allocate temporary array for vector potential contributions\n");
				return;
			}
			memset(vectorPotentialContributions, 0, numBytes);

#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
			for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
				
				int idxSourceStart =  idxThread      * nSourcePerThread;
				int idxSourceEnd   = (idxThread + 1) * nSourcePerThread;
				if (idxThread < nSourceRemainder) {
					idxSourceStart += idxThread;
					idxSourceEnd   += idxThread + 1;
				} else {
					idxSourceStart += nSourceRemainder;
					idxSourceEnd   += nSourceRemainder;
				}

				kernelVectorPotentialPolygonFilament(
						vertexSupplier, current,
						evalPos,
						&(vectorPotentialContributions[idxThread * numEvalPos * 3]),
						idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
						useCompensatedSummation);
			}

			// sum up contributions from source chunks
			if (useCompensatedSummation) {
				for (int i=0; i<numEvalPos; ++i) {
					double sumX[3] = { 0.0, 0.0, 0.0 };
					double sumY[3] = { 0.0, 0.0, 0.0 };
					double sumZ[3] = { 0.0, 0.0, 0.0 };
					for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
						// TODO: bad memory access pattern here --> potential bottleneck !!!
						compAdd(vectorPotentialContributions[(idxThread * numEvalPos + i) * 3 + 0], sumX);
						compAdd(vectorPotentialContributions[(idxThread * numEvalPos + i) * 3 + 1], sumY);
						compAdd(vectorPotentialContributions[(idxThread * numEvalPos + i) * 3 + 2], sumZ);
					}
					vectorPotential[i * 3 + 0] += sumX[0] + sumX[1] + sumX[2];
					vectorPotential[i * 3 + 1] += sumY[0] + sumY[1] + sumY[2];
					vectorPotential[i * 3 + 2] += sumZ[0] + sumZ[1] + sumZ[2];
				}
			} else {
				for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
					for (int i=0; i<numEvalPos; ++i) {
						vectorPotential[i * 3 + 0] += vectorPotentialContributions[(idxThread * numEvalPos + i) * 3 + 0];
						vectorPotential[i * 3 + 1] += vectorPotentialContributions[(idxThread * numEvalPos + i) * 3 + 1];
						vectorPotential[i * 3 + 2] += vectorPotentialContributions[(idxThread * numEvalPos + i) * 3 + 2];
					}
				}
			}

			free(vectorPotentialContributions);
		} else { // numEvalPos > numSegments
			// parallelize over numEvalPos

			int nThreads;
			int idxSourceStart = 0;
			int idxSourceEnd   = numSegments;
			int nEvalPerThread;
			int nEvalRemainder;
			if (numEvalPos < numProcessors) {
				nThreads = numEvalPos;

				nEvalPerThread = 1;
				nEvalRemainder = 0;
			} else {
				nThreads = numProcessors;

				nEvalPerThread = numEvalPos / nThreads;
				nEvalRemainder = numEvalPos % nThreads;
			}

#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
			for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
				
				int idxEvalStart =  idxThread      * nEvalPerThread;
				int idxEvalEnd   = (idxThread + 1) * nEvalPerThread;
				if (idxThread < nEvalRemainder) {
					idxEvalStart += idxThread;
					idxEvalEnd   += idxThread + 1;
				} else {
					idxEvalStart += nEvalRemainder;
					idxEvalEnd   += nEvalRemainder;
				}

				kernelVectorPotentialPolygonFilament(
						vertexSupplier, current,
						evalPos,
						vectorPotential,
						idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
						useCompensatedSummation);
			}
		} // parallelize over nSource or nEval
	} // parallelization
}

/**
 * Compute the magnetic field of a polygon filament
 * at a number of evaluation locations.
 *
 * @param numVertices number of points along polygon
 * @param vertices [numVertices][3: x, y, z] points along polygon; in m
 * @param current current along polygon; in A
 * @param numEvalPos number of evaluation locations
 * @param evalPos [numEvalPos][3: x, y, z] evaluation locations; in m
 * @param magneticField [numEvalPos][3: x, y, z] target array for magnetic field at evaluation locations; inT
 * @param numProcessors number of processors to use for parallelization
 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
 */
void magneticFieldPolygonFilament(
		int numVertices,
		double *vertices,
		double current,
		int numEvalPos,
		double *evalPos,
		double *magneticField,
		int numProcessors,
		bool useCompensatedSummation) {
	polygonFilamentAoS(false, numVertices, vertices, current, numEvalPos, evalPos, magneticField,
			numProcessors, useCompensatedSummation);
}


/**
 * Compute the magnetic field of a polygon filament
 * at a number of evaluation locations.
 *
 * @param numVertices number of vertices to obtain from vertexSupplier
 * @param void (*vertexSupplier)(int i, double *point): callback to put i-th current carrier polygon vertex into point as [3: x, y, z]; in m
 * @param current current along polygon; in A
 * @param numEvalPos number of evaluation locations
 * @param evalPos [numEvalPos][3: x, y, z] evaluation locations; in m
 * @param magneticField [numEvalPos][3: x, y, z] target array for magnetic field at evaluation locations; in T
 * @param numProcessors number of processors to use for parallelization
 * @param useCompensatedSummation if true, use Kahan-Babuska compensated summation to compute the superposition
 *                                of the contributions from the polygon vertices; otherwise, use standard += summation
 */
void magneticFieldPolygonFilament(
		int numVertices,
		void (*vertexSupplier)(int i, double *point),
		double current,
		int numEvalPos,
		double *evalPos,
		double *magneticField,
		int numProcessors,
		bool useCompensatedSummation) {

	if (numVertices < 2) {
		printf("need at least 2 vertices, but only got %d\n", numVertices);
		return;
	}

	// number of straight wire segments
	int numSegments = numVertices - 1;

	if (numProcessors < 1) {
		printf("need at least 1 processor, but only got %d\n", numProcessors);
		return;
	}

	if (current == 0.0) {
		return;
	}

	if (numProcessors == 1) {
		// single-threaded call
		int idxSourceStart = 0;
		int idxSourceEnd   = numSegments;
		int idxEvalStart   = 0;
		int idxEvalEnd     = numEvalPos;
		kernelMagneticFieldPolygonFilament(
				vertexSupplier, current,
				evalPos,
				magneticField,
				idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
				useCompensatedSummation);
	} else {
		// use multithreading

		if (numSegments > numEvalPos) {
			// parallelize over numSegments

			// Note that each thread needs its own copy of the magneticField array,
			// so this approach might need quite some memory in case the number of
			// threads and the number of evaluation points is large.

			int nThreads;
			int nSourcePerThread;
			int nSourceRemainder;
			if (numSegments < numProcessors) {
				nThreads = numSegments;

				nSourcePerThread = 1;
				nSourceRemainder = 0;
			} else {
				nThreads = numProcessors;

				nSourcePerThread = numSegments / nThreads;
				nSourceRemainder = numSegments % nThreads;
			}
			
			int idxEvalStart   = 0;
			int idxEvalEnd     = numEvalPos;

			int numBytes = nThreads * numEvalPos * 3 * sizeof(double);
			double *magneticFieldContributions = (double *) malloc(numBytes);
			if (magneticFieldContributions == NULL) {
				printf("failed to allocate temporary array for magnetic field contributions\n");
				return;
			}
			memset(magneticFieldContributions, 0, numBytes);

#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
			for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
				
				int idxSourceStart =  idxThread      * nSourcePerThread;
				int idxSourceEnd   = (idxThread + 1) * nSourcePerThread;
				if (idxThread < nSourceRemainder) {
					idxSourceStart += idxThread;
					idxSourceEnd   += idxThread + 1;
				} else {
					idxSourceStart += nSourceRemainder;
					idxSourceEnd   += nSourceRemainder;
				}

				kernelMagneticFieldPolygonFilament(
						vertexSupplier, current,
						evalPos,
						&(magneticFieldContributions[idxThread * numEvalPos * 3]),
						idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
						useCompensatedSummation);
			}

			// sum up contributions from source chunks
			if (useCompensatedSummation) {
				for (int i = 0; i < numEvalPos; ++i) {
					double sumX[3] = { 0.0, 0.0, 0.0 };
					double sumY[3] = { 0.0, 0.0, 0.0 };
					double sumZ[3] = { 0.0, 0.0, 0.0 };
					for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
						// TODO: bad memory access pattern here --> potential bottleneck !!!
						compAdd(magneticFieldContributions[(idxThread * numEvalPos + i) * 3 + 0], sumX);
						compAdd(magneticFieldContributions[(idxThread * numEvalPos + i) * 3 + 1], sumY);
						compAdd(magneticFieldContributions[(idxThread * numEvalPos + i) * 3 + 2], sumZ);
					}
					magneticField[i * 3 + 0] += sumX[0] + sumX[1] + sumX[2];
					magneticField[i * 3 + 1] += sumY[0] + sumY[1] + sumY[2];
					magneticField[i * 3 + 2] += sumZ[0] + sumZ[1] + sumZ[2];
				}
			} else {
				for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
					for (int i = 0; i < numEvalPos; ++i) {
						magneticField[i * 3 + 0] += magneticFieldContributions[(idxThread * numEvalPos + i) * 3 + 0];
						magneticField[i * 3 + 1] += magneticFieldContributions[(idxThread * numEvalPos + i) * 3 + 1];
						magneticField[i * 3 + 2] += magneticFieldContributions[(idxThread * numEvalPos + i) * 3 + 2];
					}
				}
			}

			free(magneticFieldContributions);
		} else { // numEvalPos > numSegments
			// parallelize over numEvalPos

			int nThreads;
			int idxSourceStart = 0;
			int idxSourceEnd   = numSegments;
			int nEvalPerThread;
			int nEvalRemainder;
			if (numEvalPos < numProcessors) {
				nThreads = numEvalPos;

				nEvalPerThread = 1;
				nEvalRemainder = 0;
			} else {
				nThreads = numProcessors;

				nEvalPerThread = numEvalPos / nThreads;
				nEvalRemainder = numEvalPos % nThreads;
			}

#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
			for (int idxThread = 0; idxThread < nThreads; ++idxThread) {
				
				int idxEvalStart =  idxThread      * nEvalPerThread;
				int idxEvalEnd   = (idxThread + 1) * nEvalPerThread;
				if (idxThread < nEvalRemainder) {
					idxEvalStart += idxThread;
					idxEvalEnd   += idxThread + 1;
				} else {
					idxEvalStart += nEvalRemainder;
					idxEvalEnd   += nEvalRemainder;
				}

				kernelMagneticFieldPolygonFilament(
						vertexSupplier, current,
						evalPos,
						magneticField,
						idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd,
						useCompensatedSummation);
			}
		} // parallelize over numSegments or numEvalPos
	} // parallelization
}

/////// struct-of-arrays (SoA) polygon filament routines

#define SOA_BLOCK 256

// range of normalized coordinates in which sqrt(a*a + b*b) can replace hypot(a, b)
#define SOA_MIN 1.0e-100
#define SOA_MAX 1.0e100

// GCC does not vectorize selects of trapping FP operations (-ftrapping-math).
#if defined(__GNUC__) && !defined(__clang__)
#define SOA_KERNEL __attribute__((optimize("no-trapping-math")))
#else
#define SOA_KERNEL
#endif

SOA_KERNEL static inline void compAddSoA(double c, double &s, double &cs, double &ccs) {
	double t = s + c;
	double c1 = (s - t) + c;
	double c2 = (c - t) + s;
	double d = (fabs(s) >= fabs(c)) ? c1 : c2;
	s = t;
	double t2 = cs + d;
	double d1 = (cs - t2) + d;
	double d2 = (d - t2) + cs;
	cs = t2;
	ccs += (fabs(cs) >= fabs(d)) ? d1 : d2;
}

// Fast path: sqrt-based far-field formulas on blocks of evaluation locations;
// everything else (near field, axis, end planes, extreme magnitudes) uses the scalar routines.
template <bool kMagneticField>
SOA_KERNEL static void kernelPolygonFilamentSoA(
		const double *vx, const double *vy, const double *vz, double current,
		const double *ex, const double *ey, const double *ez,
		double *ox, double *oy, double *oz,
		int idxSourceStart, int idxSourceEnd, int idxEvalStart, int idxEvalEnd,
		bool useCompensatedSummation) {

	const double prefactor = kMagneticField ? MU_0_BY_4_PI * current : MU_0_BY_2_PI * current;

	for (int i0 = idxEvalStart; i0 < idxEvalEnd; i0 += SOA_BLOCK) {
		const int n = (idxEvalEnd - i0 < SOA_BLOCK) ? (idxEvalEnd - i0) : SOA_BLOCK;
		const double *bx = ex + i0, *by = ey + i0, *bz = ez + i0;

		double sx[SOA_BLOCK] = {}, csx[SOA_BLOCK] = {}, ccsx[SOA_BLOCK] = {};
		double sy[SOA_BLOCK] = {}, csy[SOA_BLOCK] = {}, ccsy[SOA_BLOCK] = {};
		double sz[SOA_BLOCK] = {}, csz[SOA_BLOCK] = {}, ccsz[SOA_BLOCK] = {};
		double cx[SOA_BLOCK], cy[SOA_BLOCK], cz[SOA_BLOCK], val[SOA_BLOCK];
		bool slow[SOA_BLOCK];

		for (int k = idxSourceStart; k < idxSourceEnd; ++k) {
			const double x_i = vx[k], y_i = vy[k], z_i = vz[k];
			const double dx = vx[k + 1] - x_i, dy = vy[k + 1] - y_i, dz = vz[k + 1] - z_i;
			const double l2 = dx * dx + dy * dy + dz * dz;
			if (l2 == 0.0) {
				continue;
			}
			const double l = sqrt(l2);
			const double eX = dx / l, eY = dy / l, eZ = dz / l;
			const double scale = kMagneticField ? prefactor / l : prefactor;

#ifdef _OPENMP
#pragma omp simd
#endif // _OPENMP
			for (int j = 0; j < n; ++j) {
				const double r0x = bx[j] - x_i, r0y = by[j] - y_i, r0z = bz[j] - z_i;
				const double alignedZ = eX * r0x + eY * r0y + eZ * r0z;
				const double pX = r0x - alignedZ * eX, pY = r0y - alignedZ * eY, pZ = r0z - alignedZ * eZ;
				const double alignedRSq = pX * pX + pY * pY + pZ * pZ;
				const double alignedR = sqrt(alignedRSq);
				const double rhoP = alignedR / l, zP = alignedZ / l, omz = 1 - zP;

				// bitwise operators keep the loop free of control flow
				const double mi = (rhoP > fabs(zP)) ? rhoP : fabs(zP);
				const double mf = (rhoP > fabs(omz)) ? rhoP : fabs(omz);
				const bool far = kMagneticField ? ((rhoP >= zP) | (rhoP >= omz) | (zP < 0.0) | (zP > 1.0))
						: ((rhoP >= 1.0) | (zP <= -1.0) | (zP > 2.0));
				const bool fast = (kMagneticField ? (alignedRSq > 0.0) : (rhoP != 0.0)) & (zP != 0.0) & (zP != 1.0) & far
						& (mi >= SOA_MIN) & (mi <= SOA_MAX) & (mf >= SOA_MIN) & (mf <= SOA_MAX);

				const double r_i = sqrt(rhoP * rhoP + zP * zP), r_f = sqrt(rhoP * rhoP + omz * omz);
				// B_phi is zero on the axis, so only locations off the axis need the scalar path
				slow[j] = (!fast) & (!kMagneticField | (alignedRSq > 0.0));
				if (kMagneticField) {
					const double bPhi = fast ? scale * (rhoP * (1 / r_i + 1 / r_f)) / (rhoP * rhoP - zP * omz + r_i * r_f) : 0.0;
					const double inv = fast ? 1 / alignedR : 0.0;
					cx[j] = bPhi * (eY * pZ - eZ * pY) * inv;
					cy[j] = bPhi * (eZ * pX - eX * pZ) * inv;
					cz[j] = bPhi * (eX * pY - eY * pX) * inv;
				} else {
					val[j] = fast ? 1 / (r_i + r_f) : 0.5;
				}
			}

			for (int j = 0; j < n; ++j) {
				if (!kMagneticField && !slow[j]) {
					val[j] = scale * atanh(val[j]); // not vectorizable
				}
				if (!slow[j]) {
					continue;
				}
				const double r0x = bx[j] - x_i, r0y = by[j] - y_i, r0z = bz[j] - z_i;
				const double alignedZ = eX * r0x + eY * r0y + eZ * r0z;
				const double pX = r0x - alignedZ * eX, pY = r0y - alignedZ * eY, pZ = r0z - alignedZ * eZ;
				const double alignedR = sqrt(pX * pX + pY * pY + pZ * pZ);
				if (kMagneticField) {
					const double bPhi = scale * straightWireSegment_B_phi(alignedR / l, alignedZ / l) / alignedR;
					cx[j] = bPhi * (eY * pZ - eZ * pY);
					cy[j] = bPhi * (eZ * pX - eX * pZ);
					cz[j] = bPhi * (eX * pY - eY * pX);
				} else {
					val[j] = scale * straightWireSegment_A_z(alignedR / l, alignedZ / l);
				}
			}

			if (!kMagneticField) {
				for (int j = 0; j < n; ++j) {
					cx[j] = val[j] * eX;
					cy[j] = val[j] * eY;
					cz[j] = val[j] * eZ;
				}
			}

			if (useCompensatedSummation) {
#ifdef _OPENMP
#pragma omp simd
#endif // _OPENMP
				for (int j = 0; j < n; ++j) {
					compAddSoA(cx[j], sx[j], csx[j], ccsx[j]);
					compAddSoA(cy[j], sy[j], csy[j], ccsy[j]);
					compAddSoA(cz[j], sz[j], csz[j], ccsz[j]);
				}
			} else {
				for (int j = 0; j < n; ++j) {
					sx[j] += cx[j];
					sy[j] += cy[j];
					sz[j] += cz[j];
				}
			}
		}

		for (int j = 0; j < n; ++j) {
			ox[i0 + j] += sx[j] + csx[j] + ccsx[j];
			oy[i0 + j] += sy[j] + csy[j] + ccsy[j];
			oz[i0 + j] += sz[j] + csz[j] + ccsz[j];
		}
	}
}

void kernelMagneticFieldPolygonFilament(
		const double *vx, const double *vy, const double *vz, double current,
		const double *ex, const double *ey, const double *ez,
		double *bx, double *by, double *bz,
		int idxSourceStart, int idxSourceEnd, int idxEvalStart, int idxEvalEnd,
		bool useCompensatedSummation) {
	kernelPolygonFilamentSoA<true>(vx, vy, vz, current, ex, ey, ez, bx, by, bz,
			idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd, useCompensatedSummation);
}

void kernelVectorPotentialPolygonFilament(
		const double *vx, const double *vy, const double *vz, double current,
		const double *ex, const double *ey, const double *ez,
		double *ax, double *ay, double *az,
		int idxSourceStart, int idxSourceEnd, int idxEvalStart, int idxEvalEnd,
		bool useCompensatedSummation) {
	kernelPolygonFilamentSoA<false>(vx, vy, vz, current, ex, ey, ez, ax, ay, az,
			idxSourceStart, idxSourceEnd, idxEvalStart, idxEvalEnd, useCompensatedSummation);
}

static void polygonFilamentSoA(
		bool potential, int numVertices,
		const double *vx, const double *vy, const double *vz, double current,
		int numEvalPos, const double *ex, const double *ey, const double *ez,
		double *ox, double *oy, double *oz,
		int numProcessors, bool useCompensatedSummation) {

	if (numVertices < 2) {
		printf("need at least 2 vertices, but only got %d\n", numVertices);
		return;
	}
	if (numProcessors < 1) {
		printf("need at least 1 processor, but only got %d\n", numProcessors);
		return;
	}
	if (current == 0.0) {
		return;
	}

	const int numSegments = numVertices - 1;
	auto kernel = potential ? kernelPolygonFilamentSoA<false> : kernelPolygonFilamentSoA<true>;

	if (numProcessors > 1 && numSegments > numEvalPos) {
		// parallelize over segments, with one output array per thread
		const int nThreads = numSegments < numProcessors ? numSegments : numProcessors;
		const size_t stride = (size_t) numEvalPos;
		double *partial = (double *) calloc(3 * nThreads * stride, sizeof(double));
		if (partial == NULL) {
			printf("failed to allocate temporary array for contributions\n");
			return;
		}
#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
		for (int t = 0; t < nThreads; ++t) {
			double *p = partial + 3 * t * stride;
			kernel(vx, vy, vz, current, ex, ey, ez, p, p + stride, p + 2 * stride,
					(int) ((long) numSegments * t / nThreads), (int) ((long) numSegments * (t + 1) / nThreads),
					0, numEvalPos, useCompensatedSummation);
		}
		for (int i = 0; i < numEvalPos; ++i) {
			double sum[3][3] = {};
			for (int t = 0; t < nThreads; ++t) {
				for (int c = 0; c < 3; ++c) {
					compAdd(partial[(3 * t + c) * stride + i], sum[c]);
				}
			}
			ox[i] += sum[0][0] + sum[0][1] + sum[0][2];
			oy[i] += sum[1][0] + sum[1][1] + sum[1][2];
			oz[i] += sum[2][0] + sum[2][1] + sum[2][2];
		}
		free(partial);
	} else {
		// parallelize over evaluation locations
		const int nThreads = numEvalPos < numProcessors ? numEvalPos : numProcessors;
#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
		for (int t = 0; t < nThreads; ++t) {
			kernel(vx, vy, vz, current, ex, ey, ez, ox, oy, oz, 0, numSegments,
					(int) ((long) numEvalPos * t / nThreads), (int) ((long) numEvalPos * (t + 1) / nThreads),
					useCompensatedSummation);
		}
	}
}

void magneticFieldPolygonFilament(
		int numVertices, const double *vx, const double *vy, const double *vz, double current,
		int numEvalPos, const double *ex, const double *ey, const double *ez,
		double *bx, double *by, double *bz, int numProcessors, bool useCompensatedSummation) {
	polygonFilamentSoA(false, numVertices, vx, vy, vz, current, numEvalPos, ex, ey, ez, bx, by, bz,
			numProcessors, useCompensatedSummation);
}

void vectorPotentialPolygonFilament(
		int numVertices, const double *vx, const double *vy, const double *vz, double current,
		int numEvalPos, const double *ex, const double *ey, const double *ez,
		double *ax, double *ay, double *az, int numProcessors, bool useCompensatedSummation) {
	polygonFilamentSoA(true, numVertices, vx, vy, vz, current, numEvalPos, ex, ey, ez, ax, ay, az,
			numProcessors, useCompensatedSummation);
}

} // namespace abscab
