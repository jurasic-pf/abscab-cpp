#ifndef ABSCAB_H
#define ABSCAB_H

#ifdef _OPENMP
#include <omp.h>
#endif // _OPENMP

namespace abscab {

/** vacuum magnetic permeability in Vs/Am (CODATA-2018) */
const double MU_0 = 1.25663706212e-6;

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
double cel(double k_c, double p, double a, double b);

/////// A_z of straight wire segment

/**
 * Compute the normalized axial component of magnetic vector potential of straight wire segment,
 * evaluated along axis of wire segment (rho = 0).
 * This is a special case for points away from the wire ("far-field") for zP < -1 or zP >= 2.
 *
 * @param zP normalized axial coordinate of evaluation location; must not be in [0, 1] (on wire segment)
 * @return normalized axial component of magnetic vector potential
 */
double sws_A_z_ax_f(double zP);

/**
 * Compute the normalized axial component of magnetic vector potential of straight wire segment,
 * evaluated along axis of wire segment (rhoP = 0).
 * This is a special case for points close to the wire ("near-field") for -1 <= zP < 2.
 *
 * @param zP normalized axial coordinate of evaluation location; must not be in [0, 1] (on wire segment)
 * @return normalized axial component of magnetic vector potential
 */
double sws_A_z_ax_n(double zP);

/**
 * Compute the normalized axial component of magnetic vector potential of straight wire segment,
 * evaluated along axis of wire segment (rho = 0).
 *
 * @param zP normalized axial coordinate of evaluation location; must not be in [0, 1] (on wire segment)
 * @return normalized axial component of magnetic vector potential
 */
double sws_A_z_ax(double zP);

/**
 * Compute the normalized axial component of the magnetic vector potential of a straight wire segment,
 * evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
 * This is a special case for points away from the wire ("far-field") for rhoP > 1.
 *
 * @param rhoP normalized radial coordinate of evaluation location; must not be zero (on wire segment)
 * @return normalized axial component of magnetic vector potential
 */
double sws_A_z_rad_f(double rhoP);

/**
 * Compute the normalized axial component of the magnetic vector potential of a straight wire segment,
 * evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
 * This is a special case for points close to the wire ("near-field") for rhoP <= 1.
 *
 * @param rhoP normalized radial coordinate of evaluation location; must not be zero (on wire segment)
 * @return normalized axial component of magnetic vector potential
 */
double sws_A_z_rad_n(double rhoP);

/**
 * Compute the normalized axial component of the magnetic vector potential of a straight wire segment,
 * evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
 *
 * @param rhoP normalized radial coordinate of evaluation location; must not be zero (on wire segment)
 * @return normalized axial component of magnetic vector potential
 */
double sws_A_z_rad(double rhoP);

/**
 * Compute the normalized axial component of the magnetic vector potential of a straight wire segment.
 * This formulation is useful for points away from the wire ("far-field")
 * at rhoP >= 1 or zP <= -1 or zP > 2.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized axial component of magnetic vector potential
 */
double sws_A_z_f(double rhoP, double zP);

/**
 * Compute the normalized axial component of the magnetic vector potential of a straight wire segment.
 * This formulation is useful for points close to the wire ("near-field")
 * at rhoP < 1 and -1 < zP <= 2.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized axial component of magnetic vector potential
 */
double sws_A_z_n(double rhoP, double zP);

/////// B_phi of straight wire segment

/**
 * Compute the normalized tangential component of the magnetic field of a straight wire segment,
 * evaluated radially along the endpoints of the wire segment (zP = 0 or zP = 1).
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @return normalized tangential component of magnetic field
 */
double sws_B_phi_rad(double rhoP);

/**
 * Compute the normalized tangential component of the magnetic field of a straight wire segment.
 * This formulation is useful for points away from the wire ("far-field")
 * at rhoP >= 1 or zP <= 0 or zP >= 1 or rhoP/(1-zP) >= 1 or rhoP/zP >= 1.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized tangential component of magnetic field
 */
double sws_B_phi_f(double rhoP, double zP);

/**
 * Compute the normalized tangential component of the magnetic field of a straight wire segment.
 * This formulation is useful for points close to the wire ("near-field")
 * at rhoP < 1 and 0 < zP < 1 and rhoP/(1-zP) < 1 and rhoP/zP < 1.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized tangential component of magnetic field
 */
double sws_B_phi_n(double rhoP, double zP);

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
double cwl_A_phi_f(double rhoP, double zP);

/**
 * Compute the normalized tangential component of the magnetic vector potential of a circular wire loop.
 * This formulation is useful for points close to the wire ("near-field")
 * at 1/2 <= rhoP <= 2 and |zP| < 1.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized tangential component of magnetic vector potential
 */
double cwl_A_phi_n(double rhoP, double zP);

/**
 * Compute the normalized tangential component of the magnetic vector potential of a circular wire loop.
 * This formulation is useful for points along rhoP=1 with |zP| < 1.
 *
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized tangential component of magnetic vector potential
 */
double cwl_A_phi_v(double zP);

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
double cwl_B_rho_f(double rhoP, double zP);

/**
 * Compute the normalized radial component of the magnetic field of a circular wire loop.
 * This formulation is useful for points close to the wire ("near-field")
 * at 1/2 <= rhoP <= 2 and |zP| < 1.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized radial component of magnetic field
 */
double cwl_B_rho_n(double rhoP, double zP);

/**
 * Compute the normalized radial component of the magnetic field of a circular wire loop.
 * This formulation is useful for points along rhoP=1 with |zP| < 1.
 *
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized radial component of magnetic field
 */
double cwl_B_rho_v(double zP);

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
double cwl_B_z_f1(double rhoP, double zP);

/**
 * Compute the normalized vertical component of the magnetic field of a circular wire loop.
 * This formulation is useful for certain other points away from the wire ("far-field")
 * at rhoP > 2.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized vertical component of magnetic field
 */
double cwl_B_z_f2(double rhoP, double zP);

/**
 * Compute the normalized vertical component of the magnetic field of a circular wire loop.
 * This formulation is useful for points close to the wire ("near-field")
 * at 1/2 <= rhoP <= 2, but not rhoP=1, and |zP| <= 1.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized vertical component of magnetic field
 */
double cwl_B_z_n(double rhoP, double zP);

/**
 * Compute the normalized vertical component of the magnetic field of a circular wire loop.
 * This formulation is useful for points along rhoP=1 with |zP| <= 1.
 *
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized vertical component of magnetic field
 */
double cwl_B_z_v(double zP);

// --------------------------------------------------

/**
 * Compute the normalized axial component of the magnetic vector potential of a straight wire segment.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized axial component of magnetic vector potential
 */
double straightWireSegment_A_z(double rhoP, double zP);

/**
 * Compute the normalized tangential component of the magnetic field of a straight wire segment.
 *
 * @param rhoP normalized radial coordinate of evaluation location
 * @param zP normalized axial coordinate of evaluation location
 * @return normalized tangential component of magnetic field
 */
double straightWireSegment_B_phi(double rhoP, double zP);

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
double circularWireLoop_A_phi(double rhoP, double zP);

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
double circularWireLoop_B_rho(double rhoP, double zP);

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
double circularWireLoop_B_z(double rhoP, double zP);

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
		double current, int nEvalPos, double *evalPos, double *vectorPotential);

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
		double current, int nEvalPos, double *evalPos, double *magneticField);

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
		bool useCompensatedSummation);

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
		bool useCompensatedSummation);

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
		bool useCompensatedSummation);

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
		bool useCompensatedSummation);

// --------------------------------------------------

/**
 * Compute the magnetic vector potential of a polygon filament
 * at a number of evaluation locations.
 *
 * @param vertices [3: x, y, z][numVertices] points along polygon; in m
 * @param current current along polygon; in A
 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
 * @param vectorPotential [3: x, y, z][numEvalPos] target array for magnetic vector potential at evaluation locations; in Tm
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
#ifdef _OPENMP
		int numProcessors = omp_get_max_threads(),
#else
		int numProcessors = 1,
#endif
		bool useCompensatedSummation = true);

/**
 * Compute the magnetic vector potential of a polygon filament
 * at a number of evaluation locations.
 *
 * @param void (*vertexSupplier)(int i, double *point): callback to put i-th current carrier polygon vertex into point as [3: x, y, z]; in m
 * @param current current along polygon; in A
 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
 * @param vectorPotential [3: x, y, z][numEvalPos] target array for magnetic vector potential at evaluation locations; in Tm
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
#ifdef _OPENMP
		int numProcessors = omp_get_max_threads(),
#else
		int numProcessors = 1,
#endif
		bool useCompensatedSummation = true);

/**
 * Compute the magnetic field of a polygon filament
 * at a number of evaluation locations.
 *
 * @param vertices [3: x, y, z][numVertices] points along polygon; in m
 * @param current current along polygon; in A
 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
 * @param magneticField [3: x, y, z][numEvalPos] target array for magnetic field at evaluation locations; inT
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
#ifdef _OPENMP
		int numProcessors = omp_get_max_threads(),
#else
		int numProcessors = 1,
#endif
		bool useCompensatedSummation = true);


/**
 * Compute the magnetic field of a polygon filament
 * at a number of evaluation locations.
 *
 * @param void (*vertexSupplier)(int i, double *point): callback to put i-th current carrier polygon vertex into point as [3: x, y, z]; in m
 * @param current current along polygon; in A
 * @param evalPos [3: x, y, z][numEvalPos] evaluation locations; in m
 * @param magneticField [3: x, y, z][numEvalPos] target array for magnetic field at evaluation locations; in T
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
#ifdef _OPENMP
		int numProcessors = omp_get_max_threads(),
#else
		int numProcessors = 1,
#endif
		bool useCompensatedSummation = true);

/**
 * Struct-of-arrays variants of the polygon filament routines above: each coordinate
 * is a separate array, [numVertices] for vx, vy, vz and [numEvalPos] for the rest.
 * Results are added to the targets (bx, by, bz or ax, ay, az).
 * Generic far-field locations are vectorized and use sqrt(a*a + b*b) instead of hypot(a, b);
 * all others use the scalar routines.
 */
void magneticFieldPolygonFilament(
		int numVertices,
		const double *vx, const double *vy, const double *vz,
		double current,
		int numEvalPos,
		const double *ex, const double *ey, const double *ez,
		double *bx, double *by, double *bz,
#ifdef _OPENMP
		int numProcessors = omp_get_max_threads(),
#else
		int numProcessors = 1,
#endif
		bool useCompensatedSummation = true);

void vectorPotentialPolygonFilament(
		int numVertices,
		const double *vx, const double *vy, const double *vz,
		double current,
		int numEvalPos,
		const double *ex, const double *ey, const double *ez,
		double *ax, double *ay, double *az,
#ifdef _OPENMP
		int numProcessors = omp_get_max_threads(),
#else
		int numProcessors = 1,
#endif
		bool useCompensatedSummation = true);

/**
 * Struct-of-arrays kernels of the routines above for the segments [idxSourceStart, idxSourceEnd)
 * and the evaluation locations [idxEvalStart, idxEvalEnd); no parallelization.
 */
void kernelMagneticFieldPolygonFilament(
		const double *vx, const double *vy, const double *vz,
		double current,
		const double *ex, const double *ey, const double *ez,
		double *bx, double *by, double *bz,
		int idxSourceStart,
		int idxSourceEnd,
		int idxEvalStart,
		int idxEvalEnd,
		bool useCompensatedSummation);

void kernelVectorPotentialPolygonFilament(
		const double *vx, const double *vy, const double *vz,
		double current,
		const double *ex, const double *ey, const double *ez,
		double *ax, double *ay, double *az,
		int idxSourceStart,
		int idxSourceEnd,
		int idxEvalStart,
		int idxEvalEnd,
		bool useCompensatedSummation);

} // namespace abscab

#endif // ABSCAB_H
