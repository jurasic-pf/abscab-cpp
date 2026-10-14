
#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "abscab.hh"

#include "abscab-cpp/test_util/test_util.hh"

namespace abscab {

/** test case for cel() implementation as described in section 4.2 of the 1969 Bulirsch article */
TEST(TestCel, CheckCel) {
	double tolerance = 1.0e-15;

	double k_c = 0.1;
	double p1 =  4.1;
	double p2 = -4.1;
	double a = 1.2;
	double b = 1.1;

	double cel1 =  1.5464442694017956;
	double cel2 = -6.7687378198360556e-1;

	double c1 = cel(k_c, p1, a, b);
	double c2 = cel(k_c, p2, a, b);

//	double ra1 = fabs(cel1 - c1)/(1.0 + fabs(cel1));
//	double ra2 = fabs(cel2 - c2)/(1.0 + fabs(cel2));
//	printf("case 1: rel/abs deviation = %g\n", ra1);
//	printf("case 2: rel/abs deviation = %g\n", ra2);

	EXPECT_EQ(assertRelAbsEquals(cel1, c1, tolerance), 0);
	EXPECT_EQ(assertRelAbsEquals(cel2, c2, tolerance), 0);
}

/** test method for Straight Wire Segment methods */
TEST(TestAbscab, CheckStraightWireSegment) {

	int rowsRp = 0;
	int colsRp = 0;
	double **test_points_rp = loadColumnsFromFile("abscab-cpp/abscab/resources/testPointsRpStraightWireSegment.dat", &rowsRp, &colsRp);
	if (rowsRp < 1) {
		printf("error: need at least one row of test point coordinates for rho'\n");
		FAIL();
	}
	if (colsRp < 1) {
		printf("error: need at least one column of test point coordinates for rho'\n");
		FAIL();
	}

	int rowsZp = 0;
	int colsZp = 0;
	double **test_points_zp = loadColumnsFromFile("abscab-cpp/abscab/resources/testPointsZpStraightWireSegment.dat", &rowsZp, &colsZp);
	if (rowsZp < 1) {
		printf("error: need at least one row of test point coordinates for z'\n");
		FAIL();
	}
	if (colsZp < 1) {
		printf("error: need at least one column of test point coordinates for z'\n");
		FAIL();
	}

	if (rowsRp != rowsZp) {
		printf("error: number of rows of test point coordinates has to agree between rho' (%d) and z' (%d)\n", rowsRp, rowsZp);
		FAIL();
	}
	if (colsRp != colsZp) {
		printf("error: number of columns of test point coordinates has to agree between rho' (%d) and z' (%d)\n", colsRp, colsZp);
		FAIL();
	}

	int numCases = rowsRp;

	int rowsAZRef = 0;
	int colsAZRef = 0;
	double **aZRef = loadColumnsFromFile("abscab-cpp/abscab/resources/StraightWireSegment_A_z_ref.dat", &rowsAZRef, &colsAZRef);
	if (rowsAZRef < 1) {
		printf("error: need at least one row of reference values for A_z\n");
		FAIL();
	}
	if (colsAZRef < 1) {
		printf("error: need at least one column of reference values for A_z\n");
		FAIL();
	}

	if (numCases != rowsAZRef) {
		printf("error: number of reference values for A_z (%d) has to match number of test cases (%d)\n", rowsAZRef, rowsRp);
		FAIL();
	}

	int rowsBPhiRef = 0;
	int colsBPhiRef = 0;
	double **bPhiRef = loadColumnsFromFile("abscab-cpp/abscab/resources/StraightWireSegment_B_phi_ref.dat", &rowsBPhiRef, &colsBPhiRef);
	if (rowsBPhiRef < 1) {
		printf("error: need at least one row of reference values for B_phi\n");
		FAIL();
	}
	if (colsBPhiRef < 1) {
		printf("error: need at least one column of reference values for B_phi\n");
		FAIL();
	}

	if (numCases != rowsBPhiRef) {
		printf("error: number of reference values for B_phi (%d) has to match number of test cases (%d)\n", rowsBPhiRef, rowsRp);
		FAIL();
	}

	double toleranceAZ   = 1.0e-15;
	double toleranceBPhi = 1.0e-15;

	for (int i = 0; i < numCases; ++i) {
		int status = 0;

		double rp = test_points_rp[0][i];
		double zp = test_points_zp[0][i];

		// compute values using C implementation to test
		double aZ   = straightWireSegment_A_z(rp, zp);
		double bPhi = straightWireSegment_B_phi(rp, zp);

		int aZStatus = assertRelAbsEquals(aZRef[0][i], aZ, toleranceAZ);
		EXPECT_EQ(aZStatus, 0);
		if (aZStatus) {
			printf("error: mismatch at Straight Wire Segment A_z test case %d\n", i);
			printf("  rho' = %.17e\n", rp);
			printf("    z' = %.17e\n", zp);
			printf("  ref A_z = %+.17e\n", aZRef[0][i]);
			printf("  act A_z = %+.17e\n", aZ);
		}
		status |= aZStatus;

		int bPhiStatus = assertRelAbsEquals(bPhiRef[0][i], bPhi, toleranceBPhi);
		EXPECT_EQ(aZStatus, 0);
		if (bPhiStatus) {
			printf("error: mismatch at Straight Wire Segment B_phi test case %d\n", i);
			printf("  rho' = %.17e\n", rp);
			printf("    z' = %.17e\n", zp);
			printf("  ref B_phi = %+.17e\n", bPhiRef[0][i]);
			printf("  act B_phi = %+.17e\n", bPhi);
		}
		status |= bPhiStatus;

		if (status) {
			break;
		}
	}

	// free data loaded from text files
	free(test_points_rp[0]); free(test_points_rp);
	free(test_points_zp[0]); free(test_points_zp);
	free(aZRef[0]);          free(aZRef);
	free(bPhiRef[0]);        free(bPhiRef);
}

/** test method for Circular Wire Loop methods */
TEST(TestAbscab, CheckCircularWireLoop) {

	int rowsRp = 0;
	int colsRp = 0;
	double **test_points_rp = loadColumnsFromFile("abscab-cpp/abscab/resources/testPointsRpCircularWireLoop.dat", &rowsRp, &colsRp);
	if (rowsRp < 1) {
		printf("error: need at least one row of test point coordinates for rho'\n");
		FAIL();
	}
	if (colsRp < 1) {
		printf("error: need at least one column of test point coordinates for rho'\n");
		FAIL();
	}

	int rowsZp = 0;
	int colsZp = 0;
	double **test_points_zp = loadColumnsFromFile("abscab-cpp/abscab/resources/testPointsZpCircularWireLoop.dat", &rowsZp, &colsZp);
	if (rowsZp < 1) {
		printf("error: need at least one row of test point coordinates for z'\n");
		FAIL();
	}
	if (colsZp < 1) {
		printf("error: need at least one column of test point coordinates for z'\n");
		FAIL();
	}

	if (rowsRp != rowsZp) {
		printf("error: number of rows of test point coordinates has to agree between rho' (%d) and z' (%d)\n", rowsRp, rowsZp);
		FAIL();
	}
	if (colsRp != colsZp) {
		printf("error: number of columns of test point coordinates has to agree between rho' (%d) and z' (%d)\n", colsRp, colsZp);
		FAIL();
	}

	int numCases = rowsRp;

	int rowsAPhiRef = 0;
	int colsAPhiRef = 0;
	double **aPhiRef = loadColumnsFromFile("abscab-cpp/abscab/resources/CircularWireLoop_A_phi_ref.dat", &rowsAPhiRef, &colsAPhiRef);
	if (rowsAPhiRef < 1) {
		printf("error: need at least one row of reference values for A_phi\n");
		FAIL();
	}
	if (colsAPhiRef < 1) {
		printf("error: need at least one column of reference values for A_phi\n");
		FAIL();
	}

	if (numCases != rowsAPhiRef) {
		printf("error: number of reference values for A_phi (%d) has to match number of test cases (%d)\n", rowsAPhiRef, rowsRp);
		FAIL();
	}

	int rowsBRhoRef = 0;
	int colsBRhoRef = 0;
	double **bRhoRef = loadColumnsFromFile("abscab-cpp/abscab/resources/CircularWireLoop_B_rho_ref.dat", &rowsBRhoRef, &colsBRhoRef);
	if (rowsBRhoRef < 1) {
		printf("error: need at least one row of reference values for B_rho\n");
		FAIL();
	}
	if (colsBRhoRef < 1) {
		printf("error: need at least one column of reference values for B_rho\n");
		FAIL();
	}

	if (numCases != rowsBRhoRef) {
		printf("error: number of reference values for B_rho (%d) has to match number of test cases (%d)\n", rowsBRhoRef, rowsRp);
		FAIL();
	}

	int rowsBZRef = 0;
	int colsBZRef = 0;
	double **bZRef = loadColumnsFromFile("abscab-cpp/abscab/resources/CircularWireLoop_B_z_ref.dat", &rowsBZRef, &colsBZRef);
	if (rowsBZRef < 1) {
		printf("error: need at least one row of reference values for B_z\n");
		FAIL();
	}
	if (colsBZRef < 1) {
		printf("error: need at least one column of reference values for B_z\n");
		FAIL();
	}

	if (numCases != rowsBZRef) {
		printf("error: number of reference values for B_z (%d) has to match number of test cases (%d)\n", rowsBZRef, rowsRp);
		FAIL();
	}

	double toleranceAPhi = 1.0e-15;
	double toleranceBRho = 1.0e-13;
	double toleranceBZ   = 1.0e-14;

	for (int i = 0; i < numCases; ++i) {
		int status = 0;

		double rp = test_points_rp[0][i];
		double zp = test_points_zp[0][i];

		// compute values using C implementation to test
		double aPhi = circularWireLoop_A_phi(rp, zp);
		double bRho = circularWireLoop_B_rho(rp, zp);
		double bZ   = circularWireLoop_B_z(rp, zp);

		int aPhiStatus = assertRelAbsEquals(aPhiRef[0][i], aPhi, toleranceAPhi);
		EXPECT_EQ(aPhiStatus, 0);
		if (aPhiStatus) {
			printf("error: mismatch at Circular Wire Loop A_phi test case %d\n", i);
			printf("  rho' = %.17e\n", rp);
			printf("    z' = %.17e\n", zp);
			printf("  ref A_phi = %+.17e\n", aPhiRef[0][i]);
			printf("  act A_phi = %+.17e\n", aPhi);
		}
		status |= aPhiStatus;

		int bRhoStatus = assertRelAbsEquals(bRhoRef[0][i], bRho, toleranceBRho);
		EXPECT_EQ(bRhoStatus, 0);
		if (bRhoStatus) {
			printf("error: mismatch at Circular Wire Loop B_rho test case %d\n", i);
			printf("  rho' = %.17e\n", rp);
			printf("    z' = %.17e\n", zp);
			printf("  ref B_rho = %+.17e\n", bRhoRef[0][i]);
			printf("  act B_rho = %+.17e\n", bRho);
		}
		status |= bRhoStatus;

		int bZStatus = assertRelAbsEquals(bZRef[0][i], bZ, toleranceBZ);
		EXPECT_EQ(bZStatus, 0);
		if (bZStatus) {
			printf("error: mismatch at Circular Wire Loop B_z test case %d\n", i);
			printf("  rho' = %.17e\n", rp);
			printf("    z' = %.17e\n", zp);
			printf("  ref B_z = %+.17e\n", bZRef[0][i]);
			printf("  act B_z = %+.17e\n", bZ);
		}
		status |= bZStatus;
	}

	// free data loaded from text files
	free(test_points_rp[0]); free(test_points_rp);
	free(test_points_zp[0]); free(test_points_zp);
	free(aPhiRef[0]);        free(aPhiRef);
	free(bRhoRef[0]);        free(bRhoRef);
	free(bZRef[0]);          free(bZRef);
}

TEST(TestAbscab, CheckMagneticFieldInfiniteLineFilament) {
	double tolerance = 1.0e-15;

	// Demtroeder 2, Sec. 3.2.2 ("Magnetic field of a straight wire")
	// B(r) = mu_0 * I / (2 pi r)
	// Test this here with:
	// I = 123.0 A
	// r = 0.132 m
	// => B = 0.186 mT
	double current = 123.0;
	double r = 0.132;
	double bPhiRef = MU_0 * current / (2.0 * M_PI * r);
//	printf("ref bPhi = %.5e\n", bPhiRef);

	double vertices[] = {
			0.0, 0.0, -1.0e6,
			0.0, 0.0,  1.0e6
	};

	double evalPos[] = {
			r, 0.0, 0.0
	};

	// y component is B_phi
	double magneticField[3];
	magneticFieldPolygonFilament(2, vertices, current, 1, evalPos, magneticField);
	double bPhi = magneticField[1];
//	printf("act bPhi = %.5e\n", bPhi);

//	double relAbsErr = fabs(bPhi - bPhiRef) / (1.0 + fabs(bPhiRef));
//	printf("raErr = %.5e\n", relAbsErr);

	EXPECT_EQ(assertRelAbsEquals(bPhiRef, bPhi, tolerance), 0);
}

TEST(TestAbscab, CheckBPhiInfiniteLineFilament) {
	double tolerance = 1.0e-15;

	// Demtroeder 2, Sec. 3.2.2 ("Magnetic field of a straight wire")
	// B(r) = mu_0 * I / (2 pi r)
	// Test this here with:
	// I = 123.0 A
	// r = 0.132 m
	// => B = 0.186 mT
	double current = 123.0;
	double r = 0.132;
	double bPhiRef = MU_0 * current / (2.0 * M_PI * r);
//	printf("ref bPhi = %.5e\n", bPhiRef);

	// half the length of the wire segment
	double halfL = 1e6;
	double L = 2*halfL;
	double rhoP = r / L;
	double zP = halfL / L;
	double bPhi = MU_0 * current / (4.0 * M_PI * L) * straightWireSegment_B_phi(rhoP, zP);
//	printf("act bPhi = %.5e\n", bPhi);

//	double relAbsErr = fabs(bPhi - bPhiRef) / (1.0 + fabs(bPhiRef));
//	printf("raErr = %.5e\n", relAbsErr);

	EXPECT_EQ(assertRelAbsEquals(bPhiRef, bPhi, tolerance), 0);
}

TEST(TestAbscab, CheckMagneticFieldInsideLongCoil) {
	double tolerance = 1.0e-4;

	// Demtroeder 2, Sec. 3.2.3 ("Magnetic field of a long coil")
	// B_z = mu_0 * n * I
	// where n is the winding density: n = N / L
	// of a coil of N windings over a length L
	// Example (which is tested here):
	// n = 1e3 m^{-1}
	// I = 10 A
	// => B = 0.0126T
	double bZRef = 0.0126;

	int N = 50000; // windings
	double L = 50.0; // total length of coil in m
	double n = N/L;

	double current = 10.0; // A
	double radius = 1.0; // m

	double bZ = 0.0;
	for (int i = 0; i < N; ++i) {

		// axial position of coil
		double z0 = -L/2.0 + (i + 0.5) / n;

		// compute magnetic field
		double prefac = MU_0 * current / (M_PI * radius);
		double bZContrib = prefac * circularWireLoop_B_z(0.0, z0);

//		printf("coil %d at z0 = % .3e => contrib = %.3e\n", i, z0, bZContrib);

		bZ += bZContrib;
	}
//	printf("B_z = %.5e\n", bZ);

//	double relAbsErr = fabs(bZ - bZRef) / (1.0 + fabs(bZRef));
//	printf("raErr = %.5e\n", relAbsErr);

	EXPECT_EQ(assertRelAbsEquals(bZRef, bZ, tolerance), 0);
}

// check that previous field contributions are preserved
// for circular wire loop methods
TEST(TestAbscab, CheckPreseveranceCircularWireLoop) {

	double tolerance = 1.0e-15;

	double center[3] = { 1.0, 2.0, 3.0 };
	double normal[3] = { 4.0, 5.0, 6.0 };
	double radius = 7.89;

	double current_1 = 12.3;
	double current_2 = 45.6;

	double evalPos[3] = { 7.0, 8.0, 9.0 };

	// reference: add two contributions
	double vectorPotential_ref[3] = { 0.0, 0.0, 0.0 };
	vectorPotentialCircularFilament(center, normal, radius, current_1 + current_2, 1, evalPos, vectorPotential_ref);

	// to be tested: add 
	double vectorPotential[3] = { 0.0, 0.0, 0.0 };
	vectorPotentialCircularFilament(center, normal, radius, current_1, 1, evalPos, vectorPotential);
	vectorPotentialCircularFilament(center, normal, radius, current_2, 1, evalPos, vectorPotential);

	EXPECT_EQ(assertRelAbsEquals(vectorPotential_ref[0], vectorPotential[0], tolerance), 0);
	EXPECT_EQ(assertRelAbsEquals(vectorPotential_ref[1], vectorPotential[1], tolerance), 0);
	EXPECT_EQ(assertRelAbsEquals(vectorPotential_ref[2], vectorPotential[2], tolerance), 0);

	// and again for the magnetic field ...
	double magneticField_ref[3] = { 0.0, 0.0, 0.0 };
	magneticFieldCircularFilament(center, normal, radius, current_1 + current_2, 1, evalPos, magneticField_ref);

	double magneticField[3] = { 0.0, 0.0, 0.0 };
	magneticFieldCircularFilament(center, normal, radius, current_1, 1, evalPos, magneticField);
	magneticFieldCircularFilament(center, normal, radius, current_2, 1, evalPos, magneticField);

	EXPECT_EQ(assertRelAbsEquals(magneticField_ref[0], magneticField[0], tolerance), 0);
	EXPECT_EQ(assertRelAbsEquals(magneticField_ref[1], magneticField[1], tolerance), 0);
	EXPECT_EQ(assertRelAbsEquals(magneticField_ref[2], magneticField[2], tolerance), 0);
}

TEST(TestAbscab, CheckPreseverancePolygonFilamentArrays) {

	double tolerance = 1.0e-15;

	constexpr int numVertices = 4;
	double vertices[numVertices * 3] = {
		0.0, 0.0, 0.0,
		1.0, 2.0, 3.0,
		4.0, 5.0, 6.0,
		7.0, 8.0, 9.0
	};

	double current_1 = 1.23;
	double current_2 = 4.56;

	for (int numProcessors = 1; numProcessors < 3; ++numProcessors) {
		for (int compensatedSummationCase = 0; compensatedSummationCase < 2; ++compensatedSummationCase) {
			bool useCompensatedSummation = (compensatedSummationCase == 0);
			
			// trigger both cases of numEvalPos < numSegments and numEvalPos > numSegments
			for (int numEvalPos = 1; numEvalPos < 6; numEvalPos += 4) {

				int numBytes = numEvalPos * 3 * sizeof(double);

				double *evalPos = (double *) malloc(numBytes);
				if (evalPos == NULL) {
					FAIL() << "could not allocate evalPos";
				}
				for (int idxEvalPos = 0; idxEvalPos < numEvalPos; ++idxEvalPos) {
					evalPos[idxEvalPos * 3 + 0] = -(idxEvalPos + 1);
					evalPos[idxEvalPos * 3 + 1] = -(idxEvalPos + 1);
					evalPos[idxEvalPos * 3 + 2] = -(idxEvalPos + 1);
				}

				// vector potential
				double *vectorPotential_ref = (double *) malloc(numBytes);
				if (vectorPotential_ref == NULL) {
					FAIL() << "could not allocate vectorPotential_ref";
				}
				memset(vectorPotential_ref, 0, numBytes);

				vectorPotentialPolygonFilament(numVertices, vertices, current_1 + current_2, numEvalPos, evalPos, vectorPotential_ref, numProcessors, useCompensatedSummation);

				double *vectorPotential = (double *) malloc(numBytes);
				if (vectorPotential == NULL) {
					FAIL() << "could not allocate vectorPotential";
				}
				memset(vectorPotential, 0, numBytes);

				vectorPotentialPolygonFilament(numVertices, vertices, current_1, numEvalPos, evalPos, vectorPotential, numProcessors, useCompensatedSummation);
				vectorPotentialPolygonFilament(numVertices, vertices, current_2, numEvalPos, evalPos, vectorPotential, numProcessors, useCompensatedSummation);

				for (int idxEvalPos = 0; idxEvalPos < numEvalPos; ++idxEvalPos) {
					EXPECT_EQ(assertRelAbsEquals(vectorPotential_ref[idxEvalPos * 3 + 0], vectorPotential[idxEvalPos * 3 + 0], tolerance), 0);
					EXPECT_EQ(assertRelAbsEquals(vectorPotential_ref[idxEvalPos * 3 + 1], vectorPotential[idxEvalPos * 3 + 1], tolerance), 0);
					EXPECT_EQ(assertRelAbsEquals(vectorPotential_ref[idxEvalPos * 3 + 2], vectorPotential[idxEvalPos * 3 + 2], tolerance), 0);
				}

				free(vectorPotential);
				free(vectorPotential_ref);

				// magnetic field
				double *magneticField_ref = (double *) malloc(numBytes);
				if (magneticField_ref == NULL) {
					FAIL() << "could not allocate magneticField_ref";
				}
				memset(magneticField_ref, 0, numBytes);

				magneticFieldPolygonFilament(numVertices, vertices, current_1 + current_2, numEvalPos, evalPos, magneticField_ref, numProcessors, useCompensatedSummation);

				double *magneticField = (double *) malloc(numBytes);
				if (magneticField == NULL) {
					FAIL() << "could not allocate magneticField";
				}
				memset(magneticField, 0, numBytes);

				magneticFieldPolygonFilament(numVertices, vertices, current_1, numEvalPos, evalPos, magneticField, numProcessors, useCompensatedSummation);
				magneticFieldPolygonFilament(numVertices, vertices, current_2, numEvalPos, evalPos, magneticField, numProcessors, useCompensatedSummation);

				for (int idxEvalPos = 0; idxEvalPos < numEvalPos; ++idxEvalPos) {
					EXPECT_EQ(assertRelAbsEquals(magneticField_ref[idxEvalPos * 3 + 0], magneticField[idxEvalPos * 3 + 0], tolerance), 0);
					EXPECT_EQ(assertRelAbsEquals(magneticField_ref[idxEvalPos * 3 + 1], magneticField[idxEvalPos * 3 + 1], tolerance), 0);
					EXPECT_EQ(assertRelAbsEquals(magneticField_ref[idxEvalPos * 3 + 2], magneticField[idxEvalPos * 3 + 2], tolerance), 0);
				}

				free(magneticField);
				free(magneticField_ref);

				free(evalPos);
			}
		}
	}
}

TEST(TestAbscab, CheckPreseverancePolygonFilamentVertexSupplier) {

	double tolerance = 1.0e-15;

	constexpr int numVertices = 4;
	void (*vertexSupplier)(int, double *) = [](int i, double *point) {
		switch (i) {
		case 0:
			point[0] = 0.0; point[1] = 0.0; point[2] = 0.0;
			break;
		case 1:
			point[0] = 1.0; point[1] = 2.0; point[2] = 3.0;
			break;
		case 2:
			point[0] = 4.0; point[1] = 5.0; point[2] = 6.0;
			break;
		case 3:
			point[0] = 7.0; point[1] = 8.0; point[2] = 9.0;
			break;
		default:
			std::stringstream msg;
			msg << "vertexSupplier only provides " << numVertices << " vertices";
			FAIL() << msg.str();
			break;
		}
	};

	double current_1 = 1.23;
	double current_2 = 4.56;

	for (int numProcessors = 1; numProcessors < 3; ++numProcessors) {
		for (int compensatedSummationCase = 0; compensatedSummationCase < 2; ++compensatedSummationCase) {
			bool useCompensatedSummation = (compensatedSummationCase == 0);

			// trigger both cases of numEvalPos < numSegments and numEvalPos > numSegments
			for (int numEvalPos = 1; numEvalPos < 6; numEvalPos += 4) {

				int numBytes = numEvalPos * 3 * sizeof(double);

				double *evalPos = (double *) malloc(numBytes);
				if (evalPos == NULL) {
					FAIL() << "could not allocate evalPos";
				}
				for (int idxEvalPos = 0; idxEvalPos < numEvalPos; ++idxEvalPos) {
					evalPos[idxEvalPos * 3 + 0] = -(idxEvalPos + 1);
					evalPos[idxEvalPos * 3 + 1] = -(idxEvalPos + 1);
					evalPos[idxEvalPos * 3 + 2] = -(idxEvalPos + 1);
				}

				// vector potential
				double *vectorPotential_ref = (double *) malloc(numBytes);
				if (vectorPotential_ref == NULL) {
					FAIL() << "could not allocate vectorPotential_ref";
				}
				memset(vectorPotential_ref, 0, numBytes);

				vectorPotentialPolygonFilament(numVertices, vertexSupplier, current_1 + current_2, numEvalPos, evalPos, vectorPotential_ref, numProcessors, useCompensatedSummation);

				double *vectorPotential = (double *) malloc(numBytes);
				if (vectorPotential == NULL) {
					FAIL() << "could not allocate vectorPotential";
				}
				memset(vectorPotential, 0, numBytes);

				vectorPotentialPolygonFilament(numVertices, vertexSupplier, current_1, numEvalPos, evalPos, vectorPotential, numProcessors, useCompensatedSummation);
				vectorPotentialPolygonFilament(numVertices, vertexSupplier, current_2, numEvalPos, evalPos, vectorPotential, numProcessors, useCompensatedSummation);

				for (int idxEvalPos = 0; idxEvalPos < numEvalPos; ++idxEvalPos) {
					EXPECT_EQ(assertRelAbsEquals(vectorPotential_ref[idxEvalPos * 3 + 0], vectorPotential[idxEvalPos * 3 + 0], tolerance), 0);
					EXPECT_EQ(assertRelAbsEquals(vectorPotential_ref[idxEvalPos * 3 + 1], vectorPotential[idxEvalPos * 3 + 1], tolerance), 0);
					EXPECT_EQ(assertRelAbsEquals(vectorPotential_ref[idxEvalPos * 3 + 2], vectorPotential[idxEvalPos * 3 + 2], tolerance), 0);
				}

				free(vectorPotential);
				free(vectorPotential_ref);

				// magnetic field
				double *magneticField_ref = (double *) malloc(numBytes);
				if (magneticField_ref == NULL) {
					FAIL() << "could not allocate magneticField_ref";
				}
				memset(magneticField_ref, 0, numBytes);

				magneticFieldPolygonFilament(numVertices, vertexSupplier, current_1 + current_2, numEvalPos, evalPos, magneticField_ref, numProcessors, useCompensatedSummation);

				double *magneticField = (double *) malloc(numBytes);
				if (magneticField == NULL) {
					FAIL() << "could not allocate magneticField";
				}
				memset(magneticField, 0, numBytes);

				magneticFieldPolygonFilament(numVertices, vertexSupplier, current_1, numEvalPos, evalPos, magneticField, numProcessors, useCompensatedSummation);
				magneticFieldPolygonFilament(numVertices, vertexSupplier, current_2, numEvalPos, evalPos, magneticField, numProcessors, useCompensatedSummation);

				for (int idxEvalPos = 0; idxEvalPos < numEvalPos; ++idxEvalPos) {
					EXPECT_EQ(assertRelAbsEquals(magneticField_ref[idxEvalPos * 3 + 0], magneticField[idxEvalPos * 3 + 0], tolerance), 0);
					EXPECT_EQ(assertRelAbsEquals(magneticField_ref[idxEvalPos * 3 + 1], magneticField[idxEvalPos * 3 + 1], tolerance), 0);
					EXPECT_EQ(assertRelAbsEquals(magneticField_ref[idxEvalPos * 3 + 2], magneticField[idxEvalPos * 3 + 2], tolerance), 0);
				}

				free(magneticField);
				free(magneticField_ref);

				free(evalPos);
			}
		}
	}
}



/////// struct-of-arrays polygon filament routines

namespace {

typedef std::vector<double> Vec;

// SoA output of both routines for a polygon given as (x, y, z) arrays; bx.. and ax.. start at zero
struct SoAResult { Vec bx, by, bz, ax, ay, az; };

SoAResult runSoA(const Vec &vx, const Vec &vy, const Vec &vz, const Vec &ex, const Vec &ey, const Vec &ez,
		int numProcessors = 1, bool compensated = true) {
	int ne = (int) ex.size();
	SoAResult r = {Vec(ne), Vec(ne), Vec(ne), Vec(ne), Vec(ne), Vec(ne)};
	magneticFieldPolygonFilament((int) vx.size(), vx.data(), vy.data(), vz.data(), 1.0, ne,
			ex.data(), ey.data(), ez.data(), r.bx.data(), r.by.data(), r.bz.data(), numProcessors, compensated);
	vectorPotentialPolygonFilament((int) vx.size(), vx.data(), vy.data(), vz.data(), 1.0, ne,
			ex.data(), ey.data(), ez.data(), r.ax.data(), r.ay.data(), r.az.data(), numProcessors, compensated);
	return r;
}

// largest |a - b| / |b| over the locations (all three components together)
double maxRelDiff(const Vec &ax, const Vec &ay, const Vec &az, const Vec &bx, const Vec &by, const Vec &bz) {
	double worst = 0.0;
	for (size_t i = 0; i < ax.size(); ++i) {
		double n = std::sqrt(bx[i] * bx[i] + by[i] * by[i] + bz[i] * bz[i]);
		double d = std::sqrt(pow(ax[i] - bx[i], 2) + pow(ay[i] - by[i], 2) + pow(az[i] - bz[i], 2));
		worst = std::max(worst, n > 0.0 ? d / n : d);
	}
	return worst;
}

// scalar AoS kernels (hypot-based) as the reference
SoAResult runScalar(const Vec &vx, const Vec &vy, const Vec &vz, const Vec &ex, const Vec &ey, const Vec &ez, bool compensated = true) {
	int nv = (int) vx.size(), ne = (int) ex.size();
	Vec v(3 * nv), e(3 * ne), b(3 * ne, 0.0), a(3 * ne, 0.0);
	for (int i = 0; i < nv; ++i) { v[3 * i] = vx[i]; v[3 * i + 1] = vy[i]; v[3 * i + 2] = vz[i]; }
	for (int i = 0; i < ne; ++i) { e[3 * i] = ex[i]; e[3 * i + 1] = ey[i]; e[3 * i + 2] = ez[i]; }
	kernelMagneticFieldPolygonFilament(v.data(), 1.0, e.data(), b.data(), 0, nv - 1, 0, ne, compensated);
	kernelVectorPotentialPolygonFilament(v.data(), 1.0, e.data(), a.data(), 0, nv - 1, 0, ne, compensated);
	SoAResult r = {Vec(ne), Vec(ne), Vec(ne), Vec(ne), Vec(ne), Vec(ne)};
	for (int i = 0; i < ne; ++i) {
		r.bx[i] = b[3 * i]; r.by[i] = b[3 * i + 1]; r.bz[i] = b[3 * i + 2];
		r.ax[i] = a[3 * i]; r.ay[i] = a[3 * i + 1]; r.az[i] = a[3 * i + 2];
	}
	return r;
}

} // namespace

// the reference data of CheckStraightWireSegment, evaluated through the SoA polygon routines
TEST(TestAbscab, CheckStraightWireSegmentSoA) {
	int rows, cols;
	double **rp = loadColumnsFromFile("abscab-cpp/abscab/resources/testPointsRpStraightWireSegment.dat", &rows, &cols);
	double **zp = loadColumnsFromFile("abscab-cpp/abscab/resources/testPointsZpStraightWireSegment.dat", &rows, &cols);
	double **aZRef = loadColumnsFromFile("abscab-cpp/abscab/resources/StraightWireSegment_A_z_ref.dat", &rows, &cols);
	double **bPhiRef = loadColumnsFromFile("abscab-cpp/abscab/resources/StraightWireSegment_B_phi_ref.dat", &rows, &cols);

	// unit segment along x: z' is x and rho' is y, so A_z is A_x and B_phi is B_z
	Vec ex(zp[0], zp[0] + rows), ey(rp[0], rp[0] + rows), ez(rows, 0.0);
	SoAResult r = runSoA({0.0, 1.0}, {0.0, 0.0}, {0.0, 0.0}, ex, ey, ez);

	for (int i = 0; i < rows; ++i) {
		double a = r.ax[i] / (MU_0 / (2.0 * M_PI));
		EXPECT_LT(fabs(a - aZRef[0][i]), 1.0e-15 * fabs(aZRef[0][i]) + 1.0e-300) << "A_z rho' " << ey[i] << " z' " << ex[i];
		if (ey[i] != 0.0) {
			double b = r.bz[i] / (MU_0 / (4.0 * M_PI));
			EXPECT_LT(fabs(b - bPhiRef[0][i]), 1.0e-15 * fabs(bPhiRef[0][i]) + 1.0e-300) << "B_phi rho' " << ey[i] << " z' " << ex[i];
		}
	}
	free(rp[0]); free(rp); free(zp[0]); free(zp); free(aZRef[0]); free(aZRef); free(bPhiRef[0]); free(bPhiRef);
}

// all regimes of a coil: generic, on the axis, end planes, near field, and the different parallelizations
TEST(TestAbscab, CheckPolygonFilamentSoAMatchesScalar) {
	Vec vx, vy, vz;
	for (int i = 0; i < 40; ++i) {
		double t = 2.0 * M_PI * i / 39;
		vx.push_back((1.0 + 0.3 * cos(t)) * cos(0.2 * sin(t)));
		vy.push_back((1.0 + 0.3 * cos(t)) * sin(0.2 * sin(t)));
		vz.push_back(0.3 * sin(t) + 0.01 * t);
		if (i == 10) { // zero-length segment
			vx.push_back(vx.back()); vy.push_back(vy.back()); vz.push_back(vz.back());
		}
	}
	Vec ex, ey, ez;
	for (int i = 0; i < 1000; ++i) { // generic locations; not a multiple of the block size
		ex.push_back(1.0 + 1.5 * sin(1.3 * i)); ey.push_back(1.5 * sin(2.9 * i + 1)); ez.push_back(0.8 * sin(0.7 * i + 2));
	}
	const double along[] = {-3.0, -0.5, 0.0, 0.3, 0.5, 1.0, 1.2, 1.5, 2.0};
	const double perp[] = {0.0, 1e-14, 1e-9, 1e-3, 0.2, 0.5, 3.0};
	for (int seg : {0, 5, 11, 20, 37}) {
		double d[3] = {vx[seg + 1] - vx[seg], vy[seg + 1] - vy[seg], vz[seg + 1] - vz[seg]};
		double l = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
		double px = -d[1] / hypot(d[0], d[1]), py = d[0] / hypot(d[0], d[1]); // perpendicular to the segment
		for (double a : along) {
			for (double q : perp) {
				ex.push_back(vx[seg] + a * d[0] + q * l * px);
				ey.push_back(vy[seg] + a * d[1] + q * l * py);
				ez.push_back(vz[seg] + a * d[2]);
			}
		}
	}

	for (bool compensated : {true, false}) {
		SoAResult ref = runScalar(vx, vy, vz, ex, ey, ez, compensated);
		for (int numProcessors : {1, 3}) {
			SoAResult r = runSoA(vx, vy, vz, ex, ey, ez, numProcessors, compensated);
			EXPECT_LT(maxRelDiff(r.bx, r.by, r.bz, ref.bx, ref.by, ref.bz), 2.0e-14);
			EXPECT_LT(maxRelDiff(r.ax, r.ay, r.az, ref.ax, ref.ay, ref.az), 2.0e-14);
		}
	}

	// many segments, few locations: parallelized over segments
	Vec fx(ex.begin(), ex.begin() + 5), fy(ey.begin(), ey.begin() + 5), fz(ez.begin(), ez.begin() + 5);
	SoAResult serial = runSoA(vx, vy, vz, fx, fy, fz, 1), threaded = runSoA(vx, vy, vz, fx, fy, fz, 8);
	EXPECT_LT(maxRelDiff(threaded.bx, threaded.by, threaded.bz, serial.bx, serial.by, serial.bz), 1.0e-13);

	// the public kernels on a sub-range of locations agree with the full routine
	Vec kx(5, 0.0), ky(5, 0.0), kz(5, 0.0);
	kernelMagneticFieldPolygonFilament(vx.data(), vy.data(), vz.data(), 1.0, fx.data(), fy.data(), fz.data(),
			kx.data(), ky.data(), kz.data(), 0, (int) vx.size() - 1, 1, 3, true);
	EXPECT_EQ(kx[0], 0.0);
	EXPECT_EQ(kx[3], 0.0);
	for (int i = 1; i < 3; ++i) {
		EXPECT_EQ(kx[i], serial.bx[i]);
		EXPECT_EQ(kz[i], serial.bz[i]);
	}

	// the array-of-structs routines give the same result as the SoA ones, and add to the target
	Vec v(3 * vx.size()), e(3 * 5), b(3 * 5, 1.0);
	for (size_t i = 0; i < vx.size(); ++i) { v[3 * i] = vx[i]; v[3 * i + 1] = vy[i]; v[3 * i + 2] = vz[i]; }
	for (int i = 0; i < 5; ++i) { e[3 * i] = fx[i]; e[3 * i + 1] = fy[i]; e[3 * i + 2] = fz[i]; }
	magneticFieldPolygonFilament((int) vx.size(), v.data(), 1.0, 5, e.data(), b.data(), 1, true);
	for (int i = 0; i < 5; ++i) {
		EXPECT_EQ(b[3 * i], 1.0 + serial.bx[i]);
		EXPECT_EQ(b[3 * i + 2], 1.0 + serial.bz[i]);
	}
}

// squares of these normalized coordinates under- or overflow: the scalar routines must take over
TEST(TestAbscab, CheckPolygonFilamentSoAExtremeMagnitudes) {
	for (double length : {1.0e150, 1.0e-150}) {
		Vec ex, ey, ez;
		for (double offset : {1e-12, 1e-8, 1e-3, 1.0, 1e3, 1e100}) {
			for (double x : {offset, -offset, length + offset, 0.5 * length}) {
				ex.push_back(x); ey.push_back(offset); ez.push_back(0.5 * offset);
			}
		}
		Vec vx = {0.0, length}, vy = {0.0, 0.0}, vz = {0.0, 0.0};
		SoAResult ref = runScalar(vx, vy, vz, ex, ey, ez), r = runSoA(vx, vy, vz, ex, ey, ez);
		EXPECT_LT(maxRelDiff(r.bx, r.by, r.bz, ref.bx, ref.by, ref.bz), 1.0e-13) << "length " << length;
		EXPECT_LT(maxRelDiff(r.ax, r.ay, r.az, ref.ax, ref.ay, ref.az), 1.0e-13) << "length " << length;
	}
}

} // namespace abscab


