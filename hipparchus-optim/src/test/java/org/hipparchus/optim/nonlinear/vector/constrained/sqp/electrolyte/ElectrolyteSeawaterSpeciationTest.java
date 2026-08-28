package org.hipparchus.optim.nonlinear.vector.constrained.sqp.electrolyte;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

/** ripopt electrolyte problem 13: seawater speciation. */
public class ElectrolyteSeawaterSpeciationTest {

    /** Reference objective for this benchmark problem. */
    private static final double EXPECTED_OBJECTIVE = -1.3482721229153687;

    /** Absolute tolerance used for the objective comparison. */
    private static final double OBJECTIVE_TOLERANCE = 0.0001;

    private static final int N = 15;
    private static final double[] CHARGES = {
        1.0, 1.0, 2.0, 2.0, -1.0, -2.0, -1.0, -2.0,
        1.0, -1.0, 0.0, 0.0, 1.0, -1.0, -1.0
    };
    private static final double[] DH_A = {
        4.0, 3.0, 6.0, 6.0, 3.0, 5.0, 4.0, 5.4,
        9.0, 3.5, 0.0, 0.0, 4.0, 4.0, 3.5
    };
    private static final double[] DH_B = {
        0.075, 0.015, 0.165, 0.165, 0.015, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    };

    private static final double NA_TOTAL = 0.4861;
    private static final double K_TOTAL = 0.01058;
    private static final double MG_TOTAL = 0.05474;
    private static final double CA_TOTAL = 0.01065;
    private static final double CL_TOTAL = 0.5658;
    private static final double S_TOTAL = 0.02927;
    private static final double C_TOTAL = 0.002048;

    private static double[] mu0() {
        final double[] mu = new double[N];
        mu[7] = 10.33 * ElectrolyteTestSupport.LN10;
        mu[9] = 14.0 * ElectrolyteTestSupport.LN10;
        mu[10] = -2.23 * ElectrolyteTestSupport.LN10;
        mu[11] = -2.30 * ElectrolyteTestSupport.LN10;
        mu[12] = (14.0 - 2.58) * ElectrolyteTestSupport.LN10;
        mu[13] = -0.70 * ElectrolyteTestSupport.LN10;
        mu[14] = -0.85 * ElectrolyteTestSupport.LN10;
        return mu;
    }

    private static final class Objective extends ElectrolyteTestSupport.Objective {
        Objective() { super(N); }

        @Override
        public double value(final RealVector x) {
            final double[] mu = mu0();
            final double ionicStrength = ElectrolyteTestSupport.ionicStrength(x, CHARGES);
            double value = 0.0;
            for (int i = 0; i < N; ++i) {
                final double molality = FastMath.exp(x.getEntry(i));
                final double lnGamma = CHARGES[i] == 0.0 ?
                        0.1 * ionicStrength :
                        ElectrolyteTestSupport.lnGammaDh(CHARGES[i], DH_A[i], DH_B[i], ionicStrength);
                value += molality * (mu[i] + lnGamma + x.getEntry(i));
            }
            return value;
        }

        @Override
        public RealVector gradient(final RealVector x) {
            final double[] mu = mu0();
            final double ionicStrength = ElectrolyteTestSupport.ionicStrength(x, CHARGES);
            final double[] lnGamma = new double[N];
            final double[] dLnGamma = new double[N];
            double s1 = 0.0;
            for (int i = 0; i < N; ++i) {
                if (CHARGES[i] == 0.0) {
                    lnGamma[i] = 0.1 * ionicStrength;
                    dLnGamma[i] = 0.1;
                } else {
                    lnGamma[i] = ElectrolyteTestSupport.lnGammaDh(CHARGES[i], DH_A[i], DH_B[i], ionicStrength);
                    dLnGamma[i] = ElectrolyteTestSupport.dLnGammaDhDi(CHARGES[i], DH_A[i], DH_B[i], ionicStrength);
                }
                s1 += FastMath.exp(x.getEntry(i)) * dLnGamma[i];
            }

            final RealVector gradient = new ArrayRealVector(N);
            for (int j = 0; j < N; ++j) {
                final double molality = FastMath.exp(x.getEntry(j));
                final double a = mu[j] + lnGamma[j] + x.getEntry(j) + 1.0;
                final double dIdx = 0.5 * CHARGES[j] * CHARGES[j] * molality;
                gradient.setEntry(j, molality * a + dIdx * s1);
            }
            return gradient;
        }
    }

    private static final class Equality extends EqualityConstraint {
        Equality() { super(new ArrayRealVector(8)); }
        @Override public int dim() { return N; }

        @Override
        public RealVector value(final RealVector x) {
            final double[] m = new double[N];
            for (int i = 0; i < N; ++i) { m[i] = FastMath.exp(x.getEntry(i)); }
            double charge = 0.0;
            for (int i = 0; i < N; ++i) { charge += CHARGES[i] * m[i]; }
            return new ArrayRealVector(new double[] {
                m[0] + m[13] - NA_TOTAL,
                m[1] + m[14] - K_TOTAL,
                m[2] + m[10] + m[12] - MG_TOTAL,
                m[3] + m[11] - CA_TOTAL,
                m[4] - CL_TOTAL,
                m[5] + m[10] + m[11] + m[13] + m[14] - S_TOTAL,
                m[6] + m[7] - C_TOTAL,
                charge
            }, false);
        }

        @Override
        public RealMatrix jacobian(final RealVector x) {
            final double[] m = new double[N];
            for (int i = 0; i < N; ++i) { m[i] = FastMath.exp(x.getEntry(i)); }
            final RealMatrix j = MatrixUtils.createRealMatrix(8, N);
            j.setEntry(0, 0, m[0]); j.setEntry(0, 13, m[13]);
            j.setEntry(1, 1, m[1]); j.setEntry(1, 14, m[14]);
            j.setEntry(2, 2, m[2]); j.setEntry(2, 10, m[10]); j.setEntry(2, 12, m[12]);
            j.setEntry(3, 3, m[3]); j.setEntry(3, 11, m[11]);
            j.setEntry(4, 4, m[4]);
            j.setEntry(5, 5, m[5]); j.setEntry(5, 10, m[10]); j.setEntry(5, 11, m[11]);
            j.setEntry(5, 13, m[13]); j.setEntry(5, 14, m[14]);
            j.setEntry(6, 6, m[6]); j.setEntry(6, 7, m[7]);
            for (int i = 0; i < N; ++i) { j.setEntry(7, i, CHARGES[i] * m[i]); }
            return j;
        }
    }

    @Test
    public void testSeawaterSpeciation() {
        ElectrolyteTestSupport.solve(
                "Seawater speciation",
                new Objective(),
                new Equality(),
                new double[] {
                    FastMath.log(0.48), FastMath.log(0.01), FastMath.log(0.05), FastMath.log(0.01),
                    FastMath.log(0.56), FastMath.log(0.025), FastMath.log(0.002), FastMath.log(1.0e-5),
                    FastMath.log(1.0e-8), FastMath.log(1.0e-6), FastMath.log(1.0e-3), FastMath.log(1.0e-3),
                    FastMath.log(1.0e-5), FastMath.log(1.0e-3), FastMath.log(1.0e-4)
                },
                new SimpleBounds(ElectrolyteTestSupport.filled(N, -46.0),
                                 ElectrolyteTestSupport.filled(N, 1.0)),
                EXPECTED_OBJECTIVE,
                OBJECTIVE_TOLERANCE,
                1.0e-3);
    }
}
