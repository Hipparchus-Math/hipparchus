package org.hipparchus.optim.nonlinear.vector.constrained.sqp.electrolyte;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

/** ripopt electrolyte problem 8: water-butanol-NaCl liquid-liquid equilibrium. */
public class ElectrolyteButanolWaterLleTest {

    /** Reference objective for this benchmark problem. */
    private static final double EXPECTED_OBJECTIVE = 7.8258405870940747e-10;

    /** Absolute tolerance used for the objective comparison. */
    private static final double OBJECTIVE_TOLERANCE = 1e-08;

    private static final double TAU12 = 0.50;
    private static final double TAU21 = 4.50;
    private static final double ALPHA = 0.40;
    private static final double KS = 0.19;
    private static final double M_NACL = 1.0;

    private static final class Objective extends ElectrolyteTestSupport.Objective {
        Objective() { super(2); }
        @Override public double value(final RealVector x) {
            final double d0 = x.getEntry(0) - 0.006;
            final double d1 = x.getEntry(1) - 0.48;
            return 1.0e-6 * (d0 * d0 + d1 * d1);
        }
        @Override public RealVector gradient(final RealVector x) {
            return new ArrayRealVector(new double[] {
                2.0e-6 * (x.getEntry(0) - 0.006),
                2.0e-6 * (x.getEntry(1) - 0.48)
            }, false);
        }
    }

    private static final class Equality extends EqualityConstraint {
        Equality() { super(new ArrayRealVector(2)); }
        @Override public int dim() { return 2; }
        @Override public RealVector value(final RealVector x) {
            final double xAq = x.getEntry(0);
            final double xOrg = x.getEntry(1);
            final double[] aq = ElectrolyteTestSupport.nrtlBinary(xAq, TAU12, TAU21, ALPHA);
            final double[] org = ElectrolyteTestSupport.nrtlBinary(xOrg, TAU12, TAU21, ALPHA);
            return new ArrayRealVector(new double[] {
                aq[0] + FastMath.log(xAq) + KS * M_NACL - org[0] - FastMath.log(xOrg),
                aq[1] + FastMath.log(1.0 - xAq) - org[1] - FastMath.log(1.0 - xOrg)
            }, false);
        }
        @Override public RealMatrix jacobian(final RealVector x) {
            final double xAq = x.getEntry(0);
            final double xOrg = x.getEntry(1);
            final RealMatrix j = MatrixUtils.createRealMatrix(2, 2);
            j.setEntry(0, 0, ElectrolyteTestSupport.dNrtlLnGamma1Dx1(xAq, TAU12, TAU21, ALPHA) + 1.0 / xAq);
            j.setEntry(0, 1, -ElectrolyteTestSupport.dNrtlLnGamma1Dx1(xOrg, TAU12, TAU21, ALPHA) - 1.0 / xOrg);
            j.setEntry(1, 0, ElectrolyteTestSupport.dNrtlLnGamma2Dx1(xAq, TAU12, TAU21, ALPHA) - 1.0 / (1.0 - xAq));
            j.setEntry(1, 1, -ElectrolyteTestSupport.dNrtlLnGamma2Dx1(xOrg, TAU12, TAU21, ALPHA) + 1.0 / (1.0 - xOrg));
            return j;
        }
    }

    @Test
    public void testButanolWaterLle() {
        ElectrolyteTestSupport.solve(
                "BuOH-water LLE",
                new Objective(),
                new Equality(),
                new double[] { 0.006, 0.48 },
                new SimpleBounds(new double[] { 1.0e-4, 0.3 }, new double[] { 0.05, 0.95 }),
                EXPECTED_OBJECTIVE,
                OBJECTIVE_TOLERANCE,
                1.0e-4);
    }
}
