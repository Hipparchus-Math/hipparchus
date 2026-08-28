package org.hipparchus.optim.nonlinear.vector.constrained.sqp.electrolyte;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

/** ripopt electrolyte problem 3: NaCl strong-electrolyte speciation. */
public class ElectrolyteNaClSpeciationTest {

    /** Reference objective for this benchmark problem. */
    private static final double EXPECTED_OBJECTIVE = -0.48326822740568581;

    /** Absolute tolerance used for the objective comparison. */
    private static final double OBJECTIVE_TOLERANCE = 1.0000000000000001e-05;

    private static final int N = 4;
    private static final double[] CHARGES = { 1.0, -1.0, 1.0, -1.0 };
    private static final double[] DH_A = { 4.0, 3.0, 9.0, 3.5 };
    private static final double[] DH_B = { 0.075, 0.015, 0.0, 0.0 };
    private static final double[] MU0 = { 0.0, 0.0, 0.0, 14.0 * ElectrolyteTestSupport.LN10 };

    private static final class Equality extends EqualityConstraint {
        Equality() { super(new ArrayRealVector(3)); }
        @Override public int dim() { return N; }
        @Override public RealVector value(final RealVector x) {
            final double na = FastMath.exp(x.getEntry(0));
            final double cl = FastMath.exp(x.getEntry(1));
            final double h = FastMath.exp(x.getEntry(2));
            final double oh = FastMath.exp(x.getEntry(3));
            return new ArrayRealVector(new double[] { na - 0.1, cl - 0.1, na + h - cl - oh }, false);
        }
        @Override public RealMatrix jacobian(final RealVector x) {
            final RealMatrix j = MatrixUtils.createRealMatrix(3, N);
            final double na = FastMath.exp(x.getEntry(0));
            final double cl = FastMath.exp(x.getEntry(1));
            final double h = FastMath.exp(x.getEntry(2));
            final double oh = FastMath.exp(x.getEntry(3));
            j.setEntry(0, 0, na);
            j.setEntry(1, 1, cl);
            j.setEntry(2, 0, na);
            j.setEntry(2, 1, -cl);
            j.setEntry(2, 2, h);
            j.setEntry(2, 3, -oh);
            return j;
        }
    }

    @Test
    public void testNaClSpeciation() {
        ElectrolyteTestSupport.solve(
                "NaCl speciation",
                new ElectrolyteTestSupport.GibbsObjective(MU0, CHARGES, DH_A, DH_B),
                new Equality(),
                new double[] { FastMath.log(0.1), FastMath.log(0.1),
                               FastMath.log(1.0e-7), FastMath.log(1.0e-7) },
                new SimpleBounds(ElectrolyteTestSupport.filled(N, -30.0),
                                 ElectrolyteTestSupport.filled(N, 1.0)),
                EXPECTED_OBJECTIVE,
                OBJECTIVE_TOLERANCE,
                1.0e-4);
    }
}
