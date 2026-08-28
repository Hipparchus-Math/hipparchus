package org.hipparchus.optim.nonlinear.vector.constrained.sqp.electrolyte;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

/** ripopt electrolyte problem 2: CO2-water speciation. */
public class ElectrolyteCo2WaterSpeciationTest {

    /** Reference objective for this benchmark problem. */
    private static final double EXPECTED_OBJECTIVE = -6.933670187752001e-3;

    /** Absolute tolerance used for the objective comparison. */
    private static final double OBJECTIVE_TOLERANCE = 1.0e-4;

    private static final int N = 5;
    private static final double C_TOTAL = 0.001;
    private static final double[] CHARGES = { 0.0, -1.0, -2.0, 1.0, -1.0 };
    private static final double[] DH_A = { 0.0, 4.0, 5.4, 9.0, 3.5 };
    private static final double[] DH_B = { 0.0, 0.0, 0.0, 0.0, 0.0 };
    private static final double[] MU0 = {
        0.0,
        6.35 * ElectrolyteTestSupport.LN10,
        (6.35 + 10.33) * ElectrolyteTestSupport.LN10,
        0.0,
        14.0 * ElectrolyteTestSupport.LN10
    };

    private static final class Equality extends EqualityConstraint {
        Equality() { super(new ArrayRealVector(2)); }
        @Override public int dim() { return N; }
        @Override public RealVector value(final RealVector x) {
            return new ArrayRealVector(new double[] {
                FastMath.exp(x.getEntry(0)) + FastMath.exp(x.getEntry(1)) +
                FastMath.exp(x.getEntry(2)) - C_TOTAL,
                -FastMath.exp(x.getEntry(1)) - 2.0 * FastMath.exp(x.getEntry(2)) +
                FastMath.exp(x.getEntry(3)) - FastMath.exp(x.getEntry(4))
            }, false);
        }
        @Override public RealMatrix jacobian(final RealVector x) {
            final RealMatrix j = MatrixUtils.createRealMatrix(2, N);
            j.setEntry(0, 0, FastMath.exp(x.getEntry(0)));
            j.setEntry(0, 1, FastMath.exp(x.getEntry(1)));
            j.setEntry(0, 2, FastMath.exp(x.getEntry(2)));
            j.setEntry(1, 1, -FastMath.exp(x.getEntry(1)));
            j.setEntry(1, 2, -2.0 * FastMath.exp(x.getEntry(2)));
            j.setEntry(1, 3, FastMath.exp(x.getEntry(3)));
            j.setEntry(1, 4, -FastMath.exp(x.getEntry(4)));
            return j;
        }
    }

    @Test
    public void testCo2WaterSpeciation() {
        final double base = FastMath.log(C_TOTAL / 3.0);
        ElectrolyteTestSupport.solve(
                "CO2-water speciation",
                new ElectrolyteTestSupport.GibbsObjective(MU0, CHARGES, DH_A, DH_B),
                new Equality(),
                new double[] { base, base, FastMath.log(1.0e-8),
                               FastMath.log(1.0e-6), FastMath.log(1.0e-8) },
                new SimpleBounds(ElectrolyteTestSupport.filled(N, -35.0),
                                 ElectrolyteTestSupport.filled(N, -2.0)),
                EXPECTED_OBJECTIVE,
                OBJECTIVE_TOLERANCE,
                1.0e-4);
    }
}
