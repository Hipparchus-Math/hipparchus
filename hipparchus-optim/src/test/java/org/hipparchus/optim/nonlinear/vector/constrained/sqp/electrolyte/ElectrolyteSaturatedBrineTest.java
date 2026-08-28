package org.hipparchus.optim.nonlinear.vector.constrained.sqp.electrolyte;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

/** ripopt electrolyte problem 9: saturated-brine VLE and SLE. */
public class ElectrolyteSaturatedBrineTest {

    /** Reference objective for this benchmark problem. */
    private static final double EXPECTED_OBJECTIVE = 0;

    /** Absolute tolerance used for the objective comparison. */
    private static final double OBJECTIVE_TOLERANCE = 0;

    private static final double BETA0 = 0.0765;
    private static final double BETA1 = 0.2664;
    private static final double C_PHI = 0.00127;
    private static final double LN_KSP = 3.627;
    private static final double P_W_PURE = 3.169;

    private static final class Objective extends ElectrolyteTestSupport.Objective {
        Objective() { super(3); }
        @Override public double value(final RealVector x) { return 0.0; }
        @Override public RealVector gradient(final RealVector x) { return new ArrayRealVector(3); }
    }

    private static final class Equality extends EqualityConstraint {
        Equality() { super(new ArrayRealVector(3)); }
        @Override public int dim() { return 3; }
        @Override public RealVector value(final RealVector x) {
            final double m = x.getEntry(0);
            final double aw = x.getEntry(1);
            final double pw = x.getEntry(2);
            final double phi = ElectrolyteTestSupport.pitzerOsmotic(m, BETA0, BETA1, C_PHI);
            return new ArrayRealVector(new double[] {
                2.0 * ElectrolyteTestSupport.pitzerLnGammaPm(m, BETA0, BETA1, C_PHI) +
                2.0 * FastMath.log(m) - LN_KSP,
                aw - FastMath.exp(-phi * 2.0 * m * ElectrolyteTestSupport.M_W),
                pw - aw * P_W_PURE
            }, false);
        }
        @Override public RealMatrix jacobian(final RealVector x) {
            final double m = x.getEntry(0);
            final double phi = ElectrolyteTestSupport.pitzerOsmotic(m, BETA0, BETA1, C_PHI);
            final double dPhi = ElectrolyteTestSupport.dPitzerOsmoticDm(m, BETA0, BETA1, C_PHI);
            final double exp = FastMath.exp(-phi * 2.0 * m * ElectrolyteTestSupport.M_W);
            final RealMatrix j = MatrixUtils.createRealMatrix(3, 3);
            j.setEntry(0, 0, 2.0 * ElectrolyteTestSupport.dPitzerLnGammaPmDm(m, BETA0, BETA1, C_PHI) + 2.0 / m);
            j.setEntry(1, 0, exp * (dPhi * 2.0 * m * ElectrolyteTestSupport.M_W +
                                    phi * 2.0 * ElectrolyteTestSupport.M_W));
            j.setEntry(1, 1, 1.0);
            j.setEntry(2, 1, -P_W_PURE);
            j.setEntry(2, 2, 1.0);
            return j;
        }
    }

    @Test
    public void testSaturatedBrine() {
        ElectrolyteTestSupport.solve(
                "Saturated brine",
                new Objective(),
                new Equality(),
                new double[] { 3.0, 0.8, 2.5 },
                new SimpleBounds(new double[] { 0.1, 0.5, 1.0 },
                                 new double[] { 15.0, 1.0, 3.5 }),
                EXPECTED_OBJECTIVE,
                OBJECTIVE_TOLERANCE,
                1.0e-4);
    }
}
