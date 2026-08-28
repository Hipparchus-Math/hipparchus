package org.hipparchus.optim.nonlinear.vector.constrained.sqp.electrolyte;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

/** ripopt electrolyte problem 7: NaCl solubility. */
public class ElectrolyteNaClSolubilityTest {

    /** Reference objective for this benchmark problem. */
    private static final double EXPECTED_OBJECTIVE = 0;

    /** Absolute tolerance used for the objective comparison. */
    private static final double OBJECTIVE_TOLERANCE = 0.0001;

    private static final double BETA0 = 0.0765;
    private static final double BETA1 = 0.2664;
    private static final double C_PHI = 0.00127;
    private static final double LN_KSP = 3.627;

    private static final class Objective extends ElectrolyteTestSupport.Objective {
        Objective() { super(1); }
        @Override public double value(final RealVector x) {
            final double m = x.getEntry(0);
            final double residual = 2.0 * ElectrolyteTestSupport.pitzerLnGammaPm(m, BETA0, BETA1, C_PHI) +
                                    2.0 * FastMath.log(m) - LN_KSP;
            return residual * residual;
        }
        @Override public RealVector gradient(final RealVector x) {
            final double m = x.getEntry(0);
            final double residual = 2.0 * ElectrolyteTestSupport.pitzerLnGammaPm(m, BETA0, BETA1, C_PHI) +
                                    2.0 * FastMath.log(m) - LN_KSP;
            final double derivative = 2.0 * ElectrolyteTestSupport.dPitzerLnGammaPmDm(m, BETA0, BETA1, C_PHI) +
                                      2.0 / m;
            return new ArrayRealVector(new double[] { 2.0 * residual * derivative }, false);
        }
    }

    @Test
    public void testNaClSolubility() {
        ElectrolyteTestSupport.solve(
                "NaCl solubility",
                new Objective(),
                new double[] { 3.0 },
                new SimpleBounds(new double[] { 0.1 }, new double[] { 15.0 }),
                EXPECTED_OBJECTIVE,
                OBJECTIVE_TOLERANCE);
    }
}
