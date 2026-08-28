package org.hipparchus.optim.nonlinear.vector.constrained.sqp.electrolyte;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

/** ripopt electrolyte problem 1: water autoionization. */
public class ElectrolyteWaterAutoionizationTest {

    /** Reference objective for this benchmark problem. */
    private static final double EXPECTED_OBJECTIVE = 3.799326928221924e-7;

    /** Absolute tolerance used for the objective comparison. */
    private static final double OBJECTIVE_TOLERANCE = 1.0e-4;

    private static final class Objective extends ElectrolyteTestSupport.Objective {
        Objective() { super(1); }

        @Override
        public double value(final RealVector x) {
            final double m = FastMath.exp(x.getEntry(0));
            final double lgH = ElectrolyteTestSupport.lnGammaDh(1.0, 9.0, 0.0, m);
            final double lgOh = ElectrolyteTestSupport.lnGammaDh(1.0, 3.5, 0.0, m);
            return m * (lgH + x.getEntry(0)) +
                   m * (-ElectrolyteTestSupport.LN_KW + lgOh + x.getEntry(0));
        }

        @Override
        public RealVector gradient(final RealVector x) {
            final double xv = x.getEntry(0);
            final double m = FastMath.exp(xv);
            final double lgH = ElectrolyteTestSupport.lnGammaDh(1.0, 9.0, 0.0, m);
            final double lgOh = ElectrolyteTestSupport.lnGammaDh(1.0, 3.5, 0.0, m);
            final double dlgH = ElectrolyteTestSupport.dLnGammaDhDi(1.0, 9.0, 0.0, m);
            final double dlgOh = ElectrolyteTestSupport.dLnGammaDhDi(1.0, 3.5, 0.0, m);
            final double dfDm = lgH + lgOh - ElectrolyteTestSupport.LN_KW +
                                2.0 * xv + 2.0 + m * (dlgH + dlgOh);
            return new ArrayRealVector(new double[] { m * dfDm }, false);
        }
    }

    @Test
    public void testWaterAutoionization() {
        ElectrolyteTestSupport.solve(
                "Water autoionization",
                new Objective(),
                new double[] { FastMath.log(1.0e-5) },
                new SimpleBounds(new double[] { -23.0 }, new double[] { -6.9 }),
                EXPECTED_OBJECTIVE,
                OBJECTIVE_TOLERANCE);
    }
}
