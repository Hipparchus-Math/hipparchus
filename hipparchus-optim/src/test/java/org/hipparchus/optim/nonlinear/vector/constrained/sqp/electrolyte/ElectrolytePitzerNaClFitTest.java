package org.hipparchus.optim.nonlinear.vector.constrained.sqp.electrolyte;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.SimpleBounds;
import org.junit.jupiter.api.Test;

/** ripopt electrolyte problem 10: Pitzer NaCl parameter fit. */
public class ElectrolytePitzerNaClFitTest {

    /** Reference objective for this benchmark problem. */
    private static final double EXPECTED_OBJECTIVE = 0;

    /** Absolute tolerance used for the objective comparison. */
    private static final double OBJECTIVE_TOLERANCE = 9.9999999999999995e-07;

    private static final double[] TRUE_PARAMS = { 0.0765, 0.2664, 0.00127 };
    private static final double[] MOLALITIES = { 0.1, 0.2, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0 };

    private static double[] data() {
        final double[] data = new double[MOLALITIES.length];
        for (int i = 0; i < MOLALITIES.length; ++i) {
            data[i] = ElectrolyteTestSupport.pitzerOsmotic(
                    MOLALITIES[i], TRUE_PARAMS[0], TRUE_PARAMS[1], TRUE_PARAMS[2]);
        }
        return data;
    }

    private static final class Objective extends ElectrolyteTestSupport.Objective {
        Objective() { super(3); }
        @Override public double value(final RealVector x) {
            final double[] observations = data();
            double value = 0.0;
            for (int i = 0; i < MOLALITIES.length; ++i) {
                final double residual = ElectrolyteTestSupport.pitzerOsmotic(
                        MOLALITIES[i], x.getEntry(0), x.getEntry(1), x.getEntry(2)) - observations[i];
                value += residual * residual;
            }
            return value;
        }
        @Override public RealVector gradient(final RealVector x) {
            final double[] observations = data();
            final RealVector gradient = new ArrayRealVector(3);
            for (int i = 0; i < MOLALITIES.length; ++i) {
                final double m = MOLALITIES[i];
                final double residual = ElectrolyteTestSupport.pitzerOsmotic(
                        m, x.getEntry(0), x.getEntry(1), x.getEntry(2)) - observations[i];
                final double[] partials = ElectrolyteTestSupport.pitzerOsmoticPartials(m);
                for (int j = 0; j < 3; ++j) {
                    gradient.addToEntry(j, 2.0 * residual * partials[j]);
                }
            }
            return gradient;
        }
    }

    @Test
    public void testPitzerNaClFit() {
        ElectrolyteTestSupport.solve(
                "Pitzer NaCl fit",
                new Objective(),
                new double[] { 0.2, 0.5, 0.005 },
                new SimpleBounds(new double[] { -1.0, -1.0, -0.1 },
                                 new double[] { 1.0, 2.0, 0.1 }),
                EXPECTED_OBJECTIVE,
                OBJECTIVE_TOLERANCE);
    }
}
