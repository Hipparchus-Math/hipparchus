package org.hipparchus.optim.nonlinear.vector.constrained.sqp.electrolyte;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

/** ripopt electrolyte problem 11: simultaneous Debye-Huckel fit for three salts. */
public class ElectrolyteMultiSaltDhFitTest {

    /** Reference objective for this benchmark problem. */
    private static final double EXPECTED_OBJECTIVE = 0;

    /** Absolute tolerance used for the objective comparison. */
    private static final double OBJECTIVE_TOLERANCE = 0.0001;

    private static final int N = 8;
    private static final double[] TRUE_PARAMS = { 4.0, 0.075, 3.0, 0.015, 6.0, 0.165, 3.0, 0.015 };
    private static final double[] MOLALITIES = { 0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0 };

    private static double nacl(final RealVector x, final double m) {
        return 0.5 * (ElectrolyteTestSupport.lnGammaDh(1.0, x.getEntry(0), x.getEntry(1), m) +
                      ElectrolyteTestSupport.lnGammaDh(1.0, x.getEntry(6), x.getEntry(7), m));
    }

    private static double kcl(final RealVector x, final double m) {
        return 0.5 * (ElectrolyteTestSupport.lnGammaDh(1.0, x.getEntry(2), x.getEntry(3), m) +
                      ElectrolyteTestSupport.lnGammaDh(1.0, x.getEntry(6), x.getEntry(7), m));
    }

    private static double cacl2(final RealVector x, final double m) {
        final double ionicStrength = 3.0 * m;
        return (ElectrolyteTestSupport.lnGammaDh(2.0, x.getEntry(4), x.getEntry(5), ionicStrength) +
                2.0 * ElectrolyteTestSupport.lnGammaDh(1.0, x.getEntry(6), x.getEntry(7), ionicStrength)) / 3.0;
    }

    private static RealVector trueParameters() {
        return new ArrayRealVector(TRUE_PARAMS, false);
    }

    private static final class Objective extends ElectrolyteTestSupport.Objective {
        Objective() { super(N); }
        @Override public double value(final RealVector x) {
            final RealVector truth = trueParameters();
            double value = 0.0;
            for (final double m : MOLALITIES) {
                double residual = nacl(x, m) - nacl(truth, m);
                value += residual * residual;
            }
            for (final double m : MOLALITIES) {
                double residual = kcl(x, m) - kcl(truth, m);
                value += residual * residual;
            }
            for (final double m : MOLALITIES) {
                double residual = cacl2(x, m) - cacl2(truth, m);
                value += residual * residual;
            }
            return value;
        }
        @Override public RealVector gradient(final RealVector x) {
            // Source-faithful forward-difference gradient.
            final RealVector gradient = new ArrayRealVector(N);
            final double f0 = value(x);
            for (int i = 0; i < N; ++i) {
                final double h = 1.0e-7 * FastMath.max(FastMath.abs(x.getEntry(i)), 1.0e-5);
                final RealVector xp = x.copy();
                xp.addToEntry(i, h);
                gradient.setEntry(i, (value(xp) - f0) / h);
            }
            return gradient;
        }
    }

    @Test
    public void testMultiSaltDhFit() {
        final double[] initial = new double[N];
        final double[] lower = new double[N];
        final double[] upper = new double[N];
        for (int i = 0; i < 4; ++i) {
            initial[2 * i] = 3.0;
            initial[2 * i + 1] = 0.0;
            lower[2 * i] = 1.0;
            upper[2 * i] = 10.0;
            lower[2 * i + 1] = -0.5;
            upper[2 * i + 1] = 0.5;
        }
        ElectrolyteTestSupport.solve(
                "Multi-salt DH fit",
                new Objective(),
                initial,
                new SimpleBounds(lower, upper),
                EXPECTED_OBJECTIVE,
                OBJECTIVE_TOLERANCE);
    }
}
