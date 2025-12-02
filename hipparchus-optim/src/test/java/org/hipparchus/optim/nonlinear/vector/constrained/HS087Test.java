package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class HS087Test {

    /** Objective piecewise: F1(x1) + F2(x2) come in TP87 MODE=2. */
    private static class HS087Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 6; }

        @Override
        public double value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);

            // F1
            final double F1 = (x1 < 300.0) ? 30.0 * x1 : 31.0 * x1;

            // F2
            final double F2;
            if (x2 < 100.0) {
                F2 = 28.0 * x2;
            } else if (x2 < 200.0) {
                F2 = 29.0 * x2;
            } else {
                F2 = 30.0 * x2;
            }

            return F1 + F2;
        }

        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    /**
     * 
     *  g1 = -x1 + 300 - (x3*x4/A)*cos(B - x6) + (C*x3^2/A)*D
     *  g2 = -x2 - (x3*x4/A)*cos(B + x6) + (C*x4^2/A)*D
     *  g3 = -x5 - (x3*x4/A)*sin(B + x6) + (C*x4^2/A)*E
     *  g4 = 200 - (x3*x4/A)*sin(B - x6) + (C*x3^2/A)*E
     */
    private static class HS087Eq extends EqualityConstraint {
        // Costanti dal MODE=0 (inizio subroutine)
        private static final double A = 131.078;
        private static final double B = 1.48477;
        private static final double D = FastMath.cos(1.47588);
        private static final double E = FastMath.sin(1.47588);
        private static final double Cc = 0.90798; // rinomino C→Cc per evitare clash col nome classe

        HS087Eq() { super(new ArrayRealVector(new double[]{0, 0, 0, 0})); }

        @Override public RealVector value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);
            final double x3 = x.getEntry(2);
            final double x4 = x.getEntry(3);
            final double x5 = x.getEntry(4);
            final double x6 = x.getEntry(5);

            final double x3x4_over_A = (x3 * x4) / A;
            final double term_x3_sq = (Cc * x3 * x3) / A;
            final double term_x4_sq = (Cc * x4 * x4) / A;

            final double g1 = -x1 + 300.0 - x3x4_over_A * FastMath.cos(B - x6) + term_x3_sq * D;
            final double g2 = -x2 - x3x4_over_A * FastMath.cos(B + x6) + term_x4_sq * D;
            final double g3 = -x5 - x3x4_over_A * FastMath.sin(B + x6) + term_x4_sq * E;
            final double g4 = 200.0 - x3x4_over_A * FastMath.sin(B - x6) + term_x3_sq * E;

            return new ArrayRealVector(new double[]{ g1, g2, g3, g4 });
        }

        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 6; }
    }

    @Test
    public void testHS087() {
        // Guess MODE=1
        InitialGuess guess = new InitialGuess(new double[]{
                390.0, 1000.0, 419.5, 340.5, 198.175, 0.5
        });

        // Bounds MODE=1
        SimpleBounds bounds = new SimpleBounds(
                new double[]{ 0.0, 0.0, 340.0, 340.0, -1000.0, 0.0 },
                new double[]{ 400.0, 1000.0, 420.0, 420.0, 1000.0, 0.5236 }
        );

        SQPOptimizerS2 optimizer = new SQPOptimizerS2();
        optimizer.setDebugPrinter(System.out::println); // richiesto

        double expected = 0.892759773493e4; // FEX

        LagrangeSolution sol = optimizer.optimize(
                guess,
                new ObjectiveFunction(new HS087Obj()),
                new HS087Eq(),
                bounds
        );

        assertEquals(expected, sol.getValue(), 1e-3);
    }
}
