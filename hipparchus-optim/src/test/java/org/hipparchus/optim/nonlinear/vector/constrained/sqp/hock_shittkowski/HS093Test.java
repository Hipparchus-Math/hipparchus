package org.hipparchus.optim.nonlinear.vector.constrained.sqp.hock_shittkowski;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * HS093 (TP93) — 6 vars, 2 nonlinear inequalities (Fortran G(i) <= 0).
 * Bounds handled separately: x_i >= 0.
 */
public class HS093Test {

    private static final double TOL = 1e-6;

    /** f(x) from TP93. */
    private static final class TP93Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 6; }

        @Override public double value(RealVector x) {
            double x1 = x.getEntry(0), x2 = x.getEntry(1), x3 = x.getEntry(2);
            double x4 = x.getEntry(3), x5 = x.getEntry(4), x6 = x.getEntry(5);

            double v1 = x1 + x2 + x3;
            double v2 = x1 + 1.57 * x2 + x4;
            double v3 = x1 * x4;
            double v4 = x3 * x2;

            return 0.0204 * v3 * v1
                 + 0.0187 * v4 * v2
                 + 0.0607 * v3 * v1 * x5 * x5
                 + 0.0437 * v4 * v2 * x6 * x6;
        }

        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x)  { throw new UnsupportedOperationException(); }
    }

    /**
     * Two nonlinear inequalities, using Fortran's sign convention: G(x) <= 0.
     * G1 = 1e-3*x1*x2*x3*x4*x5*x6 - 2.07
     * G2 = 1 - 6.2e-4*x1*x4*x5^2*(x1+x2+x3) - 5.8e-4*x2*x3*x6^2*(x1+1.57*x2+x4)
     */
    private static final class TP93Ineq extends InequalityConstraint {
        TP93Ineq() {
            // vector sized for 2 inequalities
            super(new ArrayRealVector(new double[]{0.0, 0.0}));
        }
        @Override public RealVector value(RealVector x) {
            double x1 = x.getEntry(0), x2 = x.getEntry(1), x3 = x.getEntry(2);
            double x4 = x.getEntry(3), x5 = x.getEntry(4), x6 = x.getEntry(5);

            double g1 = 1.0e-3 * x1 * x2 * x3 * x4 * x5 * x6 - 2.07;

            double v3 = (x1 + x2 + x3);
            double v4 = (x1 + 1.57 * x2 + x4);
            double g2 = 1.0
                    - 6.2e-4 * x1 * x4 * x5 * x5 * v3
                    - 5.8e-4 * x2 * x3 * x6 * x6 * v4;

            return new ArrayRealVector(new double[]{ g1, g2 });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public int dim() { return 6; }
    }

    @Test
    public void testTP93() {
        // Fortran initial guess
        double[] x0 = { 5.54, 4.4, 12.02, 11.82, 0.702, 0.852 };

        // Lower bounds 0, no upper bounds
        double[] lo = {0,0,0,0,0,0};
        double[] hi = {
                Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY,
                Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY
        };

        // Reported solution (XEX) and FEX from the Fortran block
        double[] xex = {
                0.533266639884e+01,
                0.465674439073e+01,
                0.104329901123e+02,
                0.120823085893e+02,
                0.752607369745e+00,
                0.878650836850e+00
        };
        double fex = 0.135075961229e+03; // 135.075961229

        // Print objective at reported XEX to verify
        double fAtXex = new TP93Obj().value(new ArrayRealVector(xex));
        System.out.printf("TP93 objective at XEX = %.12f (expected %.12f)%n", fAtXex, fex);
        assertEquals(fex, fAtXex, 1e-5);

        SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = optimizer.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new TP93Obj()),
                new TP93Ineq(),
                new SimpleBounds(lo, hi)
        );

        System.out.printf("TP93 optimizer f* = %.12f%n", sol.getValue());
        // Tightish tolerance (these problems are well-scaled)
        HSProblemTestUtils.assertExpectedObjective(fex, sol);
    }
}

