/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements...
 */
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
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

/** HS TP83. 5 vars, 6 nonlinear inequalities (box constraints on V1,V2,V3). */
public class HS083Test {

    // ---- Costanti dal COMMON/DATA83 ----
    private static final double A  = 5.3578547;
    private static final double B  = 0.8356891;
    private static final double Cc = 37.293239;
    private static final double D  = 4.0792141e4;

    private static final double A1  = 85.334407;
    private static final double A2  = 5.6858e-3;
    private static final double A3  = 6.262e-4;
    private static final double A4  = 2.2053e-3;

    private static final double A5  = 80.51249;
    private static final double A6  = 7.1317e-3;
    private static final double A7  = 2.9955e-3;
    private static final double A8  = 2.1813e-3;

    private static final double A9  = 9.300961;
    private static final double A10 = 4.7026e-3;
    private static final double A11 = 1.2547e-3;
    private static final double A12 = 1.9085e-3;

    // Bounds: x1∈[78,102], x2∈[33,45], x3..x5∈[27,45]
    private static final double[] LB = { 78.0, 33.0, 27.0, 27.0, 27.0 };
    private static final double[] UB = { 102.0, 45.0, 45.0, 45.0, 45.0 };

    /** f(x) = A*x3^2 + B*x1*x5 + C*x1 - D. */
    private static class TP83Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 5; }
        @Override public double value(RealVector X) {
            final double x1 = X.getEntry(0);
            final double x3 = X.getEntry(2);
            final double x5 = X.getEntry(4);
            return A*x3*x3 + B*x1*x5 + Cc*x1 - D;
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x)  { throw new UnsupportedOperationException(); }
    }

    /**
     * 6 disuguaglianze in forma g(x) >= 0:
     * { V1, V2, V3, 92 - V1, 20 - V2, 5 - V3 } >= 0
     */
    private static class TP83Ineq extends InequalityConstraint {
        TP83Ineq() { super(new ArrayRealVector(new double[]{ 0,0,0,0,0,0 })); }
        @Override public int dim() { return 5; }
        @Override public RealVector value(RealVector X) {
            final double x1 = X.getEntry(0), x2 = X.getEntry(1),
                         x3 = X.getEntry(2), x4 = X.getEntry(3), x5 = X.getEntry(4);

            final double V1 = A1 + A2*x2*x5 + A3*x1*x4 - A4*x3*x5;
            final double V2 = A5 + A6*x2*x5 + A7*x1*x2 + A8*x3*x3 - 90.0;
            final double V3 = A9 + A10*x3*x5 + A11*x1*x3 + A12*x3*x4 - 20.0;

            final double g1 = V1;
            final double g2 = V2;
            final double g3 = V3;
            final double g4 = 92.0 - V1;
            final double g5 = 20.0 - V2;
            final double g6 = 5.0  - V3;

            return new ArrayRealVector(new double[]{ g1, g2, g3, g4, g5, g6 });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    @Test
    public void testHS083() {
        final InitialGuess guess = new InitialGuess(new double[]{ 78.0, 33.0, 27.0, 27.0, 27.0 });
        final SimpleBounds bounds = new SimpleBounds(LB, UB);

        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();

        final LagrangeSolution sol = optimizer.optimize(
                guess,
                new ObjectiveFunction(new TP83Obj()),
                new TP83Ineq(),
                bounds
        );

        // FEX (TP83): -0.306655386717D+05
        HSProblemTestUtils.assertExpectedObjective(-30665.5386717, sol);
    }
}
