/*
 * Licensed to the Hipparchus project under one or more contributor license agreements...
 */
package org.hipparchus.optim.nonlinear.vector.constrained.sqp.hock_shittkowski;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

/** HS TP114. 10 vars, 8 inequalities (G1..G8), 3 equalities (G9..G11). */
public class HS114Test {

    // Bounds dal MODE=1
    private static final double[] LB = {
        1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 85.0, 90.0, 3.0, 1.2, 145.0
    };
    private static final double[] UB = {
        2000.0, 16000.0, 120.0, 5000.0, 2000.0, 93.0, 95.0, 12.0, 4.0, 162.0
    };

    /** f(x) = 5.04*x1 + 0.035*x2 + 10*x3 + 3.36*x5 - 0.063*x4*x7. */
    private static class TP114Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 10; }
        @Override public double value(RealVector X) {
            final double x1=X.getEntry(0), x2=X.getEntry(1), x3=X.getEntry(2),
                         x4=X.getEntry(3), x5=X.getEntry(4), x7=X.getEntry(6);
            return 5.04*x1 + 0.035*x2 + 10.0*x3 + 3.36*x5 - 0.063*x4*x7;
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x)  { throw new UnsupportedOperationException(); }
    }

    /** 8 disuguaglianze (G1..G8) in forma g(x) >= 0. */
    private static class TP114Ineq extends InequalityConstraint {
        TP114Ineq() { super(new ArrayRealVector(new double[]{0,0,0,0,0,0,0,0})); }
        @Override public int dim() { return 10; }
        @Override public RealVector value(RealVector X) {
            final double x1=X.getEntry(0), x2=X.getEntry(1), x3=X.getEntry(2),
                         x4=X.getEntry(3), x6=X.getEntry(5), x7=X.getEntry(6),
                         x8=X.getEntry(7), x9=X.getEntry(8), x10=X.getEntry(9);

            final double g1 = 35.82 - 0.222*x10 - 0.9*x9;
            final double g2 = -133.0 + 3.0*x7 - 0.99*x10;
            final double g3 = -35.82 + 0.222*x10 + (10.0/9.0)*x9;
            final double g4 = 133.0 - 3.0*x7 + x10/0.99;

            final double g5 = 1.12*x1 + 0.13167*x1*x8 - 6.67e-3*x1*x8*x8 - 0.99*x4;
            final double g6 = 57.425 + 1.098*x8 - 0.038*x8*x8 + 0.325*x6 - 0.99*x7;
            final double g7 = -1.12*x1 - 0.13167*x1*x8 + 6.67e-3*x1*x8*x8 + x4/0.99;
            final double g8 = -57.425 - 1.098*x8 + 0.038*x8*x8 - 0.325*x6 + x7/0.99;

            return new ArrayRealVector(new double[]{ g1,g2,g3,g4,g5,g6,g7,g8 });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    /** 3 uguaglianze (G9..G11) in forma h(x) = 0. */
    private static class TP114Eq extends EqualityConstraint {
        TP114Eq() { super(new ArrayRealVector(new double[]{0.0, 0.0, 0.0})); }
        @Override public int dim() { return 10; }
        @Override public RealVector value(RealVector X) {
            final double x1=X.getEntry(0), x2=X.getEntry(1), x3=X.getEntry(2),
                         x4=X.getEntry(3), x5=X.getEntry(4), x6=X.getEntry(5),
                         x8=X.getEntry(7), x9=X.getEntry(8), x10=X.getEntry(9);

            // G9: 1.22*x4 - x1 - x5 = 0
            final double h1 = 1.22*x4 - x1 - x5;

            // G10: 9.8e4*x3/(x4*x9 + 1e3*x3) - x6 = 0
            final double denom = x4*x9 + 1.0e3*x3;
            final double h2 = 9.8e4*x3/denom - x6;

            // G11: (x2 + x5)/x1 - x8 = 0
            final double h3 = (x2 + x5)/x1 - x8;

            return new ArrayRealVector(new double[]{ h1, h2, h3 });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    @Test
    public void testHS114() {
        final InitialGuess guess = new InitialGuess(new double[]{
            1.745e3, 1.2e4, 110.0, 3.048e3, 1.974e3, 89.2, 92.8, 8.0, 3.6, 145.0
        });

        final SimpleBounds bounds = new SimpleBounds(LB, UB);

        final SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        final LagrangeSolution sol = opt.optimize(
            guess,
            new ObjectiveFunction(new TP114Obj()),
            new TP114Eq(),    // 3 equalities
            new TP114Ineq(),  // 8 inequalities
            bounds
        );

        // FEX = -0.176880696344D+04
        HSProblemTestUtils.assertExpectedObjective(-1768.80696344, sol);
    }
}

