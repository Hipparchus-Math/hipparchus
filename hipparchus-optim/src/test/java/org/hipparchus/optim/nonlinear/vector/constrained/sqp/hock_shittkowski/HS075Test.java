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
import org.hipparchus.optim.nonlinear.vector.constrained.EqualityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.InequalityConstraint;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

/** HS TP75 (same model as TP74, with A=0.48). */
public class HS075Test {

    private static final double A = 0.48;

    // Bounds: x1,x2 in [0,1200]; x3,x4 in [-A, +A]
    private static final double[] LB = { 0.0, 0.0, -A, -A };
    private static final double[] UB = { 1200.0, 1200.0,  A,  A };

    /** f(x) = 3*x1 + 1e-6*x1^3 + 2*x2 + (2e-6/3)*x2^3. */
    private static class TP75Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 4; }
        @Override public double value(RealVector X) {
            final double x1 = X.getEntry(0);
            final double x2 = X.getEntry(1);
            return 3.0*x1 + 1.0e-6*x1*x1*x1
                 + 2.0*x2 + (2.0e-6/3.0)*x2*x2*x2;
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x)  { throw new UnsupportedOperationException(); }
    }

    /** 2 disuguaglianze: g(x) ≥ 0. */
    private static class TP75Ineq extends InequalityConstraint {
        TP75Ineq() { super(new ArrayRealVector(new double[]{ 0.0, 0.0 })); }
        @Override public int dim() { return 4; }
        @Override public RealVector value(RealVector X) {
            final double x3 = X.getEntry(2), x4 = X.getEntry(3);
            final double g1 =  x4 - x3 + A;
            final double g2 =  x3 - x4 + A;
            return new ArrayRealVector(new double[]{ g1, g2 });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    /** 3 uguaglianze: h(x)=0. */
    private static class TP75Eq extends EqualityConstraint {
        TP75Eq() { super(new ArrayRealVector(new double[]{ 0.0, 0.0, 0.0 })); }
        @Override public int dim() { return 4; }
        @Override public RealVector value(RealVector X) {
            final double x1=X.getEntry(0), x2=X.getEntry(1),
                         x3=X.getEntry(2), x4=X.getEntry(3);

            final double h1 = 1000.0*(FastMath.sin(-x3 - 0.25) + FastMath.sin(-x4 - 0.25))
                            + 894.8-x1 ;

            final double h2 = 1000.0*(FastMath.sin( x3 - 0.25) + FastMath.sin( x3 - x4 - 0.25))
                            + 894.8-x2 ;

            final double h3 = 1000.0*(FastMath.sin( x4 - 0.25) + FastMath.sin( x4 - x3 - 0.25))
                            + 1294.8;

            return new ArrayRealVector(new double[]{ h1, h2, h3 });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    @Test
    public void testHS075() {
        final InitialGuess guess = new InitialGuess(new double[]{ 0.0, 0.0, 0.0, 0.0 });
        final SimpleBounds bounds = new SimpleBounds(LB, UB);

        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();

        final LagrangeSolution sol = optimizer.optimize(
                guess,
                new ObjectiveFunction(new TP75Obj()),
                new TP75Eq(),
                new TP75Ineq(),
                bounds
        );

        // FEX (TP75): 0.517441288686D+04
        HSProblemTestUtils.assertExpectedObjective(5174.41288686, sol);
    }
}