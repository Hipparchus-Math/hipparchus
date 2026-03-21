package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

/** HS TP74 (Schittkowski). KN1=1 -> A=0.55. 2 ineq lineari, 3 eq non lineari. */
public class HS074Test {

    private static final double A = 0.55; // A(1)

    // Bounds: x1,x2 in [0,1200]; x3,x4 in [-A, +A]
    private static final double[] LB = { 0.0, 0.0, -A, -A };
    private static final double[] UB = { 1200.0, 1200.0,  A,  A };

    /** f(x) = 3*x1 + 1e-6*x1^3 + 2*x2 + (2e-6/3)*x2^3. */
    private static class TP74Obj extends TwiceDifferentiableFunction {
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
    private static class TP74Ineq extends InequalityConstraint {
        TP74Ineq() { super(new ArrayRealVector(new double[]{ 0.0, 0.0 })); }
        @Override public int dim() { return 4; }
        @Override public RealVector value(RealVector X) {
            final double x3 = X.getEntry(2), x4 = X.getEntry(3);
            final double g1 =  x4 - x3 + A;
            final double g2 =  x3 - x4 + A;
            return new ArrayRealVector(new double[]{ g1, g2 });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    /** 3 uguaglianze accorpate: h(x) = 0. */
    private static class TP74Eq extends EqualityConstraint {
        TP74Eq() { super(new ArrayRealVector(new double[]{ 0.0, 0.0, 0.0 })); }
        @Override public int dim() { return 4; }
        @Override public RealVector value(RealVector X) {
            final double x1=X.getEntry(0), x2=X.getEntry(1),
                         x3=X.getEntry(2), x4=X.getEntry(3);

            final double h1 = 1000.0*(FastMath.sin(-x3 - 0.25) + FastMath.sin(-x4 - 0.25))
                            + 894.8 - x1;

            final double h2 = 1000.0*(FastMath.sin( x3 - 0.25) + FastMath.sin( x3 - x4 - 0.25))
                            + 894.8 - x2;

            final double h3 = 1000.0*(FastMath.sin( x4 - 0.25) + FastMath.sin( x4 - x3 - 0.25))
                            - 1294.8;

            return new ArrayRealVector(new double[]{ h1, h2, h3 });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    @Test
    public void testHS074() {
        final InitialGuess guess = new InitialGuess(new double[]{ 0.0, 0.0, 0.0, 0.0 });
        final SimpleBounds bounds = new SimpleBounds(LB, UB);

        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();

        final LagrangeSolution sol = optimizer.optimize(
                guess,
                new ObjectiveFunction(new TP74Obj()),
                new TP74Eq(),      // 3 equalities
                new TP74Ineq(),    // 2 inequalities (±A on x4-x3)
                bounds
        );

        // FEX (TP74, KN1=1): 0.512649810934D+04
        HSProblemTestUtils.assertExpectedObjective(5126.49810934, sol);
    }
}