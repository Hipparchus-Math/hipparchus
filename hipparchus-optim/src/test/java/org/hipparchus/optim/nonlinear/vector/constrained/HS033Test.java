
package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS033Test {

    // Bounds: x1 >= 0, x2 >= 0, 0 <= x3 <= 5 ; niente upper per x1,x2 -> mettiamo BIG
    private static final double BIG = 1.0e12;
    private static final double[] LB = { 0.0, 0.0, 0.0 };
    private static final double[] UB = { BIG, BIG, 5.0 };

    /** f(x) = (x1-1)(x1-2)(x1-3) + x3. */
    private static class TP33Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 3; }
        @Override public double value(RealVector X) {
            final double x1 = X.getEntry(0);
            final double x3 = X.getEntry(2);
            return (x1 - 1.0)*(x1 - 2.0)*(x1 - 3.0) + x3;
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x)  { throw new UnsupportedOperationException(); }
    }

    /** Due disequazioni in forma g(x) >= 0. */
    private static class TP33Ineq extends InequalityConstraint {
        TP33Ineq() { super(new ArrayRealVector(new double[]{0.0, 0.0})); }
        @Override public int dim() { return 3; }

        @Override public RealVector value(RealVector X) {
            final double x1 = X.getEntry(0);
            final double x2 = X.getEntry(1);
            final double x3 = X.getEntry(2);

            final double g1 = x3*x3 - x1*x1 - x2*x2;                 // >= 0
            final double g2 = x1*x1 + x2*x2 + x3*x3 - 4.0;           // >= 0

            return new ArrayRealVector(new double[]{ g1, g2 });
        }

        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    @Test
    public void testHS033() {
        final InitialGuess guess = new InitialGuess(new double[]{ 0.1, 0.1, 3.0 });
        final SimpleBounds bounds = new SimpleBounds(LB, UB);

        final SQPOptimizerS2 opt = new SQPOptimizerS2();
        if (Boolean.getBoolean("hipparchus.debug.sqp")) {
            opt.setDebugPrinter(System.out::println);
        }
        SQPOption sqpOption=new SQPOption();
        sqpOption.setMaxLineSearchIteration(20);
        sqpOption.setEps(10e-11);
        final LagrangeSolution sol = opt.optimize(
            guess,
            new ObjectiveFunction(new TP33Obj()),
            new TP33Ineq(),
            bounds,
            sqpOption
        );

        final double expected = Math.sqrt(2.0) - 6.0;
        assertEquals(expected, sol.getValue(), 1e-4);
    }
}
