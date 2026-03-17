/*
 * Licensed to the Hipparchus project under one or more contributor license agreements...
 */
package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

/** HS TP84 (hs84). 5 variabili, 3 vincoli a intervallo scritti come 6 inequality. */
public class HS084Test {

    // a[1..21] (uso indicizzazione 1-based: a[0] dummy)
    private static final double[] a = {
        0.0,
        -24345.0,
        -8720288.849,
         150512.5253,
        -156.6950325,
         476470.3222,
         729482.8271,
        -145421.402,
         2931.1506,
        -40.427932,
         5106.192,
         15711.36,
        -155011.1084,
         4360.53352,
         12.9492344,
         10236.884,
         13176.786,
        -326669.5104,
         7390.68412,
        -27.8986976,
         16643.076,
         30988.146
    };

    // Bounds (l/u nel tuo modello)
    private static final double[] LB = { 0.0, 1.2, 20.0, 9.0, 6.5 };
    private static final double[] UB = { 1000.0, 2.4, 60.0, 9.3, 7.0 };

    // Guess xi
    private static final double[] X0 = { 2.52, 2.0, 37.5, 9.25, 6.8 };

    /** f(x) = -a1 - a2*x1 - a3*x1*x2 - a4*x1*x3 - a5*x1*x4 - a6*x1*x5. */
    private static class TP84Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 5; }
        @Override public double value(RealVector X) {
            final double x1=X.getEntry(0), x2=X.getEntry(1),
                         x3=X.getEntry(2), x4=X.getEntry(3), x5=X.getEntry(4);
            return -a[1]
                   - a[2]*x1
                   - a[3]*x1*x2
                   - a[4]*x1*x3
                   - a[5]*x1*x4
                   - a[6]*x1*x5;
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x)  { throw new UnsupportedOperationException(); }
    }

    /**
     * Vincoli a intervallo scritti come g(x) >= 0 (la tua convenzione):
     * 0 <= V1 <= 294000,  0 <= V2 <= 294000,  0 <= V3 <= 277200
     * con:
     *  V1 = a7*x1 + a8*x1*x2 + a9*x1*x3 + a10*x1*x4 + a11*x1*x5
     *  V2 = a12*x1 + a13*x1*x2 + a14*x1*x3 + a15*x1*x4 + a16*x1*x5
     *  V3 = a17*x1 + a18*x1*x2 + a19*x1*x3 + a20*x1*x4 + a21*x1*x5
     */
    private static class TP84Ineq extends InequalityConstraint {
        TP84Ineq() { super(new ArrayRealVector(new double[]{0,0,0,0,0,0})); }
        @Override public int dim() { return 5; }
        @Override public RealVector value(RealVector X) {
            final double x1=X.getEntry(0), x2=X.getEntry(1),
                         x3=X.getEntry(2), x4=X.getEntry(3), x5=X.getEntry(4);

            final double V1 = x1*(a[7]  + a[8]*x2  + a[9]*x3  + a[10]*x4 + a[11]*x5);
            final double V2 = x1*(a[12] + a[13]*x2 + a[14]*x3 + a[15]*x4 + a[16]*x5);
            final double V3 = x1*(a[17] + a[18]*x2 + a[19]*x3 + a[20]*x4 + a[21]*x5);

            return new ArrayRealVector(new double[]{
                V1, V2, V3,
                294000.0 - V1,
                294000.0 - V2,
                277200.0 - V3
            });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
    }

    @Test
    public void testHS084() {
        final InitialGuess guess = new InitialGuess(X0);
        final SimpleBounds bounds = new SimpleBounds(LB, UB);

        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        optimizer.setDebugPrinter(System.out::println);

        final LagrangeSolution sol = optimizer.optimize(
            guess,
            new ObjectiveFunction(new TP84Obj()),
            new TP84Ineq(),
            bounds
        );

        // best known objective = -5280335.133
        HSProblemTestUtils.assertExpectedObjective(-5280335.133, sol);
    }
}
