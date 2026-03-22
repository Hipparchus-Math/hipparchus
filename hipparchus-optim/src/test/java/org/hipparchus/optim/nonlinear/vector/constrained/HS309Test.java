package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS309Test {

    static final class HS309Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 2; }

       @Override
public double value(RealVector x) {
    final double x1 = x.getEntry(0);
    final double x2 = x.getEntry(1);
    double p = 1.41  * Math.pow(x1, 4)   // .141D+1
             - 12.76 * Math.pow(x1, 3)   // .1276D+2
             + 39.91 * Math.pow(x1, 2)   // .3991D+2
             - 51.93 * x1                // .5193D+2
             + 24.37;                    // .2437D+2
    double q = x2 - 3.9;                 // .39D+1
    return p + q*q;
}

@Override
public RealVector gradient(RealVector x) {
    final double x1 = x.getEntry(0);
    final double x2 = x.getEntry(1);
    double g1 = 5.64  * Math.pow(x1, 3)  // .564D+1
              - 38.28 * Math.pow(x1, 2)  // .3828D+2
              + 79.82 * x1               // .7982D+2
              - 51.93;                   // .5193D+2
    double g2 = 2.0 * x2 - 7.8;          // 2*(x2 - 3.9)  -> .2D+1, .78D+1
    return new ArrayRealVector(new double[]{g1, g2}, false);
}


        @Override
        public RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided");
        }
    }

    private LagrangeSolution solve() {
        SQPOptimizerS2 opt = new SQPOptimizerS2();
        opt.setDebugPrinter(System.out::println); // richiesto

        // Start: X(1)=0, X(2)=0 
        double[] x0 = {0.0, 0.0};

       
        return opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS309Obj())
        );
    }

    @Test
    public void testHS309() {
        final double fExpected = -3.9871708;
        LagrangeSolution sol = solve();
        double f = sol.getValue();
        assertEquals(fExpected, f, 1.0e-4 * (Math.abs(fExpected) + 1.0), "objective mismatch");
    }
}
