package org.hipparchus.optim.nonlinear.vector.constrained.sqp.hock_shittkowski;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.TwiceDifferentiableFunction;



import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

public class HS314Test {
    
    
    static final class HS314Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 2; }

        @Override public double value(RealVector x) {
           
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            
            double a = x1 - 2.0; // A=X(1)-.2D+1
            double b = x2 - 1.0; // B=X(2)-.1D+1
            
            // G1=(X(1)**2)/(-.4D+1)-X(2)**2+.1D+1
            double g1 = (x1 * x1 / -4.0) - (x2 * x2) + 1.0;
            
            // H1=X(1)-.2D+1*X(2)+.1D+1
            double h1 = x1 - 2.0 * x2 + 1.0;
            
            // FX=A**2+B**2+.4D-1/G1+H1**2/.2D+0
            return a * a + b * b + (0.04 / g1) + (h1 * h1 / 0.2);
        }

        @Override public RealVector gradient(RealVector x) {
           
            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);
            
            double g1 = (x1 * x1 / -4.0) - (x2 * x2) + 1.0;
            double h1 = x1 - 2.0 * x2 + 1.0;

            // GF(1)
            double g1_val = 2.0 * (x1 - 2.0 + (x1 * 0.01 / (g1 * g1)) + 5.0 * h1);
            
            // GF(2)
            double g2_val = 2.0 * (x2 - 1.0 + (x2 * 0.04 / (g1 * g1)) - 10.0 * h1);
            
            return new ArrayRealVector(new double[]{g1_val, g2_val}, false);
        }

        @Override public org.hipparchus.linear.RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided");
        }
    }

    private static double[] start() { 
       
        return new double[]{2.0, 2.0}; 
    }
    
    private static SimpleBounds bounds() {
        
        return new SimpleBounds(new double[]{1.0, 1.0},new double[]{Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY});
    }

    @Test
    public void testHS314() {
        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer() ;

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(start()),
                bounds(),
                new ObjectiveFunction(new HS314Obj())
        );

        double f = sol.getValue();
        
        final double fExpected = 0.16904;
        
        assertEquals(fExpected, f, 1.0e-3, "objective mismatch");
        
       
    }
}
