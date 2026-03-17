/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
 */
package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.special.Erf;
import org.hipparchus.util.FastMath;
import static org.junit.jupiter.api.Assertions.assertEquals;
import org.junit.jupiter.api.Test;

/** HS TP69 (Schittkowski). KN1=2 → A=0.1, B=1000, Z=4. */
public class HS069Test {

    // Costanti dal COMMON /D68/ con KN1=2
    private static final double a = 0.1;      // A(2)
    private static final double b = 1000.0;   // B(2)
    private static final double z = 4.0;      // Z(2)

    // Bounds MODE=1: XL(1)=1e-4, XL(2)=0, XL(3)=0, XL(4)=0; XU(1,2)=100; XU(3,4)=2
    private static final double[] LB = { 1.0e-4, 0.0, 0.0, 1.0e-4 };
    private static final double[] UB = { 1.0, 100.0, 2.0, 2.0 };

    /** Φ(z): CDF normale standard = 0.5*(1+erf(z/√2)). */
    private static double phi(double t) {
        return 0.5 * (1.0 + Erf.erf(t / FastMath.sqrt(2.0)));
    }

    /** f(x) = (a*z - x4*(b*(exp(x1)-1) - x3)/(exp(x1)-1 + x4)) / x1. */
    private static class TP69Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 4; }
        @Override public double value(RealVector X) {
            final double x1 = X.getEntry(0);
            final double x3 = X.getEntry(2);
            final double x4 = X.getEntry(3);
            final double v1 = FastMath.exp(x1) - 1.0;     // V1 = exp(x1) - 1
            return (a * z - x4 * (b * v1 - x3) / (v1 + x4)) / x1;
        }
        @Override public RealVector gradient(RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(RealVector x)  { throw new UnsupportedOperationException(); }
    }

    /** 2 uguaglianze (INDEX1): g1(x)=0, g2(x)=0 con MDNORD ≡ Φ. */
    private static class TP69Eq extends EqualityConstraint {
        TP69Eq() { super(new ArrayRealVector(new double[]{ 0.0, 0.0 })); }
        @Override public int dim() { return 4; }
        @Override public RealVector value(RealVector X) {
            final double x2 = X.getEntry(1);
            final double x3 = X.getEntry(2);
            final double x4 = X.getEntry(3);
            // g1 = x3 - 2*Phi(-x2)
            final double g1 = x3 - 2.0 * phi(-x2);
            // g2 = x4 - Phi(-x2 + sqrt(z)) - Phi(-x2 - sqrt(z)), con z=4 => sqrt(z)=2
            final double rtZ = FastMath.sqrt(z);
            final double g2 = x4 - phi(-x2 + rtZ) - phi(-x2 - rtZ);
            return new ArrayRealVector(new double[]{ g1, g2 });
        }
        @Override public RealMatrix jacobian(RealVector x) { throw new UnsupportedOperationException(); }
    }
    
    @Test
    public void testHS069() {
        final InitialGuess guess = new InitialGuess(new double[]{ 1.0, 1.0, 1.0, 1.0 });
        final SimpleBounds bounds = new SimpleBounds(LB, UB);
         SQPOption option=new SQPOption();
        option.setGradientMode(GradientMode.FORWARD);
        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();

        final LagrangeSolution sol = optimizer.optimize(
                option,
                guess,
                new ObjectiveFunction(new TP69Obj()),
                new TP69Eq(),
               bounds
        );

        // FEX (TP69): -0.956712887064D+03
        HSProblemTestUtils.assertExpectedObjective(-956.712887064, sol);
    }
}
