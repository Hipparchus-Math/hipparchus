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

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

public class HS055Test {

    private static final class HS055Obj extends TwiceDifferentiableFunction {
        @Override public int dim() { return 6; }

        @Override
        public double value(final RealVector x) {
          
      final double x1=x.getEntry(0);
      final double x2=x.getEntry(1);
      final double x3=x.getEntry(2);
      final double x4=x.getEntry(3);
      final double x5=x.getEntry(4);
      final double x6=x.getEntry(5);
        // X14 = X(1)*X(4)
      //IF (X14.GT.1.0D1) X14 = 1.0D1   
      //FX=X(1)+2.D0*X(2)+4.D0*X(5)+DEXP(X14)   
     
       final double x14 = FastMath.min(10.0, x1*x4);
            return x1 + 2.0 * x2 + 4.0 * x5 + FastMath.exp(x14);
        }

        @Override public RealVector gradient(final RealVector x) { throw new UnsupportedOperationException(); }
        @Override public RealMatrix hessian(final RealVector x) { throw new UnsupportedOperationException(); }
    }

    private static final class HS055Eq extends EqualityConstraint {
        HS055Eq() { super(new ArrayRealVector(new double[6])); }

        @Override
        public RealVector value(final RealVector x) {
            return new ArrayRealVector(new double[] {
                x.getEntry(0) + 2.0 * x.getEntry(1) + 5.0 * x.getEntry(4) - 6.0,
                x.getEntry(0) + x.getEntry(1) + x.getEntry(2) - 3.0,
                x.getEntry(3) + x.getEntry(4) + x.getEntry(5) - 2.0,
                x.getEntry(0) + x.getEntry(3) - 1.0,
                x.getEntry(1) + x.getEntry(4) - 2.0,
                x.getEntry(2) + x.getEntry(5) - 2.0
            }, false);
        }

        @Override
        public RealMatrix jacobian(final RealVector x) {
            return new Array2DRowRealMatrix(new double[][] {
                {1.0, 2.0, 0.0, 0.0, 5.0, 0.0},
                {1.0, 1.0, 1.0, 0.0, 0.0, 0.0},
                {0.0, 0.0, 0.0, 1.0, 1.0, 1.0},
                {1.0, 0.0, 0.0, 1.0, 0.0, 0.0},
                {0.0, 1.0, 0.0, 0.0, 1.0, 0.0},
                {0.0, 0.0, 1.0, 0.0, 0.0, 1.0}
            }, false);
        }

        @Override public int dim() { return 6; }
    }

    @Test
    void testHS055() {
        final SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        final LagrangeSolution sol = optimizer.optimize(
                new InitialGuess(new double[] {1.0, 2.0, 0.0, 0.0, 0.0, 2.0}),
                new ObjectiveFunction(new HS055Obj()),
                new HS055Eq(),
                new SimpleBounds(
                        new double[] {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                        new double[] {1.0, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, 1.0,
                                      Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY}
                )
        );

        HSProblemTestUtils.assertExpectedObjective(19.0 / 3.0, sol);
    }
}
