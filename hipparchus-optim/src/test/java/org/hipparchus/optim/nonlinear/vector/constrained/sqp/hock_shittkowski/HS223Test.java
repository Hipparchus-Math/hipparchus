/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License,
 * Version 2.0.
 */

package org.hipparchus.optim.nonlinear.vector.constrained.sqp.hock_shittkowski;

import org.hipparchus.linear.Array2DRowRealMatrix;
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
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * HS223 (TP223)
 *
 * N    = 2
 * NILI = 0
 * NINL = 2
 * NELI = 0
 * NENL = 0
 *
 * Bounds:
 *   0 <= x1 <= 1
 *   0 <= x2 <= 10
 *
 * Objective:
 *   f(x) = -x1
 *
 * Inequalities (Fortran G, kept with same sign convention):
 *
 *   G1(x) = exp(exp(x1))
 *   G2(x) = x2 - exp(exp(x1))
 *
 * Reference solution (LEX = .TRUE.):
 *
 *   x* = ( log(log 10), 10 )
 *   f* = -log(log 10)
 */
public class HS223Test {

    private static final int DIM      = 2;
    private static final int NUM_INEQ = 2;

    // -------------------------------------------------------------------------
    // Objective f = -x1
    // -------------------------------------------------------------------------
    private static class HS223Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            return -x.getEntry(0);
        }

        @Override
        public RealVector gradient(RealVector x) {
            // df/dx1 = -1, df/dx2 = 0
            return new ArrayRealVector(new double[]{-1.0, 0.0}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            // Zero Hessian – BFGS will build curvature
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Inequalities:
    //   G1 = exp(exp(x1))
    //   G2 = x2 - exp(exp(x1))
    // -------------------------------------------------------------------------
    private static class HS223Ineq extends InequalityConstraint {

        HS223Ineq() {
            // RHS = 0 for both inequalities; we pass G(x) directly.
            super(new ArrayRealVector(new double[NUM_INEQ]));
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {
            final double x1 = x.getEntry(0);
            final double x2 = x.getEntry(1);

            final double e1 = FastMath.exp(x1);
            final double e2 = FastMath.exp(e1); // exp(exp(x1))

            double g1 = e2;
            double g2 = x2 - e2;

            return new ArrayRealVector(new double[]{g1, g2}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {
            final double x1 = x.getEntry(0);

            final double e1 = FastMath.exp(x1);
            final double e2 = FastMath.exp(e1); // exp(exp(x1))
            final double dCommon = e1 * e2;     // d/dx1 exp(exp(x1))

            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);

            // G1 = exp(exp(x1))
            // dG1/dx1 = exp(x1) * exp(exp(x1)), dG1/dx2 = 0
            J.setEntry(0, 0, dCommon);
            J.setEntry(0, 1, 0.0);

            // G2 = x2 - exp(exp(x1))
            // dG2/dx1 = -exp(x1) * exp(exp(x1)), dG2/dx2 = 1
            J.setEntry(1, 0, -dCommon);
            J.setEntry(1, 1, 1.0);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS223_optimization() {

        // Initial guess (Fortran mode 1): x1 = 0.1, x2 = 3.3
        double[] x0 = new double[]{0.1, 3.3};

        // Bounds: 0 <= x1 <= 1, 0 <= x2 <= 10
        double[] lower = new double[]{0.0, 0.0};
        double[] upper = new double[]{1.0, 10.0};
        SimpleBounds bounds = new SimpleBounds(lower, upper);

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS223Obj()),
                null,              // no equalities
                new HS223Ineq(),   // 2 inequalities
                bounds
        );

        double f = sol.getValue();

        final double fExpected = -FastMath.log(FastMath.log(10.0));
        final double tol = 1.0e-6 * (FastMath.abs(fExpected) + 1.0);

        assertEquals(fExpected, f, tol);
    }
}
