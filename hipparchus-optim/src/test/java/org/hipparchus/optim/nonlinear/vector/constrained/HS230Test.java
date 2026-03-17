/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 *
 * The Hipparchus project licenses this file to You under the Apache License,
 * Version 2.0.
 */

package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/**
 * HS230 (TP230)
 *
 * N    = 2
 * NILI = 0
 * NINL = 2 (two nonlinear inequalities)
 * NELI = 0
 * NENL = 0
 *
 * Objective:
 *   f(x) = x2
 *
 * Fortran constraints (G(x) ≥ 0):
 *   G1(x) = -2 x1^2 + x1^3 + x2
 *   G2(x) = -2 (1 - x1)^2 + (1 - x1)^3 + x2
 *
 * Reference solution (MODE=1):
 *   x*  = (0.5, 0.375)
 *   f*  = 0.375
 */
public class HS230Test {

    private static final int DIM      = 2;
    private static final int NUM_INEQ = 2;
    private static final int NUM_EQ   = 0;

    // -------------------------------------------------------------------------
    // Objective
    // -------------------------------------------------------------------------
    private static class HS230Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public double value(RealVector x) {
            return x.getEntry(1); // x2
        }

        @Override
        public RealVector gradient(RealVector x) {
            return new ArrayRealVector(new double[]{0.0, 1.0}, false);
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            return new Array2DRowRealMatrix(DIM, DIM);
        }
    }

    // -------------------------------------------------------------------------
    // Inequality constraints G(x) ≥ 0
    // -------------------------------------------------------------------------
    private static class HS230Ineq extends InequalityConstraint {

        HS230Ineq() {
            // RHS = 0 for both inequalities (G >= 0)
            super(new ArrayRealVector(new double[NUM_INEQ]));
        }

        @Override
        public int dim() {
            return DIM;
        }

        @Override
        public RealVector value(RealVector x) {

            double x1 = x.getEntry(0);
            double x2 = x.getEntry(1);

            // Fortran:
            // G1 = -2*x1^2 + x1^3 + x2
            double G1 = -2.0 * x1 * x1 + x1 * x1 * x1 + x2;

            // G2 = -2*(1-x1)^2 + (1-x1)^3 + x2
            double y  = 1.0 - x1;
            double G2 = -2.0 * y * y + y * y * y + x2;

            return new ArrayRealVector(new double[]{G1, G2}, false);
        }

        @Override
        public RealMatrix jacobian(RealVector x) {

            double x1 = x.getEntry(0);

            // dG1/dx1 = -4*x1 + 3*x1^2
            double dG1dx1 = -4.0 * x1 + 3.0 * x1 * x1;
            // dG1/dx2 = 1
            double dG1dx2 = 1.0;

            double y = 1.0 - x1;
            // dG2/dx1 = 4*(1-x1) - 3*(1-x1)^2
            double dG2dx1 = 4.0 * y - 3.0 * y * y;
            // dG2/dx2 = 1
            double dG2dx2 = 1.0;

            RealMatrix J = new Array2DRowRealMatrix(NUM_INEQ, DIM);
            J.setEntry(0, 0, dG1dx1);
            J.setEntry(0, 1, dG1dx2);

            J.setEntry(1, 0, dG2dx1);
            J.setEntry(1, 1, dG2dx2);

            return J;
        }
    }

    // -------------------------------------------------------------------------
    // Test
    // -------------------------------------------------------------------------
    @Test
    public void testHS230_optimization() {

        // Initial guess (MODE=1): X(1)=0, X(2)=0
        double[] x0 = new double[]{0.0, 0.0};

        SQPOptimizerS2 opt = HSProblemTestUtils.newOptimizer();

        LagrangeSolution sol = opt.optimize(
                new InitialGuess(x0),
                new ObjectiveFunction(new HS230Obj()),
                null,              // no equalities
                new HS230Ineq(),   // 2 nonlinear inequalities
                null               // no explicit bounds
        );

        double f = sol.getValue();

        final double val = 0.375;
        HSProblemTestUtils.assertExpectedObjective(val, sol);
    }
}
