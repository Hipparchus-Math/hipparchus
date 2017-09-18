/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */

package org.hipparchus.analysis.solvers;

import org.hipparchus.analysis.QuinticFunction;
import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.differentiation.DSFactory;
import org.hipparchus.analysis.differentiation.DerivativeStructure;
import org.hipparchus.analysis.differentiation.UnivariateDifferentiableFunction;
import org.hipparchus.analysis.solvers.AllowedSolution;
import org.hipparchus.analysis.solvers.BracketingNthOrderBrentSolver;
import org.hipparchus.analysis.solvers.NewtonRaphsonSolver;
import org.hipparchus.analysis.solvers.UnivariateSolver;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.junit.Assert;
import org.junit.Test;

/**
 * Test case for {@link BracketingNthOrderBrentSolver bracketing n<sup>th</sup> order Brent} solver.
 *
 */
public final class BracketingNthOrderBrentSolverTest extends BaseSecantSolverAbstractTest {
    /** {@inheritDoc} */
    @Override
    protected UnivariateSolver getSolver() {
        return new BracketingNthOrderBrentSolver();
    }

    /** {@inheritDoc} */
    @Override
    protected int[] getQuinticEvalCounts() {
        return new int[] {1, 3, 8, 1, 9, 4, 8, 1, 12, 1, 16};
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testInsufficientOrder1() {
        new BracketingNthOrderBrentSolver(1.0e-10, 1);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testInsufficientOrder2() {
        new BracketingNthOrderBrentSolver(1.0e-10, 1.0e-10, 1);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testInsufficientOrder3() {
        new BracketingNthOrderBrentSolver(1.0e-10, 1.0e-10, 1.0e-10, 1);
    }

    @Test
    public void testConstructorsOK() {
        Assert.assertEquals(2, new BracketingNthOrderBrentSolver(1.0e-10, 2).getMaximalOrder());
        Assert.assertEquals(2, new BracketingNthOrderBrentSolver(1.0e-10, 1.0e-10, 2).getMaximalOrder());
        Assert.assertEquals(2, new BracketingNthOrderBrentSolver(1.0e-10, 1.0e-10, 1.0e-10, 2).getMaximalOrder());
    }

    @Test
    public void testConvergenceOnFunctionAccuracy() {
        BracketingNthOrderBrentSolver solver =
                new BracketingNthOrderBrentSolver(1.0e-12, 1.0e-10, 0.001, 3);
        QuinticFunction f = new QuinticFunction();
        double result = solver.solve(20, f, 0.2, 0.9, 0.4, AllowedSolution.BELOW_SIDE);
        Assert.assertEquals(0, f.value(result), solver.getFunctionValueAccuracy());
        Assert.assertTrue(f.value(result) <= 0);
        Assert.assertTrue(result - 0.5 > solver.getAbsoluteAccuracy());
        result = solver.solve(20, f, -0.9, -0.2,  -0.4, AllowedSolution.ABOVE_SIDE);
        Assert.assertEquals(0, f.value(result), solver.getFunctionValueAccuracy());
        Assert.assertTrue(f.value(result) >= 0);
        Assert.assertTrue(result + 0.5 < -solver.getAbsoluteAccuracy());
    }

    @Test
    public void testIssue716() {
        BracketingNthOrderBrentSolver solver =
                new BracketingNthOrderBrentSolver(1.0e-12, 1.0e-10, 1.0e-22, 5);
        UnivariateFunction sharpTurn = new UnivariateFunction() {
            public double value(double x) {
                return (2 * x + 1) / (1.0e9 * (x + 1));
            }
        };
        double result = solver.solve(100, sharpTurn, -0.9999999, 30, 15, AllowedSolution.RIGHT_SIDE);
        Assert.assertEquals(0, sharpTurn.value(result), solver.getFunctionValueAccuracy());
        Assert.assertTrue(sharpTurn.value(result) >= 0);
        Assert.assertEquals(-0.5, result, 1.0e-10);
    }

    @Test
    public void testToleranceLessThanUlp() {
        // function that is never zero
        UnivariateFunction f = (x) -> x < 2.1 ? -1 : 1;
        // tolerance less than 1 ulp(x)
        UnivariateSolver solver = new BracketingNthOrderBrentSolver(0, 1e-18, 0, 5);

        // make sure it doesn't throw a maxIterations exception
        double result = solver.solve(100, f, 0.0, 5.0);
        Assert.assertEquals(2.1, result, 0.0);
    }

    @Test
    public void testFasterThanNewton() {
        // the following test functions come from Beny Neta's paper:
        // "Several New Methods for solving Equations"
        // intern J. Computer Math Vol 23 pp 265-282
        // available here: http://www.math.nps.navy.mil/~bneta/SeveralNewMethods.PDF
        // the reference roots have been computed by the Dfp solver to more than
        // 80 digits and checked with emacs (only the first 20 digits are reproduced here)
        compare(new TestFunction(0.0, -2, 2) {
            @Override
            public DerivativeStructure value(DerivativeStructure x) {
                return x.sin().subtract(x.multiply(0.5));
            }
        });
        compare(new TestFunction(6.3087771299726890947, -5, 10) {
            @Override
            public DerivativeStructure value(DerivativeStructure x) {
                return x.pow(5).add(x).subtract(10000);
            }
        });
        compare(new TestFunction(9.6335955628326951924, 0.001, 10) {
            @Override
            public DerivativeStructure value(DerivativeStructure x) {
                return x.sqrt().subtract(x.reciprocal()).subtract(3);
            }
        });
        compare(new TestFunction(2.8424389537844470678, -5, 5) {
            @Override
            public DerivativeStructure value(DerivativeStructure x) {
                return x.exp().add(x).subtract(20);
            }
        });
        compare(new TestFunction(8.3094326942315717953, 0.001, 10) {
            @Override
            public DerivativeStructure value(DerivativeStructure x) {
                return x.log().add(x.sqrt()).subtract(5);
            }
        });
        compare(new TestFunction(1.4655712318767680266, -0.5, 1.5) {
            @Override
            public DerivativeStructure value(DerivativeStructure x) {
                return x.subtract(1).multiply(x).multiply(x).subtract(1);
            }
        });

    }

    private void compare(TestFunction f) {
        compare(f, f.getRoot(), f.getMin(), f.getMax());
    }

    private void compare(final UnivariateDifferentiableFunction f,
                         double root, double min, double max) {
        NewtonRaphsonSolver newton = new NewtonRaphsonSolver(1.0e-12);
        BracketingNthOrderBrentSolver bracketing =
                new BracketingNthOrderBrentSolver(1.0e-12, 1.0e-12, 1.0e-18, 5);
        double resultN;
        try {
            resultN = newton.solve(100, f, min, max);
        } catch (MathIllegalStateException tmee) {
            resultN = Double.NaN;
        }
        double resultB;
        try {
            resultB = bracketing.solve(100, f, min, max);
        } catch (MathIllegalStateException tmee) {
            resultB = Double.NaN;
        }
        Assert.assertEquals(root, resultN, newton.getAbsoluteAccuracy());
        Assert.assertEquals(root, resultB, bracketing.getAbsoluteAccuracy());

        // bracketing solver evaluates only function value, we set the weight to 1
        final int weightedBracketingEvaluations = bracketing.getEvaluations();

        // Newton-Raphson solver evaluates both function value and derivative, we set the weight to 2
        final int weightedNewtonEvaluations = 2 * newton.getEvaluations();

        Assert.assertTrue(weightedBracketingEvaluations < weightedNewtonEvaluations);

    }

    private static abstract class TestFunction implements UnivariateDifferentiableFunction {

        private final double root;
        private final double min;
        private final double max;

        protected TestFunction(final double root, final double min, final double max) {
            this.root = root;
            this.min  = min;
            this.max  = max;
        }

        public double getRoot() {
            return root;
        }

        public double getMin() {
            return min;
        }

        public double getMax() {
            return max;
        }

        public double value(final double x) {
            return value(new DSFactory(0, 0).constant(x)).getValue();
        }

        public abstract DerivativeStructure value(final DerivativeStructure t);

    }

}
