/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
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
package org.hipparchus.linear;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.IterationEvent;
import org.hipparchus.util.IterationListener;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.assertThrows;

public class ConjugateGradientTest {

    @Test
    public void testNonSquareOperator() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final Array2DRowRealMatrix a = new Array2DRowRealMatrix(2, 3);
            final IterativeLinearSolver solver;
            solver = new ConjugateGradient(10, 0., false);
            final ArrayRealVector b = new ArrayRealVector(a.getRowDimension());
            final ArrayRealVector x = new ArrayRealVector(a.getColumnDimension());
            solver.solve(a, b, x);
        });
    }

    @Test
    public void testDimensionMismatchRightHandSide() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final Array2DRowRealMatrix a = new Array2DRowRealMatrix(3, 3);
            final IterativeLinearSolver solver;
            solver = new ConjugateGradient(10, 0., false);
            final ArrayRealVector b = new ArrayRealVector(2);
            final ArrayRealVector x = new ArrayRealVector(3);
            solver.solve(a, b, x);
        });
    }

    @Test
    public void testDimensionMismatchSolution() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final Array2DRowRealMatrix a = new Array2DRowRealMatrix(3, 3);
            final IterativeLinearSolver solver;
            solver = new ConjugateGradient(10, 0., false);
            final ArrayRealVector b = new ArrayRealVector(3);
            final ArrayRealVector x = new ArrayRealVector(2);
            solver.solve(a, b, x);
        });
    }

    @Test
    public void testNonPositiveDefiniteLinearOperator() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final Array2DRowRealMatrix a = new Array2DRowRealMatrix(2, 2);
            a.setEntry(0, 0, -1.);
            a.setEntry(0, 1, 2.);
            a.setEntry(1, 0, 3.);
            a.setEntry(1, 1, 4.);
            final IterativeLinearSolver solver;
            solver = new ConjugateGradient(10, 0., true);
            final ArrayRealVector b = new ArrayRealVector(2);
            b.setEntry(0, -1.);
            b.setEntry(1, -1.);
            final ArrayRealVector x = new ArrayRealVector(2);
            solver.solve(a, b, x);
        });
    }

    @Test
    public void testUnpreconditionedSolution() {
        final int n = 5;
        final int maxIterations = 100;
        final RealLinearOperator a = new HilbertMatrix(n);
        final InverseHilbertMatrix ainv = new InverseHilbertMatrix(n);
        final IterativeLinearSolver solver;
        solver = new ConjugateGradient(maxIterations, 1E-10, true);
        final RealVector b = new ArrayRealVector(n);
        for (int j = 0; j < n; j++) {
            b.set(0.);
            b.setEntry(j, 1.);
            final RealVector x = solver.solve(a, b);
            for (int i = 0; i < n; i++) {
                final double actual = x.getEntry(i);
                final double expected = ainv.getEntry(i, j);
                final double delta = 1E-10 * FastMath.abs(expected);
                final String msg = String.format("entry[%d][%d]", i, j);
                Assertions.assertEquals(expected, actual, delta, msg);
            }
        }
    }

    @Test
    public void testUnpreconditionedInPlaceSolutionWithInitialGuess() {
        final int n = 5;
        final int maxIterations = 100;
        final RealLinearOperator a = new HilbertMatrix(n);
        final InverseHilbertMatrix ainv = new InverseHilbertMatrix(n);
        final IterativeLinearSolver solver;
        solver = new ConjugateGradient(maxIterations, 1E-10, true);
        final RealVector b = new ArrayRealVector(n);
        for (int j = 0; j < n; j++) {
            b.set(0.);
            b.setEntry(j, 1.);
            final RealVector x0 = new ArrayRealVector(n);
            x0.set(1.);
            final RealVector x = solver.solveInPlace(a, b, x0);
            Assertions.assertSame(x0, x, "x should be a reference to x0");
            for (int i = 0; i < n; i++) {
                final double actual = x.getEntry(i);
                final double expected = ainv.getEntry(i, j);
                final double delta = 1E-10 * FastMath.abs(expected);
                final String msg = String.format("entry[%d][%d)", i, j);
                Assertions.assertEquals(expected, actual, delta, msg);
            }
        }
    }

    @Test
    public void testUnpreconditionedSolutionWithInitialGuess() {
        final int n = 5;
        final int maxIterations = 100;
        final RealLinearOperator a = new HilbertMatrix(n);
        final InverseHilbertMatrix ainv = new InverseHilbertMatrix(n);
        final IterativeLinearSolver solver;
        solver = new ConjugateGradient(maxIterations, 1E-10, true);
        final RealVector b = new ArrayRealVector(n);
        for (int j = 0; j < n; j++) {
            b.set(0.);
            b.setEntry(j, 1.);
            final RealVector x0 = new ArrayRealVector(n);
            x0.set(1.);
            final RealVector x = solver.solve(a, b, x0);
            Assertions.assertNotSame(x0, x, "x should not be a reference to x0");
            for (int i = 0; i < n; i++) {
                final double actual = x.getEntry(i);
                final double expected = ainv.getEntry(i, j);
                final double delta = 1E-10 * FastMath.abs(expected);
                final String msg = String.format("entry[%d][%d]", i, j);
                Assertions.assertEquals(expected, actual, delta, msg);
                Assertions.assertEquals(1., x0.getEntry(i), Math.ulp(1.), msg);
            }
        }
    }

    /**
     * Check whether the estimate of the (updated) residual corresponds to the
     * exact residual. This fails to be true for a large number of iterations,
     * due to the loss of orthogonality of the successive search directions.
     * Therefore, in the present test, the number of iterations is limited.
     */
    @Test
    public void testUnpreconditionedResidual() {
        final int n = 10;
        final int maxIterations = n;
        final RealLinearOperator a = new HilbertMatrix(n);
        final ConjugateGradient solver;
        solver = new ConjugateGradient(maxIterations, 1E-15, true);
        final RealVector r = new ArrayRealVector(n);
        final RealVector x = new ArrayRealVector(n);
        final IterationListener listener = new IterationListener() {

            public void terminationPerformed(final IterationEvent e) {
                // Do nothing
            }

            public void iterationStarted(final IterationEvent e) {
                // Do nothing
            }

            public void iterationPerformed(final IterationEvent e) {
                final IterativeLinearSolverEvent evt;
                evt = (IterativeLinearSolverEvent) e;
                RealVector v = evt.getResidual();
                r.setSubVector(0, v);
                v = evt.getSolution();
                x.setSubVector(0, v);
            }

            public void initializationPerformed(final IterationEvent e) {
                // Do nothing
            }
        };
        solver.getIterationManager().addIterationListener(listener);
        final RealVector b = new ArrayRealVector(n);
        for (int j = 0; j < n; j++) {
            b.set(0.);
            b.setEntry(j, 1.);

            boolean caught = false;
            try {
                solver.solve(a, b);
            } catch (MathIllegalStateException e) {
                caught = true;
                final RealVector y = a.operate(x);
                for (int i = 0; i < n; i++) {
                    final double actual = b.getEntry(i) - y.getEntry(i);
                    final double expected = r.getEntry(i);
                    final double delta = 1E-6 * FastMath.abs(expected);
                    final String msg = String
                        .format("column %d, residual %d", i, j);
                    Assertions.assertEquals(expected, actual, delta, msg);
                }
            }
            Assertions
                .assertTrue(caught,
                            "MathIllegalStateException should have been caught");
        }
    }

    @Test
    public void testNonSquarePreconditioner() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final Array2DRowRealMatrix a = new Array2DRowRealMatrix(2, 2);
            final RealLinearOperator m = new RealLinearOperator() {

                @Override
                public RealVector operate(final RealVector x) {
                    throw new UnsupportedOperationException();
                }

                @Override
                public int getRowDimension() {
                    return 2;
                }

                @Override
                public int getColumnDimension() {
                    return 3;
                }
            };
            final PreconditionedIterativeLinearSolver solver;
            solver = new ConjugateGradient(10, 0d, false);
            final ArrayRealVector b = new ArrayRealVector(a.getRowDimension());
            solver.solve(a, m, b);
        });
    }

    @Test
    public void testMismatchedOperatorDimensions() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final Array2DRowRealMatrix a = new Array2DRowRealMatrix(2, 2);
            final RealLinearOperator m = new RealLinearOperator() {

                @Override
                public RealVector operate(final RealVector x) {
                    throw new UnsupportedOperationException();
                }

                @Override
                public int getRowDimension() {
                    return 3;
                }

                @Override
                public int getColumnDimension() {
                    return 3;
                }
            };
            final PreconditionedIterativeLinearSolver solver;
            solver = new ConjugateGradient(10, 0d, false);
            final ArrayRealVector b = new ArrayRealVector(a.getRowDimension());
            solver.solve(a, m, b);
        });
    }

    @Test
    public void testNonPositiveDefinitePreconditioner() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final Array2DRowRealMatrix a = new Array2DRowRealMatrix(2, 2);
            a.setEntry(0, 0, 1d);
            a.setEntry(0, 1, 2d);
            a.setEntry(1, 0, 3d);
            a.setEntry(1, 1, 4d);
            final RealLinearOperator m = new RealLinearOperator() {

                @Override
                public RealVector operate(final RealVector x) {
                    final ArrayRealVector y = new ArrayRealVector(2);
                    y.setEntry(0, -x.getEntry(0));
                    y.setEntry(1, x.getEntry(1));
                    return y;
                }

                @Override
                public int getRowDimension() {
                    return 2;
                }

                @Override
                public int getColumnDimension() {
                    return 2;
                }
            };
            final PreconditionedIterativeLinearSolver solver;
            solver = new ConjugateGradient(10, 0d, true);
            final ArrayRealVector b = new ArrayRealVector(2);
            b.setEntry(0, -1d);
            b.setEntry(1, -1d);
            solver.solve(a, m, b);
        });
    }

    @Test
    public void testPreconditionedSolution() {
        final int n = 8;
        final int maxIterations = 100;
        final RealLinearOperator a = new HilbertMatrix(n);
        final InverseHilbertMatrix ainv = new InverseHilbertMatrix(n);
        final RealLinearOperator m = JacobiPreconditioner.create(a);
        final PreconditionedIterativeLinearSolver solver;
        solver = new ConjugateGradient(maxIterations, 1E-15, true);
        final RealVector b = new ArrayRealVector(n);
        for (int j = 0; j < n; j++) {
            b.set(0.);
            b.setEntry(j, 1.);
            final RealVector x = solver.solve(a, m, b);
            for (int i = 0; i < n; i++) {
                final double actual = x.getEntry(i);
                final double expected = ainv.getEntry(i, j);
                final double delta = 1E-6 * FastMath.abs(expected);
                final String msg = String.format("coefficient (%d, %d)", i, j);
                Assertions.assertEquals(expected, actual, delta, msg);
            }
        }
    }

    @Test
    public void testPreconditionedResidual() {
        final int n = 10;
        final int maxIterations = n;
        final RealLinearOperator a = new HilbertMatrix(n);
        final RealLinearOperator m = JacobiPreconditioner.create(a);
        final ConjugateGradient solver;
        solver = new ConjugateGradient(maxIterations, 1E-15, true);
        final RealVector r = new ArrayRealVector(n);
        final RealVector x = new ArrayRealVector(n);
        final IterationListener listener = new IterationListener() {

            public void terminationPerformed(final IterationEvent e) {
                // Do nothing
            }

            public void iterationStarted(final IterationEvent e) {
                // Do nothing
            }

            public void iterationPerformed(final IterationEvent e) {
                final IterativeLinearSolverEvent evt;
                evt = (IterativeLinearSolverEvent) e;
                RealVector v = evt.getResidual();
                r.setSubVector(0, v);
                v = evt.getSolution();
                x.setSubVector(0, v);
            }

            public void initializationPerformed(final IterationEvent e) {
                // Do nothing
            }
        };
        solver.getIterationManager().addIterationListener(listener);
        final RealVector b = new ArrayRealVector(n);

        for (int j = 0; j < n; j++) {
            b.set(0.);
            b.setEntry(j, 1.);

            boolean caught = false;
            try {
                solver.solve(a, m, b);
            } catch (MathIllegalStateException e) {
                caught = true;
                final RealVector y = a.operate(x);
                for (int i = 0; i < n; i++) {
                    final double actual = b.getEntry(i) - y.getEntry(i);
                    final double expected = r.getEntry(i);
                    final double delta = 1E-6 * FastMath.abs(expected);
                    final String msg = String.format("column %d, residual %d", i, j);
                    Assertions.assertEquals(expected, actual, delta, msg);
                }
            }
            Assertions.assertTrue(caught, "MathIllegalStateException should have been caught");
        }
    }

    @Test
    public void testPreconditionedSolution2() {
        final int n = 100;
        final int maxIterations = 100000;
        final Array2DRowRealMatrix a = new Array2DRowRealMatrix(n, n);
        double daux = 1.;
        for (int i = 0; i < n; i++) {
            a.setEntry(i, i, daux);
            daux *= 1.2;
            for (int j = i + 1; j < n; j++) {
                if (i == j) {
                } else {
                    final double value = 1.0;
                    a.setEntry(i, j, value);
                    a.setEntry(j, i, value);
                }
            }
        }
        final RealLinearOperator m = JacobiPreconditioner.create(a);
        final PreconditionedIterativeLinearSolver pcg;
        final IterativeLinearSolver cg;
        pcg = new ConjugateGradient(maxIterations, 1E-6, true);
        cg = new ConjugateGradient(maxIterations, 1E-6, true);
        final RealVector b = new ArrayRealVector(n);
        final String pattern = "preconditioned gradient (%d iterations) should"
                               + " have been faster than unpreconditioned (%d iterations)";
        String msg;
        for (int j = 0; j < 1; j++) {
            b.set(0.);
            b.setEntry(j, 1.);
            final RealVector px = pcg.solve(a, m, b);
            final RealVector x = cg.solve(a, b);
            final int npcg = pcg.getIterationManager().getIterations();
            final int ncg = cg.getIterationManager().getIterations();
            msg = String.format(pattern, npcg, ncg);
            Assertions.assertTrue(npcg < ncg, msg);
            for (int i = 0; i < n; i++) {
                msg = String.format("row %d, column %d", i, j);
                final double expected = x.getEntry(i);
                final double actual = px.getEntry(i);
                final double delta = 1E-6 * FastMath.abs(expected);
                Assertions.assertEquals(expected, actual, delta, msg);
            }
        }
    }

    @Test
    public void testEventManagement() {
        final int n = 5;
        final int maxIterations = 100;
        final RealLinearOperator a = new HilbertMatrix(n);
        final IterativeLinearSolver solver;
        /*
         * count[0] = number of calls to initializationPerformed
         * count[1] = number of calls to iterationStarted
         * count[2] = number of calls to iterationPerformed
         * count[3] = number of calls to terminationPerformed
         */
        final int[] count = new int[] {0, 0, 0, 0};
        final IterationListener listener = new IterationListener() {
            private void doTestVectorsAreUnmodifiable(final IterationEvent e) {
                final IterativeLinearSolverEvent evt;
                evt = (IterativeLinearSolverEvent) e;
                try {
                    evt.getResidual().set(0.0);
                    Assertions.fail("r is modifiable");
                } catch (MathRuntimeException exc){
                    // Expected behavior
                }
                try {
                    evt.getRightHandSideVector().set(0.0);
                    Assertions.fail("b is modifiable");
                } catch (MathRuntimeException exc){
                    // Expected behavior
                }
                try {
                    evt.getSolution().set(0.0);
                    Assertions.fail("x is modifiable");
                } catch (MathRuntimeException exc){
                    // Expected behavior
                }
            }

            public void initializationPerformed(final IterationEvent e) {
                ++count[0];
                doTestVectorsAreUnmodifiable(e);
            }

            public void iterationPerformed(final IterationEvent e) {
                ++count[2];
                Assertions.assertEquals(count[2], e.getIterations() - 1, "iteration performed");
                doTestVectorsAreUnmodifiable(e);
            }

            public void iterationStarted(final IterationEvent e) {
                ++count[1];
                Assertions.assertEquals(count[1], e.getIterations() - 1, "iteration started");
                doTestVectorsAreUnmodifiable(e);
            }

            public void terminationPerformed(final IterationEvent e) {
                ++count[3];
                doTestVectorsAreUnmodifiable(e);
            }
        };
        solver = new ConjugateGradient(maxIterations, 1E-10, true);
        solver.getIterationManager().addIterationListener(listener);
        final RealVector b = new ArrayRealVector(n);
        for (int j = 0; j < n; j++) {
            Arrays.fill(count, 0);
            b.set(0.);
            b.setEntry(j, 1.);
            solver.solve(a, b);
            String msg = String.format("column %d (initialization)", j);
            Assertions.assertEquals(1, count[0], msg);
            msg = String.format("column %d (finalization)", j);
            Assertions.assertEquals(1, count[3], msg);
        }
    }

    @Test
    public void testUnpreconditionedNormOfResidual() {
        final int n = 5;
        final int maxIterations = 100;
        final RealLinearOperator a = new HilbertMatrix(n);
        final IterativeLinearSolver solver;
        final IterationListener listener = new IterationListener() {

            private void doTestNormOfResidual(final IterationEvent e) {
                final IterativeLinearSolverEvent evt;
                evt = (IterativeLinearSolverEvent) e;
                final RealVector x = evt.getSolution();
                final RealVector b = evt.getRightHandSideVector();
                final RealVector r = b.subtract(a.operate(x));
                final double rnorm = r.getNorm();
                Assertions.assertEquals(rnorm, evt.getNormOfResidual(),
                    FastMath.max(1E-5 * rnorm, 1E-10),
                    "iteration performed (residual)");
            }

            public void initializationPerformed(final IterationEvent e) {
                doTestNormOfResidual(e);
            }

            public void iterationPerformed(final IterationEvent e) {
                doTestNormOfResidual(e);
            }

            public void iterationStarted(final IterationEvent e) {
                doTestNormOfResidual(e);
            }

            public void terminationPerformed(final IterationEvent e) {
                doTestNormOfResidual(e);
            }
        };
        solver = new ConjugateGradient(maxIterations, 1E-10, true);
        solver.getIterationManager().addIterationListener(listener);
        final RealVector b = new ArrayRealVector(n);
        for (int j = 0; j < n; j++) {
            b.set(0.);
            b.setEntry(j, 1.);
            solver.solve(a, b);
        }
    }

    @Test
    public void testPreconditionedNormOfResidual() {
        final int n = 5;
        final int maxIterations = 100;
        final RealLinearOperator a = new HilbertMatrix(n);
        final RealLinearOperator m = JacobiPreconditioner.create(a);
        final PreconditionedIterativeLinearSolver solver;
        final IterationListener listener = new IterationListener() {

            private void doTestNormOfResidual(final IterationEvent e) {
                final IterativeLinearSolverEvent evt;
                evt = (IterativeLinearSolverEvent) e;
                final RealVector x = evt.getSolution();
                final RealVector b = evt.getRightHandSideVector();
                final RealVector r = b.subtract(a.operate(x));
                final double rnorm = r.getNorm();
                Assertions.assertEquals(rnorm, evt.getNormOfResidual(),
                    FastMath.max(1E-5 * rnorm, 1E-10),
                    "iteration performed (residual)");
            }

            public void initializationPerformed(final IterationEvent e) {
                doTestNormOfResidual(e);
            }

            public void iterationPerformed(final IterationEvent e) {
                doTestNormOfResidual(e);
            }

            public void iterationStarted(final IterationEvent e) {
                doTestNormOfResidual(e);
            }

            public void terminationPerformed(final IterationEvent e) {
                doTestNormOfResidual(e);
            }
        };
        solver = new ConjugateGradient(maxIterations, 1E-10, true);
        solver.getIterationManager().addIterationListener(listener);
        final RealVector b = new ArrayRealVector(n);
        for (int j = 0; j < n; j++) {
            b.set(0.);
            b.setEntry(j, 1.);
            solver.solve(a, m, b);
        }
    }
}
