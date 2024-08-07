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

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.IterationEvent;
import org.hipparchus.util.IterationListener;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class SymmLQTest {

    public void saundersTest(final int n, final boolean goodb,
                             final boolean precon, final double shift,
                             final double pertbn) {
        final RealLinearOperator a = new RealLinearOperator() {

            @Override
            public RealVector operate(final RealVector x) {
                if (x.getDimension() != n) {
                    throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                           x.getDimension(), n);
                }
                final double[] y = new double[n];
                for (int i = 0; i < n; i++) {
                    y[i] = (i + 1) * 1.1 / n * x.getEntry(i);
                }
                return new ArrayRealVector(y, false);
            }

            @Override
            public int getRowDimension() {
                return n;
            }

            @Override
            public int getColumnDimension() {
                return n;
            }
        };
        final double shiftm = shift;
        final double pertm = FastMath.abs(pertbn);
        final RealLinearOperator minv;
        if (precon) {
            minv = new RealLinearOperator() {
                @Override
                public int getRowDimension() {
                    return n;
                }

                @Override
                public int getColumnDimension() {
                    return n;
                }

                @Override
                public RealVector operate(final RealVector x) {
                    if (x.getDimension() != n) {
                        throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                               x.getDimension(), n);
                    }
                    final double[] y = new double[n];
                    for (int i = 0; i < n; i++) {
                        double d = (i + 1) * 1.1 / n;
                        d = FastMath.abs(d - shiftm);
                        if (i % 10 == 0) {
                            d += pertm;
                        }
                        y[i] = x.getEntry(i) / d;
                    }
                    return new ArrayRealVector(y, false);
                }
            };
        } else {
            minv = null;
        }
        final RealVector xtrue = new ArrayRealVector(n);
        for (int i = 0; i < n; i++) {
            xtrue.setEntry(i, n - i);
        }
        final RealVector b = a.operate(xtrue);
        b.combineToSelf(1.0, -shift, xtrue);
        final SymmLQ solver = new SymmLQ(2 * n, 1E-12, true);
        final RealVector x = solver.solve(a, minv, b, goodb, shift);
        final RealVector y = a.operate(x);
        final RealVector r1 = new ArrayRealVector(n);
        for (int i = 0; i < n; i++) {
            final double bi = b.getEntry(i);
            final double yi = y.getEntry(i);
            final double xi = x.getEntry(i);
            r1.setEntry(i, bi - yi + shift * xi);
        }
        final double enorm = x.subtract(xtrue).getNorm() / xtrue.getNorm();
        final double etol = 1E-5;
        assertTrue(enorm <= etol, "enorm=" + enorm + ", " +
        solver.getIterationManager().getIterations());
    }

    @Test
    void testSolveSaunders1() {
        saundersTest(1, false, false, 0., 0.);
    }

    @Test
    void testSolveSaunders2() {
        saundersTest(2, false, false, 0., 0.);
    }

    @Test
    void testSolveSaunders3() {
        saundersTest(1, false, true, 0., 0.);
    }

    @Test
    void testSolveSaunders4() {
        saundersTest(2, false, true, 0., 0.);
    }

    @Test
    void testSolveSaunders5() {
        saundersTest(5, false, true, 0., 0.);
    }

    @Test
    void testSolveSaunders6() {
        saundersTest(5, false, true, 0.25, 0.);
    }

    @Test
    void testSolveSaunders7() {
        saundersTest(50, false, false, 0., 0.);
    }

    @Test
    void testSolveSaunders8() {
        saundersTest(50, false, false, 0.25, 0.);
    }

    @Test
    void testSolveSaunders9() {
        saundersTest(50, false, true, 0., 0.10);
    }

    @Test
    void testSolveSaunders10() {
        saundersTest(50, false, true, 0.25, 0.10);
    }

    @Test
    void testSolveSaunders11() {
        saundersTest(1, true, false, 0., 0.);
    }

    @Test
    void testSolveSaunders12() {
        saundersTest(2, true, false, 0., 0.);
    }

    @Test
    void testSolveSaunders13() {
        saundersTest(1, true, true, 0., 0.);
    }

    @Test
    void testSolveSaunders14() {
        saundersTest(2, true, true, 0., 0.);
    }

    @Test
    void testSolveSaunders15() {
        saundersTest(5, true, true, 0., 0.);
    }

    @Test
    void testSolveSaunders16() {
        saundersTest(5, true, true, 0.25, 0.);
    }

    @Test
    void testSolveSaunders17() {
        saundersTest(50, true, false, 0., 0.);
    }

    @Test
    void testSolveSaunders18() {
        saundersTest(50, true, false, 0.25, 0.);
    }

    @Test
    void testSolveSaunders19() {
        saundersTest(50, true, true, 0., 0.10);
    }

    @Test
    void testSolveSaunders20() {
        saundersTest(50, true, true, 0.25, 0.10);
    }

    @Test
    void testNonSquareOperator() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final Array2DRowRealMatrix a = new Array2DRowRealMatrix(2, 3);
            final IterativeLinearSolver solver;
            solver = new SymmLQ(10, 0., false);
            final ArrayRealVector b = new ArrayRealVector(a.getRowDimension());
            final ArrayRealVector x = new ArrayRealVector(a.getColumnDimension());
            solver.solve(a, b, x);
        });
    }

    @Test
    void testDimensionMismatchRightHandSide() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final Array2DRowRealMatrix a = new Array2DRowRealMatrix(3, 3);
            final IterativeLinearSolver solver;
            solver = new SymmLQ(10, 0., false);
            final ArrayRealVector b = new ArrayRealVector(2);
            solver.solve(a, b);
        });
    }

    @Test
    void testDimensionMismatchSolution() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final Array2DRowRealMatrix a = new Array2DRowRealMatrix(3, 3);
            final IterativeLinearSolver solver;
            solver = new SymmLQ(10, 0., false);
            final ArrayRealVector b = new ArrayRealVector(3);
            final ArrayRealVector x = new ArrayRealVector(2);
            solver.solve(a, b, x);
        });
    }

    @Test
    void testUnpreconditionedSolution() {
        final int n = 5;
        final int maxIterations = 100;
        final RealLinearOperator a = new HilbertMatrix(n);
        final InverseHilbertMatrix ainv = new InverseHilbertMatrix(n);
        final IterativeLinearSolver solver;
        solver = new SymmLQ(maxIterations, 1E-10, true);
        final RealVector b = new ArrayRealVector(n);
        for (int j = 0; j < n; j++) {
            b.set(0.);
            b.setEntry(j, 1.);
            final RealVector x = solver.solve(a, b);
            for (int i = 0; i < n; i++) {
                final double actual = x.getEntry(i);
                final double expected = ainv.getEntry(i, j);
                final double delta = 1E-6 * FastMath.abs(expected);
                final String msg = String.format("entry[%d][%d]", i, j);
                assertEquals(expected, actual, delta, msg);
            }
        }
    }

    @Test
    void testUnpreconditionedInPlaceSolutionWithInitialGuess() {
        final int n = 5;
        final int maxIterations = 100;
        final RealLinearOperator a = new HilbertMatrix(n);
        final InverseHilbertMatrix ainv = new InverseHilbertMatrix(n);
        final IterativeLinearSolver solver;
        solver = new SymmLQ(maxIterations, 1E-10, true);
        final RealVector b = new ArrayRealVector(n);
        for (int j = 0; j < n; j++) {
            b.set(0.);
            b.setEntry(j, 1.);
            final RealVector x0 = new ArrayRealVector(n);
            x0.set(1.);
            final RealVector x = solver.solveInPlace(a, b, x0);
            assertSame(x0, x, "x should be a reference to x0");
            for (int i = 0; i < n; i++) {
                final double actual = x.getEntry(i);
                final double expected = ainv.getEntry(i, j);
                final double delta = 1E-6 * FastMath.abs(expected);
                final String msg = String.format("entry[%d][%d)", i, j);
                assertEquals(expected, actual, delta, msg);
            }
        }
    }

    @Test
    void testUnpreconditionedSolutionWithInitialGuess() {
        final int n = 5;
        final int maxIterations = 100;
        final RealLinearOperator a = new HilbertMatrix(n);
        final InverseHilbertMatrix ainv = new InverseHilbertMatrix(n);
        final IterativeLinearSolver solver;
        solver = new SymmLQ(maxIterations, 1E-10, true);
        final RealVector b = new ArrayRealVector(n);
        for (int j = 0; j < n; j++) {
            b.set(0.);
            b.setEntry(j, 1.);
            final RealVector x0 = new ArrayRealVector(n);
            x0.set(1.);
            final RealVector x = solver.solve(a, b, x0);
            assertNotSame(x0, x, "x should not be a reference to x0");
            for (int i = 0; i < n; i++) {
                final double actual = x.getEntry(i);
                final double expected = ainv.getEntry(i, j);
                final double delta = 1E-6 * FastMath.abs(expected);
                final String msg = String.format("entry[%d][%d]", i, j);
                assertEquals(expected, actual, delta, msg);
                assertEquals(1., x0.getEntry(i), Math.ulp(1.), msg);
            }
        }
    }

    @Test
    void testNonSquarePreconditioner() {
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
            solver = new SymmLQ(10, 0., false);
            final ArrayRealVector b = new ArrayRealVector(a.getRowDimension());
            solver.solve(a, m, b);
        });
    }

    @Test
    void testMismatchedOperatorDimensions() {
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
            solver = new SymmLQ(10, 0d, false);
            final ArrayRealVector b = new ArrayRealVector(a.getRowDimension());
            solver.solve(a, m, b);
        });
    }

    @Test
    void testNonPositiveDefinitePreconditioner() {
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
                    y.setEntry(1, -x.getEntry(1));
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
            solver = new SymmLQ(10, 0d, true);
            final ArrayRealVector b = new ArrayRealVector(2);
            b.setEntry(0, -1d);
            b.setEntry(1, -1d);
            solver.solve(a, m, b);
        });
    }

    @Test
    void testPreconditionedSolution() {
        final int n = 8;
        final int maxIterations = 100;
        final RealLinearOperator a = new HilbertMatrix(n);
        final InverseHilbertMatrix ainv = new InverseHilbertMatrix(n);
        final RealLinearOperator m = JacobiPreconditioner.create(a);
        final PreconditionedIterativeLinearSolver solver;
        solver = new SymmLQ(maxIterations, 1E-15, true);
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
                assertEquals(expected, actual, delta, msg);
            }
        }
    }

    @Test
    void testPreconditionedSolution2() {
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
        final PreconditionedIterativeLinearSolver prec;
        final IterativeLinearSolver unprec;
        prec = new SymmLQ(maxIterations, 1E-15, true);
        unprec = new SymmLQ(maxIterations, 1E-15, true);
        final RealVector b = new ArrayRealVector(n);
        final String pattern = "preconditioned SymmLQ (%d iterations) should"
                               + " have been faster than unpreconditioned (%d iterations)";
        String msg;
        for (int j = 0; j < 1; j++) {
            b.set(0.);
            b.setEntry(j, 1.);
            final RealVector px = prec.solve(a, m, b);
            final RealVector x = unprec.solve(a, b);
            final int np = prec.getIterationManager().getIterations();
            final int nup = unprec.getIterationManager().getIterations();
            msg = String.format(pattern, np, nup);
            for (int i = 0; i < n; i++) {
                msg = String.format("row %d, column %d", i, j);
                final double expected = x.getEntry(i);
                final double actual = px.getEntry(i);
                final double delta = 5E-5 * FastMath.abs(expected);
                assertEquals(expected, actual, delta, msg);
            }
        }
    }

    @Test
    void testEventManagement() {
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
        final RealVector xFromListener = new ArrayRealVector(n);
        final IterationListener listener = new IterationListener() {

            public void initializationPerformed(final IterationEvent e) {
                ++count[0];
            }

            public void iterationPerformed(final IterationEvent e) {
                ++count[2];
                assertEquals(count[2],
                                    e.getIterations() - 1,
                                    "iteration performed");
            }

            public void iterationStarted(final IterationEvent e) {
                ++count[1];
                assertEquals(count[1],
                                    e.getIterations() - 1,
                                    "iteration started");
            }

            public void terminationPerformed(final IterationEvent e) {
                ++count[3];
                final IterativeLinearSolverEvent ilse;
                ilse = (IterativeLinearSolverEvent) e;
                xFromListener.setSubVector(0, ilse.getSolution());
            }
        };
        solver = new SymmLQ(maxIterations, 1E-10, true);
        solver.getIterationManager().addIterationListener(listener);
        final RealVector b = new ArrayRealVector(n);
        for (int j = 0; j < n; j++) {
            Arrays.fill(count, 0);
            b.set(0.);
            b.setEntry(j, 1.);
            final RealVector xFromSolver = solver.solve(a, b);
            String msg = String.format("column %d (initialization)", j);
            assertEquals(1, count[0], msg);
            msg = String.format("column %d (finalization)", j);
            assertEquals(1, count[3], msg);
            /*
             *  Check that solution is not "over-refined". When the last
             *  iteration has occurred, no further refinement should be
             *  performed.
             */
            for (int i = 0; i < n; i++){
                msg = String.format("row %d, column %d", i, j);
                final double expected = xFromSolver.getEntry(i);
                final double actual = xFromListener.getEntry(i);
                assertEquals(expected, actual, 0.0, msg);
            }
        }
    }

    @Test
    void testNonSelfAdjointOperator() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final RealLinearOperator a;
            a = new Array2DRowRealMatrix(new double[][]{
                {1., 2., 3.},
                {2., 4., 5.},
                {2.999, 5., 6.}
            });
            final RealVector b;
            b = new ArrayRealVector(new double[]{1., 1., 1.});
            new SymmLQ(100, 1., true).solve(a, b);
        });
    }

    @Test
    void testNonSelfAdjointPreconditioner() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final RealLinearOperator a = new Array2DRowRealMatrix(new double[][]{
                {1., 2., 3.},
                {2., 4., 5.},
                {3., 5., 6.}
            });
            final Array2DRowRealMatrix mMat;
            mMat = new Array2DRowRealMatrix(new double[][]{
                {1., 0., 1.},
                {0., 1., 0.},
                {0., 0., 1.}
            });
            final DecompositionSolver mSolver;
            mSolver = new LUDecomposition(mMat).getSolver();
            final RealLinearOperator minv = new RealLinearOperator() {

                @Override
                public RealVector operate(final RealVector x) {
                    return mSolver.solve(x);
                }

                @Override
                public int getRowDimension() {
                    return mMat.getRowDimension();
                }

                @Override
                public int getColumnDimension() {
                    return mMat.getColumnDimension();
                }
            };
            final RealVector b = new ArrayRealVector(new double[]{
                1., 1., 1.
            });
            new SymmLQ(100, 1., true).solve(a, minv, b);
        });
    }

    @Test
    void testUnpreconditionedNormOfResidual() {
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
                assertEquals(rnorm, evt.getNormOfResidual(),
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
        solver = new SymmLQ(maxIterations, 1E-10, true);
        solver.getIterationManager().addIterationListener(listener);
        final RealVector b = new ArrayRealVector(n);
        for (int j = 0; j < n; j++) {
            b.set(0.);
            b.setEntry(j, 1.);
            solver.solve(a, b);
        }
    }

    @Test
    void testPreconditionedNormOfResidual() {
        final int n = 5;
        final int maxIterations = 100;
        final RealLinearOperator a = new HilbertMatrix(n);
        final JacobiPreconditioner m = JacobiPreconditioner.create(a);
        final RealLinearOperator p = m.sqrt();
        final PreconditionedIterativeLinearSolver solver;
        final IterationListener listener = new IterationListener() {

            private void doTestNormOfResidual(final IterationEvent e) {
                final IterativeLinearSolverEvent evt;
                evt = (IterativeLinearSolverEvent) e;
                final RealVector x = evt.getSolution();
                final RealVector b = evt.getRightHandSideVector();
                final RealVector r = b.subtract(a.operate(x));
                final double rnorm = p.operate(r).getNorm();
                assertEquals(rnorm, evt.getNormOfResidual(),
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
        solver = new SymmLQ(maxIterations, 1E-10, true);
        solver.getIterationManager().addIterationListener(listener);
        final RealVector b = new ArrayRealVector(n);
        for (int j = 0; j < n; j++) {
            b.set(0.);
            b.setEntry(j, 1.);
            solver.solve(a, m, b);
        }
    }
}

