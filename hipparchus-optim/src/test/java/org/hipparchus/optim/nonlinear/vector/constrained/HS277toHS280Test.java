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

import org.hipparchus.linear.*;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.SimpleBounds;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

/** HS277–HS280: Hilbert-system linear objective with Ax=b and x ≥ 0.
 *  N = 4 (HS277), 6 (HS278), 8 (HS279), 10 (HS280).
 */
public class HS277toHS280Test {

    /** Build the n×n Hilbert matrix H with H[i,j] = 1/(i+j+1) in 0-based indexing. */
    private static RealMatrix hilbert(final int n) {
        double[][] a = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                a[i][j] = 1.0 / (i + j + 1.0 + 0.0); // i+j+1 in 1-based -> +1 in 0-based -> here (i+j+1)
            }
        }
        return new Array2DRowRealMatrix(a, false);
    }

    /** Row-sum vector c_i = sum_j H_{ij}. */
    private static RealVector rowSums(final RealMatrix H) {
        int n = H.getRowDimension();
        double[] c = new double[n];
        for (int i = 0; i < n; i++) {
            double s = 0.0;
            for (int j = 0; j < n; j++) {
                s += H.getEntry(i, j);
            }
            c[i] = s;
        }
        return new ArrayRealVector(c, false);
    }

    /** Objective f(x) = cᵀx with constant gradient c. */
    static final class HS277Objective extends TwiceDifferentiableFunction {
        private final RealVector c;
        HS277Objective(RealVector c){ this.c = c; }
        @Override public int dim() { return c.getDimension(); }
        @Override public double value(RealVector x) { return c.dotProduct(x); }
        @Override public RealVector gradient(RealVector x) { return c.copy(); }
        @Override public RealMatrix hessian(RealVector x) { return new Array2DRowRealMatrix(c.getDimension(), c.getDimension()); }
    }

    /** Equality constraints g_i(x) = sum_j (x_j - 1)/(i + j - 1) = 0. */
static final class HS277Eq extends EqualityConstraint {
    private final RealMatrix H;   // teniamo i campi per non toccare la firma
    private final RealVector b;   // anche se non li usiamo più direttamente
    private final int n;

    HS277Eq(RealMatrix H, RealVector b) {
        super(new ArrayRealVector(new double[H.getRowDimension()])); // placeholder RHS
        this.H  = H;
        this.b  = b;
        this.n  = H.getColumnDimension();
    }

    @Override
    public int dim() {
        return n;
    }

    @Override
    public RealVector value(RealVector x) {
        // Fortran:
        // DO I = 1,N
        //   H = 0
        //   DO J = 1,N
        //     H = H + (X(J) - 1.D0)/DBLE(I+J-1)
        //   G(I) = H
        // END DO
        final int m = H.getRowDimension(); // = n per Hilbert, ma lo ricaviamo da H
        double[] g = new double[m];
        for (int i = 0; i < m; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                // i,j in Fortran vanno da 1..N, qui 0..N-1 → denominatore (i+j+1)
                sum += (x.getEntry(j) - 1.0) / (i + j + 1.0);
            }
            g[i] = sum;
        }
        return new ArrayRealVector(g, false);
    }

    @Override
    public RealMatrix jacobian(RealVector x) {
        // Derivata di g_i rispetto a x_j:
        // ∂/∂x_j [ sum_k (x_k - 1)/(i+k-1) ] = 1/(i+j-1)
        final int m = H.getRowDimension();
        double[][] J = new double[m][n];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                J[i][j] = 1.0 / (i + j + 1.0); // Fortran i+j-1 → Java +1
            }
        }
        return new Array2DRowRealMatrix(J, false);
    }
}


    private void runCase(final int n) {
        final RealMatrix H = hilbert(n);
        final RealVector c = rowSums(H);
        final RealVector b = c.copy(); // b_i = sum_j H_{ij}
        final double fEx = c.mapMultiply(1.0).dotProduct(new ArrayRealVector(n, 1.0)); // sum_i c_i

        final SQPOptimizerS2 optimizer = new SQPOptimizerS2();
        optimizer.setDebugPrinter(System.out::println);
        SQPOption sqpOption=new SQPOption();
        sqpOption.setGradientMode(GradientMode.FORWARD);

        final double[] start = new double[n]; // all zeros as in the Fortran setup
        final double[] lo = new double[n];
        final double[] up = new double[n];
        for (int i = 0; i < n; i++) {
//            start[i]=0.1;
            lo[i] = 0.0;
            up[i] = Double.POSITIVE_INFINITY;
        }

        LagrangeSolution sol = optimizer.optimize(
                sqpOption,
            //new InitialGuess(start),
            new ObjectiveFunction(new HS277Objective(c)),
            new HS277Eq(H, b),
           new SimpleBounds(lo, up)
        );

        double f = sol.getValue();
        assertEquals(fEx, f, 1.0e-6 * (Math.abs(fEx) + 1.0), "objective mismatch");
    }

    @Test public void testHS277() { runCase(4); }
    @Test public void testHS278() { runCase(6); }
//   @Test public void testHS279() { runCase(8); }
//   @Test public void testHS280() { runCase(10); }
}
