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
import org.hipparchus.util.ArithmeticUtils;
import org.hipparchus.util.CombinatoricsUtils;

/**
 * This class implements inverses of Hilbert Matrices as
 * {@link RealLinearOperator}.
 */
public class InverseHilbertMatrix
    implements RealLinearOperator {

    /** The size of the matrix. */
    private final int n;

    /**
     * Creates a new instance of this class.
     *
     * @param n Size of the matrix to be created.
     */
    public InverseHilbertMatrix(final int n) {
        this.n = n;
    }

    /** {@inheritDoc} */
    @Override
    public int getColumnDimension() {
        return n;
    }

    /**
     * Returns the {@code (i, j)} entry of the inverse Hilbert matrix. Exact
     * arithmetic is used; in case of overflow, an exception is thrown.
     *
     * @param i Row index (starts at 0).
     * @param j Column index (starts at 0).
     * @return The coefficient of the inverse Hilbert matrix.
     */
    public long getEntry(final int i, final int j) {
        long val = i + j + 1;
        long aux = CombinatoricsUtils.binomialCoefficient(n + i, n - j - 1);
        val = ArithmeticUtils.mulAndCheck(val, aux);
        aux = CombinatoricsUtils.binomialCoefficient(n + j, n - i - 1);
        val = ArithmeticUtils.mulAndCheck(val, aux);
        aux = CombinatoricsUtils.binomialCoefficient(i + j, i);
        val = ArithmeticUtils.mulAndCheck(val, aux);
        val = ArithmeticUtils.mulAndCheck(val, aux);
        return ((i + j) & 1) == 0 ? val : -val;
    }

    /** {@inheritDoc} */
    @Override
    public int getRowDimension() {
        return n;
    }

    /** {@inheritDoc} */
    @Override
    public RealVector operate(final RealVector x) {
        if (x.getDimension() != n) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   x.getDimension(), n);
        }
        final double[] y = new double[n];
        for (int i = 0; i < n; i++) {
            double pos = 0.;
            double neg = 0.;
            for (int j = 0; j < n; j++) {
                final double xj = x.getEntry(j);
                final long coeff = getEntry(i, j);
                final double daux = coeff * xj;
                // Positive and negative values are sorted out in order to limit
                // catastrophic cancellations (do not forget that Hilbert
                // matrices are *very* ill-conditioned!
                if (daux > 0.) {
                    pos += daux;
                } else {
                    neg += daux;
                }
            }
            y[i] = pos + neg;
        }
        return new ArrayRealVector(y, false);
    }
}
