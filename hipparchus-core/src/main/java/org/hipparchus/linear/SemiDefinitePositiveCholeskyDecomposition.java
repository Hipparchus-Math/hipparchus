/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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

package org.hipparchus.linear;

import java.util.Arrays;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;

/**
 * Calculates the Cholesky decomposition of a positive semidefinite matrix.
 * <p>
 * The classic Cholesky decomposition ({@link CholeskyDecomposition}) applies to real
 * symmetric positive-definite matrix. This class extends the Cholesky decomposition to
 * positive semidefinite matrix. The main application is for estimation based on the
 * Unscented Kalman Filter.
 *
 * @see "J. Hartikainen, A. Solin, and S. Särkkä. Optimal ﬁltering with Kalman ﬁlters and smoothers,
 *       Dept. of Biomedica Engineering and Computational Sciences, Aalto University School of Science, Aug. 2011."
 *
 * @since 2.2
 */
public class SemiDefinitePositiveCholeskyDecomposition {

    /** Default threshold below which elements are not considered positive. */
    public static final double POSITIVITY_THRESHOLD = 1.0e-15;
    /** Cached value of L. */
    private RealMatrix cachedL;
    /** Cached value of LT. */
    private RealMatrix cachedLT;

    /**
     * Calculates the Cholesky decomposition of the given matrix.
     * @param matrix the matrix to decompose
     * @throws MathIllegalArgumentException if the matrix is not square.
     * @see #SemiDefinitePositiveCholeskyDecomposition(RealMatrix, double)
     * @see #POSITIVITY_THRESHOLD
     */
    public SemiDefinitePositiveCholeskyDecomposition(final RealMatrix matrix) {
        this(matrix, POSITIVITY_THRESHOLD);
    }

    /**
     * Calculates the Cholesky decomposition of the given matrix.
     * @param matrix the matrix to decompose
     * @param positivityThreshold threshold below which elements are not considered positive
     * @throws MathIllegalArgumentException if the matrix is not square.
     */
    public SemiDefinitePositiveCholeskyDecomposition(final RealMatrix matrix,
                                                     final double positivityThreshold) {
        if (!matrix.isSquare()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NON_SQUARE_MATRIX,
                                                   matrix.getRowDimension(), matrix.getColumnDimension());
        }

        final int order = matrix.getRowDimension();
        final double[][] lTData = matrix.getData();
        cachedL  = MatrixUtils.createRealMatrix(lTData);
        int def  = 1;

        final double[] zeroArray = new double[order];
        Arrays.fill(zeroArray, 0.);

        for (int i = 0; i < order; ++i) {
            cachedL.setColumn(i, zeroArray);
        }

        cachedLT = cachedL.transpose();

        for (int i = 0; i < order; ++i) {
            for (int j = 0; j < i + 1; j++) {
                double s = lTData[i][j];
                for (int k = 0; k < j; ++k) {
                    s = s - cachedL.getEntry(i, k) * cachedL.getEntry(j, k);
                }
                if (j < i) {
                    if (cachedL.getEntry(j, j) > FastMath.ulp(1.0)) {
                        cachedL.setEntry(i, j, s / cachedL.getEntry(j, j));
                    } else {
                        cachedL.setEntry(i, j, 0.);
                    }
                } else {
                    if (s < -positivityThreshold) {
                        s = 0;
                        def = -1;
                    } else if (s < positivityThreshold) {
                        s = 0;
                        def = FastMath.min(0, def);
                    }
                    cachedL.setEntry(j, j, FastMath.sqrt(s));
                }

            }
        }

        if (def < 0) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NEGATIVE_DEFINITE_MATRIX);
        }

    }

    /**
     * Returns the matrix L of the decomposition.
     * <p>L is an lower-triangular matrix</p>
     * @return the L matrix
     */
    public RealMatrix getL() {
        // return the cached matrix
        return cachedL;
    }

    /**
     * Returns the transpose of the matrix L of the decomposition.
     * <p>L<sup>T</sup> is an upper-triangular matrix</p>
     * @return the transpose of the matrix L of the decomposition
     */
    public RealMatrix getLT() {
        cachedLT = getL().transpose();
        // return the cached matrix
        return cachedLT;
    }

}
