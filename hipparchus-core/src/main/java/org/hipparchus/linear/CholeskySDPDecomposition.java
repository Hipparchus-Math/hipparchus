/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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

package org.hipparchus.linear;

import java.util.Arrays;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;

/**
 * Calculates the Cholesky decomposition of a positive semidefinite matrix.
 * <p>The classic Cholesky decomposition ({@link CholeskyDecomposition}) applies to real symmetric positive-definite
 * matrix. This class extends cholesky decomposition to positive semidefinite matrix. 
 * <p>This class is based on the {@code schol} Matlab function from
 * <a href="https://github.com/EEA-sensors/ekfukf">EKF-UKF Toolbox</a>
 * @see J. Hartikainen, A. Solin, and S. Särkkä “Optimal ﬁltering with Kalman ﬁlters and smoothers,” 
 * Dept. of Biomedica Engineering and Computational Sciences, Aalto University School of Science, Aug. 2011.
 * Online copy: <a href="https://www.researchgate.net/profile/Simo-Saerkkae/publication/228683456_Optimal_filtering_with_Kalman_filters_and_smoothers-a_Manual_for_Matlab_toolbox_EKFUKF/links/0deec52b885bd0d890000000/Optimal-filtering-with-Kalman-filters-and-smoothers-a-Manual-for-Matlab-toolbox-EKF-UKF.pdf?origin=publication_detail">Optimal Filtering</a>
 * @see {@code schol} function: <a href="https://github.com/EEA-sensors/ekfukf/blob/develop/schol.m">schol</a>
 */
public class CholeskySDPDecomposition {
    
    /** Row-oriented storage for L<sup>T</sup> matrix data. */
    private final double[][] lTData;
    /** Cached value of L. */
    private RealMatrix cachedL;
    /** Cached value of LT. */
    private RealMatrix cachedLT;
    /** Value 1,0,-1 denoting that A was positive definite,
        positive semidefinite or negative definite, respectively. */
    private int def;
    
    /**
     * Calculates the Cholesky decomposition of the given matrix.
     * @param matrix the matrix to decompose
     * @param matrix
     * @throws MathIllegalArgumentException if the matrix is not square.
     */
    public CholeskySDPDecomposition(final RealMatrix matrix) {
        if (!matrix.isSquare()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NON_SQUARE_MATRIX,
                                                   matrix.getRowDimension(), matrix.getColumnDimension());
        }
        
        final int order = matrix.getRowDimension();
        lTData   = matrix.getData();
        cachedL  = MatrixUtils.createRealMatrix(lTData);
        
        
        double zeroArray[] = new double[order];
        Arrays.fill(zeroArray, 0.);
        
        for (int i = 0; i < order; ++i) {
            cachedL.setColumn(i, zeroArray);
        }
        
        cachedLT = cachedL.transpose();
        
        final double e = 1.0;
        final double eps = Math.ulp(e);
        
        for (int i = 0; i < order; ++i) {
            for (int j = 0; j < i + 1; j++) {
                double s = lTData[i][j];
                for (int k = 0; k < j; ++k) {
                    s = s - cachedL.getEntry(i, k) * cachedL.getEntry(j, k);
                }
                if (j < i) {
                    if (cachedL.getEntry(j, j) > eps) {
                        cachedL.setEntry(i, j, s / cachedL.getEntry(j, j));
                    }
                    else {
                        cachedL.setEntry(i, j, 0.);
                    }
                }
                else {
                    if (s < -eps) {
                        s = 0;
                        def = -1;
                    }
                    else if (s < eps) {
                        s = 0;
                        def = Math.min(0, def);
                    }
                    cachedL.setEntry(j,j,Math.sqrt(s));
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
        
        if (cachedL == null) {
            cachedL = MatrixUtils.createRealMatrix(lTData);
        }
        return cachedL;
    }

    /**
     * Returns the transpose of the matrix L of the decomposition.
     * <p>L<sup>T</sup> is an upper-triangular matrix</p>
     * @return the transpose of the matrix L of the decomposition
     */
    public RealMatrix getLT() {
        cachedLT = cachedL.transpose();
        // return the cached matrix
        return cachedLT;
    }

    /**
     * Return the determinant of the matrix
     * @return determinant of the matrix
     */
    public double getDeterminant() {
        double determinant = 1.0;
        for (int i = 0; i < lTData.length; ++i) {
            double lTii = lTData[i][i];
            determinant *= lTii * lTii;
        }
        return determinant;
    }

}
