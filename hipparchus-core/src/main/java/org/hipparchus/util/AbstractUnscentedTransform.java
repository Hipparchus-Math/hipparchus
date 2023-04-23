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
package org.hipparchus.util;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.linear.SemiDefinitePositiveCholeskyDecomposition;

/**
 * Base class for unscented transform providers.
 * @since 2.2
 */
public abstract class AbstractUnscentedTransform implements UnscentedTransformProvider {

    /**
     * Constructor.
     * @param stateDim the dimension of the state
     */
    public AbstractUnscentedTransform(final int stateDim) {
        // Check state dimension
        if (stateDim == 0) {
            // State dimension must be different from 0
            throw new MathIllegalArgumentException(LocalizedCoreFormats.ZERO_STATE_SIZE);
        }
    }

    /**
     * {@inheritDoc}
     * <p>
     * Let n be the state dimension and Si be the ith row of the covariance matrix square root.
     * The returned array is organized as follow. Element 0 contains the process state, also
     * called the mean state. Elements from 1 to n contain the process state + Si. Finally,
     * elements from n + 1 to 2n contain the process state - Si
     */
    @Override
    public RealVector[] unscentedTransform(final RealVector state, final RealMatrix covariance) {

        // State dimensions
        final int n = state.getDimension();

        // Initialize array containing sigma points
        final RealVector[] sigmaPoints = new ArrayRealVector[(2 * n) + 1];
        sigmaPoints[0] = state;

        // Apply multiplication factor to the covariance matrix
        final double     factor = getMultiplicationFactor();
        final RealMatrix temp   = covariance.scalarMultiply(factor);

        // Compute lower triangular matrix of Cholesky decomposition
        // Note: When the estimation error covariance is propagated, it sometimes
        //       cannot maintain the positive semidefiniteness.
        //       To enhance the numerical stability of the unscented transform,
        //       the semidefinite positive Cholesky decomposition is used.
        final RealMatrix L = new SemiDefinitePositiveCholeskyDecomposition(temp).getL();

        // Compute sigma points
        for (int i = 1; i <= n; i++) {
            sigmaPoints[i]     = sigmaPoints[0].add(L.getColumnVector(i - 1));
            sigmaPoints[i + n] = sigmaPoints[0].subtract(L.getColumnVector(i - 1));
        }

        // Return sigma points
        return sigmaPoints;

    }

    /**
     * Get the factor applied to the covariance matrix during the unscented transform.
     * @return the factor applied to the covariance matrix during the unscented transform
     */
    protected abstract double getMultiplicationFactor();

}
