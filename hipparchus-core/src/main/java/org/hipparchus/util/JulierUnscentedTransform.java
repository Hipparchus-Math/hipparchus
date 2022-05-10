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
package org.hipparchus.util;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.CholeskyDecomposition;
import org.hipparchus.linear.CholeskySDPDecomposition;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;

/** Perform the unscented transform as defined by Julier and Uhlmann.
 * Particular case of {@link MerweUnscentedTransform} with alpha = 1 and beta = 0.
 * @see S. J. Julier and J. K. Uhlmann. “A New Extension of the Kalman Filter to Nonlinear Systems”.
 * Proc. SPIE 3068, Signal Processing, Sensor Fusion, and Target Recognition VI, 182 (July 28, 1997)
 * Online copy: <a href = "https://people.eecs.berkeley.edu/~pabbeel/cs287-fa19/optreadings/JulierUhlmann-UKF.pdf">https://people.eecs.berkeley.edu/~pabbeel/cs287-fa19/optreadings/JulierUhlmann-UKF.pdf</a>
 * */
public class JulierUnscentedTransform implements UnscentedTransformProvider {

    /** Default value for kappa (0.). */
    public static final double DEFAULT_KAPPA = 0;

    /** Weights for covariance matrix. */
    private final RealVector wc;

    /** Weights for mean state. */
    private final RealVector wm;

    /** Factor applied during Unscented transform for covariance matrix. */
    private final double factor;

    /**
     * Simple constructor.
     * @param stateDim the dimension of the state
     * @param kappa free parameter
     */
    public JulierUnscentedTransform(final int stateDim, final double kappa) {
        // Check state dimension
        if (stateDim == 0) {
            // State dimension must be different from 0
            throw new MathIllegalArgumentException(LocalizedCoreFormats.ZERO_STATE_SIZE);
        }



        // Initialize multiplication factor for covariance matrix
        this.factor = stateDim + kappa;

        // Initialize vectors weights
        wm = new ArrayRealVector(2 * stateDim + 1);

        // Computation of unscented kalman filter weights (See Eq. 12)
        wm.setEntry(0, kappa / (stateDim + kappa));
        for (int i = 1; i <= 2 * stateDim; i++) {
            wm.setEntry(i, 1 / (2 * (stateDim + kappa)));
        }
        wc = wm;

    }
    /**
     * Default constructor.
     * <p>
     * This constructor uses default values for <<code>kappa</code>.
     * Default kappa is equal to {@link #DEFAULT_KAPPA}.
     * </p>
     * @param stateDim state dimension
     */
    public JulierUnscentedTransform(final int stateDim) {
        this(stateDim, DEFAULT_KAPPA);

    }

    /** Compute sigma points through unscented transform (cf. Eq.12) from a state and covariance.
     * @param state state
     * @covariance covariance
     * @return sigma points.
     */
    @Override
    public RealVector[] unscentedTransform(final RealVector state, final RealMatrix covariance) {
        // State dimensions
        final int n = state.getDimension();
        // Initialize array containing sigma points
        final RealVector[] sigmaPoints = new ArrayRealVector[(2 * n) + 1];
        sigmaPoints[0] = state;
        final RealMatrix temp = covariance.scalarMultiply(factor);

        // Compute lower triangular matrix of Cholesky decomposition
        RealMatrix L;
        try {
            final CholeskyDecomposition chdecomposition = new CholeskyDecomposition(temp, 10e-14, 10e-14);
            L = chdecomposition.getL();
        }
        catch (MathIllegalArgumentException miae) {

            final CholeskySDPDecomposition csdpDecomposition = new CholeskySDPDecomposition(temp);
            L = csdpDecomposition.getL();

        }

        // Compute sigma points
        for (int i = 1; i <= n; i++) {

            // Compute Eq. 12
            sigmaPoints[i] = sigmaPoints[0].add(L.getColumnVector(i - 1));

            // Compute Eq. 12
            sigmaPoints[i + n] = sigmaPoints[0].subtract(L.getColumnVector(i - 1));

        }

        // Return sigma points
        return sigmaPoints;
    }

    /** {@inheritDoc} */
    @Override
    public RealVector getWc() {
        return wc;
    }

    /** {@inheritDoc} */
    @Override
    public RealVector getWm() {
        return wm;
    }





}
