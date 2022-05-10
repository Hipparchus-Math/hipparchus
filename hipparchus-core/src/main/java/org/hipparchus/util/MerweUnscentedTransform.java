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

/** Perform the unscented transform as defined by Wan and Merwe.
 * @see E. A. Wan and R. Van der Merwe, “The unscented Kalman filter for nonlinear estimation,”
 * in Proc. Symp. Adaptive Syst. Signal Process., Commun. Contr., Lake Louise, AB, Canada, Oct. 2000.
 * Online copy: <a href = "https://www.seas.harvard.edu/courses/cs281/papers/unscented.pdf">https://www.seas.harvard.edu/courses/cs281/papers/unscented.pdf</a>
 * */

public class MerweUnscentedTransform implements UnscentedTransformProvider {

    /** Default value for alpha (1.0E-4), see reference. */
    public static final double DEFAULT_ALPHA = 0.0001;

    /** Default value for beta (2.), see reference. */
    public static final double DEFAULT_BETA = 2;

    /** Default value for kappa, (0.), see reference. */
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
     * @param alpha scaling control
     * @param beta free parameter
     * @param kappa free parameter
     */
    public MerweUnscentedTransform(final int stateDim, final double alpha, final double beta, final double kappa) {
        // Check state dimension
        if (stateDim == 0) {
            // State dimension must be different from 0
            throw new MathIllegalArgumentException(LocalizedCoreFormats.ZERO_STATE_SIZE);
        }

        // lambda = alpha² + (n + kappa) - n (see Eq. 15)
        final double lambda = alpha * alpha * (stateDim + kappa) - stateDim;

        // Initialize multiplication factor for covariance matrix
        this.factor = stateDim + lambda;

        // Initialize vectors weights
        wm = new ArrayRealVector(2 * stateDim + 1);
        wc = new ArrayRealVector(2 * stateDim + 1);

        // Computation of unscented kalman filter weights (See Eq. 15)
        wm.setEntry(0, lambda / (stateDim + lambda));
        wc.setEntry(0, lambda / (stateDim + lambda) + 1 - alpha * alpha + beta);
        final double w = 1 / (2 * (stateDim + lambda));
        for (int i = 1; i <= 2 * stateDim; i++) {
            wm.setEntry(i, w);
            wc.setEntry(i, w);
        }
    }

    /**
     * Default constructor.
     * <p>
     * This constructor uses default values for <code>alpha</code>,
     * <code>beta</code>, and <code>kappa</code>.
     * Default alpha is equal to {@link #DEFAULT_ALPHA}, default beta is equal to {@link #DEFAULT_BETA} and default kappa is equal to {@link #DEFAULT_KAPPA}.
     * </p>
     * @param stateDim state dimension
     */
    public MerweUnscentedTransform(final int stateDim) {
        this(stateDim, DEFAULT_ALPHA, DEFAULT_BETA, DEFAULT_KAPPA);
    }


    /** Compute sigma points through unscented transform (cf. Eq.15) from a state.
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

            // Compute Eq. 13
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
