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

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;

/**
 * Unscented transform as defined by Merwe and Wan.
 * <p>
 * The unscented transform uses three parameters: alpha, beta and kappa.
 * Alpha determines the spread of the sigma points around the process state,
 * kappa is a secondary scaling parameter, and beta is used to incorporate
 * prior knowledge of the distribution of the process state.
 * </p>
 * @see "E. A. Wan and R. Van der Merwe, The unscented Kalman filter for nonlinear estimation,
 *       in Proc. Symp. Adaptive Syst. Signal Process., Commun. Contr., Lake Louise, AB, Canada, Oct. 2000."
 * @since 2.2
 */
public class MerweUnscentedTransform extends AbstractUnscentedTransform {

    /** Default value for alpha (0.5, see reference). */
    public static final double DEFAULT_ALPHA = 0.5;

    /** Default value for beta (2.0, see reference). */
    public static final double DEFAULT_BETA = 2.0;

    /** Default value for kappa, (0.0, see reference). */
    public static final double DEFAULT_KAPPA = 0.0;

    /** Weights for covariance matrix. */
    private final RealVector wc;

    /** Weights for mean state. */
    private final RealVector wm;

    /** Factor applied to the covariance matrix during the unscented transform (lambda + process state size). */
    private final double factor;

    /**
     * Default constructor.
     * <p>
     * This constructor uses default values for alpha, beta, and kappa.
     * </p>
     * @param stateDim the dimension of the state
     * @see #DEFAULT_ALPHA
     * @see #DEFAULT_BETA
     * @see #DEFAULT_KAPPA
     * @see #MerweUnscentedTransform(int, double, double, double)
     */
    public MerweUnscentedTransform(final int stateDim) {
        this(stateDim, DEFAULT_ALPHA, DEFAULT_BETA, DEFAULT_KAPPA);
    }

    /**
     * Simple constructor.
     * @param stateDim the dimension of the state
     * @param alpha scaling control parameter
     *        (determines the spread of the sigma points around the process state)
     * @param beta free parameter
     *        (used to incorporate prior knowledge of the distribution of the process state)
     * @param kappa secondary scaling factor
     *        (usually set to 0.0)
     */
    public MerweUnscentedTransform(final int stateDim, final double alpha,
                                   final double beta, final double kappa) {

        // Call super constructor
        super(stateDim);

        // lambda = alpha² + (n + kappa) - n (see Eq. 15)
        final double lambda = alpha * alpha * (stateDim + kappa) - stateDim;

        // Initialize multiplication factor for covariance matrix
        this.factor = stateDim + lambda;

        // Initialize vectors weights
        wm = new ArrayRealVector(2 * stateDim + 1);
        wc = new ArrayRealVector(2 * stateDim + 1);

        // Computation of unscented kalman filter weights (See Eq. 15)
        wm.setEntry(0, lambda / factor);
        wc.setEntry(0, lambda / factor + (1.0 - alpha * alpha + beta));
        final double w = 1.0 / (2.0 * factor);
        for (int i = 1; i <= 2 * stateDim; i++) {
            wm.setEntry(i, w);
            wc.setEntry(i, w);
        }

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

    /** {@inheritDoc} */
    @Override
    protected double getMultiplicationFactor() {
        return factor;
    }

}
