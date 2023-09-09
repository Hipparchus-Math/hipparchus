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

/** Unscented transform as defined by Julier and Uhlmann.
 * <p>
 * The unscented transform uses three parameters: alpha, beta and kappa.
 * Alpha determines the spread of the sigma points around the process state,
 * kappa is a secondary scaling parameter, and beta is used to incorporate
 * prior knowledge of the distribution of the process state.
 * <p>
 * The Julier transform is a particular case of {@link MerweUnscentedTransform} with alpha = 1 and beta = 0.
 * </p>
 * @see "S. J. Julier and J. K. Uhlmann. A New Extension of the Kalman Filter to Nonlinear Systems.
 *       Proc. SPIE 3068, Signal Processing, Sensor Fusion, and Target Recognition VI, 182 (July 28, 1997)"
 * @since 2.2
 */
public class JulierUnscentedTransform extends AbstractUnscentedTransform {

    /** Default value for kappa, (0.0, see reference). */
    public static final double DEFAULT_KAPPA = 0;

    /** Weights for covariance matrix. */
    private final RealVector wc;

    /** Weights for mean state. */
    private final RealVector wm;

    /** Factor applied to the covariance matrix during the unscented transform (lambda + process state size). */
    private final double factor;

    /**
     * Default constructor.
     * <p>
     * This constructor uses default value for kappa.
     * </p>
     * @param stateDim the dimension of the state
     * @see #DEFAULT_KAPPA
     * @see #JulierUnscentedTransform(int, double)
     */
    public JulierUnscentedTransform(final int stateDim) {
        this(stateDim, DEFAULT_KAPPA);
    }

    /**
     * Simple constructor.
     * @param stateDim the dimension of the state
     * @param kappa fscaling factor
     */
    public JulierUnscentedTransform(final int stateDim, final double kappa) {

        // Call super constructor
        super(stateDim);

        // Initialize multiplication factor for covariance matrix
        this.factor = stateDim + kappa;

        // Initialize vectors weights
        wm = new ArrayRealVector(2 * stateDim + 1);

        // Computation of unscented kalman filter weights (See Eq. 12)
        wm.setEntry(0, kappa / factor);
        for (int i = 1; i <= 2 * stateDim; i++) {
            wm.setEntry(i, 1.0 / (2.0 * factor));
        }

        // For the Julier unscented transform, there is no difference between covariance and state weights
        wc = wm;

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
