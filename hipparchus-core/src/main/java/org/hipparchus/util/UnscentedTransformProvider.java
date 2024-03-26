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
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;

/**
 * Provider for unscented transform.
 * @since 2.2
 */
public interface UnscentedTransformProvider {

    /**
     * Perform the unscented transform from a state and its covariance.
     * @param state process state
     * @param covariance covariance associated with the process state
     * @return an array containing the sigma points of the unscented transform
     */
    RealVector[] unscentedTransform(RealVector state, RealMatrix covariance);

    /**
     * Computes a weighted mean state from a given set of sigma points.
     * <p>
     * This method can be used for computing both the mean state and the mean measurement
     * in an Unscented Kalman filter.
     * </p>
     * <p>
     * It corresponds to Equation 17 of "Wan, E. A., &amp; Van Der Merwe, R. The unscented Kalman filter for nonlinear estimation"
     * </p>
     * @param sigmaPoints input samples
     * @return weighted mean state
     */
    default RealVector getUnscentedMeanState(RealVector[] sigmaPoints) {

        // Sigma point dimension
        final int sigmaPointDimension = sigmaPoints[0].getDimension();

        // Compute weighted mean
        // ---------------------

        RealVector weightedMean = new ArrayRealVector(sigmaPointDimension);

        // Compute the weight coefficients wm
        final RealVector wm = getWm();

        // Weight each sigma point and sum them
        for (int i = 0; i < sigmaPoints.length; i++) {
            weightedMean = weightedMean.add(sigmaPoints[i].mapMultiply(wm.getEntry(i)));
        }

        return weightedMean;
    }

    /** Computes the unscented covariance matrix from a weighted mean state and a set of sigma points.
     * <p>
     * This method can be used for computing both the predicted state
     * covariance matrix and the innovation covariance matrix in an Unscented Kalman filter.
     * </p>
     * <p>
     * It corresponds to Equation 18 of "Wan, E. A., &amp; Van Der Merwe, R. The unscented Kalman filter for nonlinear estimation"
     * </p>
     * @param sigmaPoints input sigma points
     * @param meanState weighted mean state
     * @return the unscented covariance matrix
     */
    default RealMatrix getUnscentedCovariance(RealVector[] sigmaPoints, RealVector meanState) {

        // State dimension
        final int stateDimension = meanState.getDimension();

        // Compute covariance matrix
        // -------------------------

        RealMatrix covarianceMatrix = MatrixUtils.createRealMatrix(stateDimension, stateDimension);

        // Compute the weight coefficients wc
        final RealVector wc = getWc();

        // Reconstruct the covariance
        for (int i = 0; i < sigmaPoints.length; i++) {
            final RealMatrix diff = MatrixUtils.createColumnRealMatrix(sigmaPoints[i].subtract(meanState).toArray());
            covarianceMatrix = covarianceMatrix.add(diff.multiplyTransposed(diff).scalarMultiply(wc.getEntry(i)));
        }

        return covarianceMatrix;
    }

    /**
     * Perform the inverse unscented transform from an array of sigma points.
     * @param sigmaPoints array containing the sigma points of the unscented transform
     * @return mean state and associated covariance
     */
    default Pair<RealVector, RealMatrix> inverseUnscentedTransform(RealVector[] sigmaPoints) {

        // Mean state
        final RealVector meanState = getUnscentedMeanState(sigmaPoints);

        // Return state and covariance
        return new Pair<>(meanState, getUnscentedCovariance(sigmaPoints, meanState));
    }

    /**
     * Get the covariance weights.
     * @return the covariance weights
     */
    RealVector getWc();

    /**
     * Get the mean weights.
     * @return the mean weights
     */
    RealVector getWm();

}
