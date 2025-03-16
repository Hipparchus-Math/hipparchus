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
package org.hipparchus.filtering.kalman;

import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.filtering.LocalizedFilterFormats;
import org.hipparchus.linear.MatrixDecomposer;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

/**
 * Kalman smoother for linear, extended or unscented filters.
 * <p>
 * This implementation is attached to a filter using the observer mechanism.  Once all measurements have been
 * processed by the filter, the smoothing method can be called.
 * </p>
 * <p>
 * For example
 * </p>
 * <pre>{@code
 *     // Kalman filter
 *     final KalmanFilter<SimpleMeasurement> filter = new LinearKalmanFilter<>(decomposer, process, initialState);
 *
 *     // Smoother observer
 *     final KalmanSmoother smoother = new KalmanSmoother(decomposer);
 *     filter.setObserver(smoother);
 *
 *     // Process measurements with filter (forwards pass)
 *     measurements.forEach(filter::estimationStep);
 *
 *     // Smooth backwards
 *     List<ProcessEstimate> smoothedStates = smoother.backwardsSmooth();
 * }</pre>
 *
 * @see "Särkkä, S. Bayesian Filtering and Smoothing. Cambridge 2013"
 */
public class KalmanSmoother implements KalmanObserver {

    /** Decomposer to use for gain calculation. */
    private final MatrixDecomposer decomposer;

    /** Storage for smoother gain matrices. */
    private final List<SmootherData> smootherData;

    /** Simple constructor.
     * @param decomposer decomposer to use for the smoother gain calculations
     */
    public KalmanSmoother(final MatrixDecomposer decomposer) {
        this.decomposer = decomposer;
        this.smootherData = new ArrayList<>();
    }

    @Override
    public void init(KalmanEstimate estimate) {
        // Add initial state to smoother data
        smootherData.add(new SmootherData(
                estimate.getCorrected().getTime(),
                null,
                null,
                estimate.getCorrected().getState(),
                estimate.getCorrected().getCovariance(),
                null
        ));

    }

    @Override
    public void updatePerformed(KalmanEstimate estimate) {
        // Smoother gain
        // We want G = D * P^(-1)
        // Calculate with G = (P^(-1) * D^T)^T
        final RealMatrix smootherGain = decomposer
                .decompose(estimate.getPredicted().getCovariance())
                .solve(estimate.getStateCrossCovariance().transpose())
                .transpose();
        smootherData.add(new SmootherData(
                estimate.getCorrected().getTime(),
                estimate.getPredicted().getState(),
                estimate.getPredicted().getCovariance(),
                estimate.getCorrected().getState(),
                estimate.getCorrected().getCovariance(),
                smootherGain
        ));
    }

    /** Backwards smooth.
     * This is a backward pass over the filtered data, recursively calculating smoothed states, using the
     * Rauch-Tung-Striebel (RTS) formulation.
     * Note that the list result is a `LinkedList`, not an `ArrayList`.
     * @return list of smoothed states
     */
    public List<ProcessEstimate> backwardsSmooth() {
        // Check for at least one measurement
        if (smootherData.size() < 2) {
            throw new MathIllegalStateException(LocalizedFilterFormats.PROCESS_AT_LEAST_ONE_MEASUREMENT);
        }

        // Initialise output
        final LinkedList<ProcessEstimate> smootherResults = new LinkedList<>();

        // Last smoothed state is the same as the filtered state
        final SmootherData lastUpdate = smootherData.get(smootherData.size() - 1);
        ProcessEstimate smoothedState = new ProcessEstimate(lastUpdate.getTime(),
                lastUpdate.getCorrectedState(), lastUpdate.getCorrectedCovariance());
        smootherResults.addFirst(smoothedState);

        // Backwards recursion on the smoothed state
        for (int i = smootherData.size() - 2; i >= 0; --i) {

            // These are from equation 8.6 in Sarkka, "Bayesian Filtering and Smoothing", Cambridge, 2013.
            final RealMatrix smootherGain = smootherData.get(i + 1).getSmootherGain();

            final RealVector smoothedMean = smootherData.get(i).getCorrectedState()
                    .add(smootherGain.operate(smoothedState.getState()
                            .subtract(smootherData.get(i + 1).getPredictedState())));

            final RealMatrix smoothedCovariance = smootherData.get(i).getCorrectedCovariance()
                    .add(smootherGain.multiply(smoothedState.getCovariance()
                                    .subtract(smootherData.get(i + 1).getPredictedCovariance()))
                            .multiplyTransposed(smootherGain));

            // Populate smoothed state
            smoothedState = new ProcessEstimate(smootherData.get(i).getTime(), smoothedMean, smoothedCovariance);
            smootherResults.addFirst(smoothedState);
        }

        return smootherResults;
    }

    /** Container for smoother data. */
    private static class SmootherData {
        /** Process time (typically the time or index of a measurement). */
        private final double time;

        /** Predicted state vector. */
        private final RealVector predictedState;

        /** Predicted covariance. */
        private final RealMatrix predictedCovariance;

        /** Corrected state vector. */
        private final RealVector correctedState;

        /** Corrected covariance. */
        private final RealMatrix correctedCovariance;

        /** Smoother gain. */
        private final RealMatrix smootherGain;

        SmootherData(final double time,
                     final RealVector predictedState,
                     final RealMatrix predictedCovariance,
                     final RealVector correctedState,
                     final RealMatrix correctedCovariance,
                     final RealMatrix smootherGain) {
            this.time = time;
            this.predictedState = predictedState;
            this.predictedCovariance = predictedCovariance;
            this.correctedState = correctedState;
            this.correctedCovariance = correctedCovariance;
            this.smootherGain = smootherGain;
        }

        /** Get the process time.
         * @return process time (typically the time or index of a measurement)
         */
        public double getTime() {
            return time;
        }

        /**
         * Get predicted state
         * @return predicted state
         */
        public RealVector getPredictedState() {
            return predictedState;
        }

        /**
         * Get predicted covariance
         * @return predicted covariance
         */
        public RealMatrix getPredictedCovariance() {
            return predictedCovariance;
        }

        /**
         * Get corrected state
         * @return corrected state
         */
        public RealVector getCorrectedState() {
            return correctedState;
        }

        /**
         * Get corrected covariance
         * @return corrected covariance
         */
        public RealMatrix getCorrectedCovariance() {
            return correctedCovariance;
        }

        /**
         * Get smoother gain (for previous time-step)
         * @return smoother gain
         */
        public RealMatrix getSmootherGain() {
            return smootherGain;
        }
    }

}
