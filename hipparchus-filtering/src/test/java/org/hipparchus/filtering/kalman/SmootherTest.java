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

import org.hipparchus.filtering.kalman.extended.ExtendedKalmanFilter;
import org.hipparchus.filtering.kalman.extended.NonLinearEvolution;
import org.hipparchus.filtering.kalman.extended.NonLinearProcess;
import org.hipparchus.filtering.kalman.linear.LinearEvolution;
import org.hipparchus.filtering.kalman.linear.LinearKalmanFilter;
import org.hipparchus.filtering.kalman.linear.LinearProcess;
import org.hipparchus.filtering.kalman.unscented.UnscentedEvolution;
import org.hipparchus.filtering.kalman.unscented.UnscentedKalmanFilter;
import org.hipparchus.filtering.kalman.unscented.UnscentedProcess;
import org.hipparchus.linear.*;
import org.hipparchus.util.MerweUnscentedTransform;
import org.hipparchus.util.UnscentedTransformProvider;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class SmootherTest {

    private static final double PROCESS_NOISE = 0.1;
    private static List<Reference> referenceData;
    private static List<SimpleMeasurement> measurements;
    private static double initialTime;
    private static ProcessEstimate initialState;
    private static MatrixDecomposer decomposer;

    @BeforeAll
    static void beforeAll() {

        // Load reference data
        referenceData = Reference.loadReferenceData(2, 1, "cv-smoother.txt");

        // Measurements (skip first one corresponding to smoothed initial state)
        final RealMatrix measurementNoise = MatrixUtils.createRealMatrix(new double[][]{{1e-3}});
        measurements = referenceData.stream()
                .skip(1)
                .map(r -> new SimpleMeasurement(r.getTime(), r.getZ(), measurementNoise))
                .collect(Collectors.toList());

        // Initial state
        initialTime = 0.0;
        final RealVector initialMean = MatrixUtils.createRealVector(new double[]{0.0, -0.5});
        final RealMatrix initialCovariance = MatrixUtils.createRealMatrix(new double[][]{{0.01, 0.0}, {0.0, 0.25}});
        initialState = new ProcessEstimate(initialTime, initialMean, initialCovariance);

        // Matrix decomposer for filter and smoother
        decomposer = new CholeskyDecomposer(1e-15, 1e-15);
    }

    // Test linear Kalman filter/smoother
    // Simple linear 1D constant velocity
    @Test
    void testKalmanSmoother() {

        // Initialise filter process
        final LinearProcess<SimpleMeasurement> process = new LinearConstantVelocity<>(initialTime, PROCESS_NOISE);

        // Smoother wrapping the filter
        final KalmanFilterSmoother<SimpleMeasurement> smoother =
                new KalmanFilterSmoother<>(new LinearKalmanFilter<>(decomposer, process, initialState), decomposer);

        // Process measurements with smoother (forwards pass)
        measurements.forEach(smoother::estimationStep);

        // Smooth backwards
        List<ProcessEstimate> smoothedStates = smoother.backwardsSmooth();

        // Check against reference
        Assertions.assertEquals(referenceData.size(), smoothedStates.size());
        for (int i = 0; i < referenceData.size(); ++i) {
            referenceData.get(i).checkState(smoothedStates.get(i).getState(), 1e-14);
            referenceData.get(i).checkCovariance(smoothedStates.get(i).getCovariance(), 1e-14);
        }
    }

    // Test extended Kalman filter/smoother
    // Simple linear 1D constant velocity
    @Test
    void testExtendedSmoother() {

        // Initialise filter
        final NonLinearProcess<SimpleMeasurement> process = new ExtendedConstantVelocity<>(initialTime, PROCESS_NOISE);

        // Smoother wrapping the filter
        final KalmanFilterSmoother<SimpleMeasurement> smoother =
                new KalmanFilterSmoother<>(new ExtendedKalmanFilter<>(decomposer, process, initialState), decomposer);

        // Process measurements with smoother (forwards pass)
        measurements.forEach(smoother::estimationStep);

        // Smooth backwards
        List<ProcessEstimate> smoothedStates = smoother.backwardsSmooth();

        // Check against reference
        Assertions.assertEquals(referenceData.size(), smoothedStates.size());
        for (int i = 0; i < referenceData.size(); ++i) {
            referenceData.get(i).checkState(smoothedStates.get(i).getState(), 1e-14);
            referenceData.get(i).checkCovariance(smoothedStates.get(i).getCovariance(), 1e-14);
        }
    }


    // Test unscented Kalman filter/smoother
    // Simple linear 1D constant velocity
    @Test
    void testUnscentedSmoother() {

        // Initialise filter
        final UnscentedProcess<SimpleMeasurement> process = new UnscentedConstantVelocity<>(initialTime, PROCESS_NOISE);

        // Smoother wrapping the filter
        final double alpha = 1.0;
        final double beta = 2.0;
        final double kappa = 0.0;
        UnscentedTransformProvider unscentedTransformProvider =
                new MerweUnscentedTransform(initialState.getState().getDimension(), alpha, beta, kappa);
        final KalmanFilterSmoother<SimpleMeasurement> smoother = new KalmanFilterSmoother<>(
                new UnscentedKalmanFilter<>(decomposer, process, initialState, unscentedTransformProvider), decomposer);

        // Process measurements with smoother (forwards pass)
        measurements.forEach(smoother::estimationStep);

        // Smooth backwards
        List<ProcessEstimate> smoothedStates = smoother.backwardsSmooth();

        // Check against reference
        Assertions.assertEquals(referenceData.size(), smoothedStates.size());
        for (int i = 0; i < referenceData.size(); ++i) {
            referenceData.get(i).checkState(smoothedStates.get(i).getState(), 1e-14);
            referenceData.get(i).checkCovariance(smoothedStates.get(i).getCovariance(), 1e-14);
        }
    }


    // Process model for the linear Kalman filter
    private static class LinearConstantVelocity<T extends Measurement> implements LinearProcess<T> {

        private double currentTime;
        private final double processNoiseScale;

        public LinearConstantVelocity(final double initialTime, final double processNoiseScale) {
            currentTime = initialTime;
            this.processNoiseScale = processNoiseScale;
        }

        @Override
        public LinearEvolution getEvolution(final T measurement) {
            final double dt = measurement.getTime() - currentTime;
            currentTime = measurement.getTime();

            return new LinearEvolution(
                    getStateTransitionMatrix(dt),
                    MatrixUtils.createRealMatrix(2, 2),
                    MatrixUtils.createRealVector(2),
                    getProcessNoise(dt, processNoiseScale),
                    getMeasurementJacobian()
            );
        }
    }

    // Process model for the extended Kalman filter
    private static class ExtendedConstantVelocity<T extends Measurement> implements NonLinearProcess<T> {

        private double currentTime;
        private final double processNoiseScale;

        public ExtendedConstantVelocity(final double initialTime, final double processNoiseScale) {
            currentTime = initialTime;
            this.processNoiseScale = processNoiseScale;
        }

        @Override
        public NonLinearEvolution getEvolution(final double previousTime,
                                               final RealVector previousState,
                                               final T measurement) {
            final double dt = measurement.getTime() - currentTime;
            currentTime = measurement.getTime();

            final RealMatrix stateTransitionMatrix = getStateTransitionMatrix(dt);
            return new NonLinearEvolution(
                    currentTime,
                    stateTransitionMatrix.operate(previousState),
                    stateTransitionMatrix,
                    getProcessNoise(dt, processNoiseScale),
                    getMeasurementJacobian()
            );
        }

        @Override
        public RealVector getInnovation(final T measurement,
                                        final NonLinearEvolution evolution,
                                        final RealMatrix innovationCovarianceMatrix) {
            return measurement.getValue().subtract(evolution.getMeasurementJacobian().operate(evolution.getCurrentState()));
        }
    }

    // Process model for the unscented Kalman filter
    private static class UnscentedConstantVelocity<T extends Measurement> implements UnscentedProcess<T> {

        private double currentTime;
        private final double processNoiseScale;

        public UnscentedConstantVelocity(final double initialTime, final double processNoiseScale) {
            currentTime = initialTime;
            this.processNoiseScale = processNoiseScale;
        }

        @Override
        public UnscentedEvolution getEvolution(final double previousTime,
                                               final RealVector[] sigmaPoints,
                                               final T measurement) {
            final double dt = measurement.getTime() - currentTime;
            currentTime = measurement.getTime();

            final RealMatrix stateTransitionMatrix = getStateTransitionMatrix(dt);

            final RealVector[] predictedSigmaPoints = Arrays.stream(sigmaPoints)
                    .map(stateTransitionMatrix::operate)
                    .toArray(RealVector[]::new);

            return new UnscentedEvolution(currentTime, predictedSigmaPoints);
        }

        @Override
        public RealMatrix getProcessNoiseMatrix(double previousTime, RealVector predictedState, T measurement) {
            final double dt = measurement.getTime() - previousTime;
            return getProcessNoise(dt, processNoiseScale);
        }

        @Override
        public RealVector[] getPredictedMeasurements(final RealVector[] predictedSigmaPoints,
                                                     final T measurement) {

            final RealMatrix measurementJacobian = getMeasurementJacobian();
            return Arrays.stream(predictedSigmaPoints)
                    .map(measurementJacobian::operate)
                    .toArray(RealVector[]::new);
        }

        @Override
        public RealVector getInnovation(final T measurement,
                                        final RealVector predictedMeasurement,
                                        final RealVector predictedState,
                                        final RealMatrix innovationCovarianceMatrix) {
            return measurement.getValue().subtract(predictedMeasurement);
        }
    }


    // Common process model quantities
    private static RealMatrix getStateTransitionMatrix(final double dt) {
        return MatrixUtils.createRealMatrix(new double[][]{{1.0, dt}, {0.0, 1.0}});
    }

    private static RealMatrix getProcessNoise(final double dt, final double processNoiseScale) {
        final double dt2 = dt * dt;
        final double dt3 = dt2 * dt;

        return MatrixUtils.createRealMatrix(new double[][]{{dt3 / 3.0, dt2 / 2.0}, {dt2 / 2.0, dt}})
                .scalarMultiply(processNoiseScale);
    }

    private static RealMatrix getMeasurementJacobian() {
        return MatrixUtils.createRealMatrix(new double[][]{{1.0, 0.0}});
    }
}
