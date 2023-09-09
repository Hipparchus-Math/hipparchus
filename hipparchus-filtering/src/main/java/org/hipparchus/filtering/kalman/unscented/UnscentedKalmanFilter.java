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
package org.hipparchus.filtering.kalman.unscented;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.filtering.kalman.KalmanFilter;
import org.hipparchus.filtering.kalman.Measurement;
import org.hipparchus.filtering.kalman.ProcessEstimate;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.MatrixDecomposer;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.util.UnscentedTransformProvider;

/**
 * Unscented Kalman filter for {@link UnscentedProcess unscented process}.
 * @param <T> the type of the measurements
 *
 * @see "Wan, E. A., & Van Der Merwe, R. (2000, October). The unscented Kalman filter for nonlinear estimation.
 *       In Proceedings of the IEEE 2000 Adaptive Systems for Signal Processing, Communications, and Control Symposium
 *       (Cat. No. 00EX373) (pp. 153-158)"
 * @since 2.2
 */
public class UnscentedKalmanFilter<T extends Measurement> implements KalmanFilter<T> {

    /** Process to be estimated. */
    private UnscentedProcess<T> process;

    /** Predicted state. */
    private ProcessEstimate predicted;

    /** Corrected state. */
    private ProcessEstimate corrected;

    /** Decompose to use for the correction phase. */
    private final MatrixDecomposer decomposer;

    /** Number of estimated parameters. */
    private final int n;

    /** Unscented transform provider. */
    private final UnscentedTransformProvider utProvider;

    /** Simple constructor.
     * @param decomposer decomposer to use for the correction phase
     * @param process unscented process to estimate
     * @param initialState initial state
     * @param utProvider unscented transform provider
     */
    public UnscentedKalmanFilter(final MatrixDecomposer decomposer,
                                 final UnscentedProcess<T> process,
                                 final ProcessEstimate initialState,
                                 final UnscentedTransformProvider utProvider) {
        this.decomposer = decomposer;
        this.process    = process;
        this.corrected  = initialState;
        this.n          = corrected.getState().getDimension();
        this.utProvider = utProvider;
        // Check state dimension
        if (n == 0) {
            // State dimension must be different from 0
            throw new MathIllegalArgumentException(LocalizedCoreFormats.ZERO_STATE_SIZE);
        }
    }

    /** {@inheritDoc} */
    @Override
    public ProcessEstimate estimationStep(final T measurement) throws MathRuntimeException {

        // Calculate sigma points
        final RealVector[] sigmaPoints = utProvider.unscentedTransform(corrected.getState(), corrected.getCovariance());

        // Perform the prediction and correction steps
        return predictionAndCorrectionSteps(measurement, sigmaPoints);

    }

    /** This method perform the prediction and correction steps of the Unscented Kalman Filter.
     * @param measurement single measurement to handle
     * @param sigmaPoints computed sigma points
     * @return estimated state after measurement has been considered
     * @throws MathRuntimeException if matrix cannot be decomposed
     */
    public ProcessEstimate predictionAndCorrectionSteps(final T measurement, final RealVector[] sigmaPoints) throws MathRuntimeException {

        // Prediction phase
        final UnscentedEvolution evolution = process.getEvolution(getCorrected().getTime(),
                                                                  sigmaPoints, measurement);

        predict(evolution.getCurrentTime(), evolution.getCurrentStates(),
                evolution.getProcessNoiseMatrix());

        // Calculate sigma points from predicted state
        final RealVector[] predictedSigmaPoints = utProvider.unscentedTransform(predicted.getState(),
                                                                                predicted.getCovariance());

        // Correction phase
        final RealVector[] predictedMeasurements = process.getPredictedMeasurements(predictedSigmaPoints, measurement);
        final RealVector   predictedMeasurement  = sum(predictedMeasurements, measurement.getValue().getDimension());
        final RealMatrix   r                     = computeInnovationCovarianceMatrix(predictedMeasurements, predictedMeasurement, measurement.getCovariance());
        final RealMatrix   crossCovarianceMatrix = computeCrossCovarianceMatrix(predictedSigmaPoints, predicted.getState(),
                                                                                predictedMeasurements, predictedMeasurement);
        final RealVector   innovation            = (r == null) ? null : process.getInnovation(measurement, predictedMeasurement, predicted.getState(), r);
        correct(measurement, r, crossCovarianceMatrix, innovation);
        return getCorrected();

    }

    /** Perform prediction step.
     * @param time process time
     * @param predictedStates predicted state vectors
     * @param noise process noise covariance matrix
     */
    private void predict(final double time, final RealVector[] predictedStates,  final RealMatrix noise) {

        // Computation of Eq. 17, weighted mean state
        final RealVector predictedState = sum(predictedStates, n);

        // Computation of Eq. 18, predicted covariance matrix
        final RealMatrix predictedCovariance = computeCovariance(predictedStates, predictedState).add(noise);

        predicted = new ProcessEstimate(time, predictedState, predictedCovariance);
        corrected = null;

    }

    /** Perform correction step.
     * @param measurement single measurement to handle
     * @param innovationCovarianceMatrix innovation covariance matrix
     * (may be null if measurement should be ignored)
     * @param crossCovarianceMatrix cross covariance matrix
     * @param innovation innovation
     * (may be null if measurement should be ignored)
     * @exception MathIllegalArgumentException if matrix cannot be decomposed
     */
    private void correct(final T measurement, final RealMatrix innovationCovarianceMatrix,
                           final RealMatrix crossCovarianceMatrix, final RealVector innovation)
        throws MathIllegalArgumentException {

        if (innovation == null) {
            // measurement should be ignored
            corrected = predicted;
            return;
        }

        // compute Kalman gain k
        // the following is equivalent to k = P_cross * (R_pred)^-1
        // we don't want to compute the inverse of a matrix,
        // we start by post-multiplying by R_pred and get
        // k.(R_pred) = P_cross
        // then we transpose, knowing that R_pred is a symmetric matrix
        // (R_pred).k^T = P_cross^T
        // then we can use linear system solving instead of matrix inversion
        final RealMatrix k = decomposer.
                             decompose(innovationCovarianceMatrix).
                             solve(crossCovarianceMatrix.transpose()).transpose();

        // correct state vector
        final RealVector correctedState = predicted.getState().add(k.operate(innovation));

        // correct covariance matrix
        final RealMatrix correctedCovariance = predicted.getCovariance().
                                               subtract(k.multiply(innovationCovarianceMatrix).multiplyTransposed(k));

        corrected = new ProcessEstimate(measurement.getTime(), correctedState, correctedCovariance,
                                        null, null, innovationCovarianceMatrix, k);

    }
    /** Get the predicted state.
     * @return predicted state
     */
    @Override
    public ProcessEstimate getPredicted() {
        return predicted;
    }

    /** Get the corrected state.
     * @return corrected state
     */
    @Override
    public ProcessEstimate getCorrected() {
        return corrected;
    }

    /** Get the unscented transform provider.
     * @return unscented transform provider
     */
    public UnscentedTransformProvider getUnscentedTransformProvider() {
        return utProvider;
    }

    /** Computes innovation covariance matrix.
     * @param predictedMeasurements predicted measurements (one per sigma point)
     * @param predictedMeasurement predicted measurements
     *        (may be null if measurement should be ignored)
     * @param r measurement covariance
     * @return innovation covariance matrix (null if predictedMeasurement is null)
     */
    private RealMatrix computeInnovationCovarianceMatrix(final RealVector[] predictedMeasurements,
                                                         final RealVector predictedMeasurement,
                                                         final RealMatrix r) {
        if (predictedMeasurement == null) {
            return null;
        }
        // Computation of the innovation covariance matrix
        final RealMatrix innovationCovarianceMatrix = computeCovariance(predictedMeasurements, predictedMeasurement);
        // Add the measurement covariance
        return innovationCovarianceMatrix.add(r);
    }

    /**
     * Computes cross covariance matrix.
     * @param predictedStates predicted states
     * @param predictedState predicted state
     * @param predictedMeasurements current measurements
     * @param predictedMeasurement predicted measurements
     * @return cross covariance matrix
     */
    private RealMatrix computeCrossCovarianceMatrix(final RealVector[] predictedStates, final RealVector predictedState,
                                                    final RealVector[] predictedMeasurements, final RealVector predictedMeasurement) {

        // Initialize the cross covariance matrix
        RealMatrix crossCovarianceMatrix = MatrixUtils.createRealMatrix(predictedState.getDimension(),
                                                                        predictedMeasurement.getDimension());

        // Covariance weights
        final RealVector wc = utProvider.getWc();

        // Compute the cross covariance matrix
        for (int i = 0; i <= 2 * n; i++) {
            final RealVector stateDiff = predictedStates[i].subtract(predictedState);
            final RealVector measDiff  = predictedMeasurements[i].subtract(predictedMeasurement);
            crossCovarianceMatrix = crossCovarianceMatrix.add(outer(stateDiff, measDiff).scalarMultiply(wc.getEntry(i)));
        }

        // Return the cross covariance
        return crossCovarianceMatrix;
    }

    /**
     * Computes a weighted mean parameter from a given samples.
     * <p>
     * This method can be used for computing both the mean state and the mean measurement.
     * <p>
     * It corresponds to the Equation 17 of "Wan, E. A., & Van Der Merwe, R.
     * The unscented Kalman filter for nonlinear estimation"
     * </p>
     * @param samples input samples
     * @param size size of the weighted mean parameter
     * @return weighted mean parameter
     */
    private RealVector sum(final RealVector[] samples, final int size) {

        // Initialize the weighted mean parameter
        RealVector mean = new ArrayRealVector(size);

        // Mean weights
        final RealVector wm = utProvider.getWm();

        // Compute weighted mean parameter
        for (int i = 0; i <= 2 * n; i++) {
            mean = mean.add(samples[i].mapMultiply(wm.getEntry(i)));
        }

        // Return the weighted mean value
        return mean;

    }

    /** Computes the covariance matrix.
     * <p>
     * This method can be used for computing both the predicted state
     * covariance matrix and the innovation covariance matrix.
     * <p>
     * It corresponds to the Equation 18 of "Wan, E. A., & Van Der Merwe, R.
     * The unscented Kalman filter for nonlinear estimation"
     * </p>
     * @param samples input samples
     * @param state weighted mean parameter
     * @return the covariance matrix
     */
    private RealMatrix computeCovariance(final RealVector[] samples,
                                         final RealVector state) {

        // Initialize the covariance matrix, by using the size of the weighted mean parameter
        final int dim = state.getDimension();
        RealMatrix covarianceMatrix = MatrixUtils.createRealMatrix(dim, dim);

        // Covariance weights
        final RealVector wc = utProvider.getWc();

        // Compute the covariance matrix
        for (int i = 0; i <= 2 * n; i++) {
            final RealVector diff = samples[i].subtract(state);
            covarianceMatrix = covarianceMatrix.add(outer(diff, diff).scalarMultiply(wc.getEntry(i)));
        }

        // Return the covariance
        return covarianceMatrix;

    }

    /** Conputes the outer product of two vectors.
     * @param a first vector
     * @param b second vector
     * @return the outer product of a and b
     */
    private RealMatrix outer(final RealVector a, final RealVector b) {

        // Initialize matrix
        final RealMatrix outMatrix = MatrixUtils.createRealMatrix(a.getDimension(), b.getDimension());

        // Fill matrix
        for (int row = 0; row < outMatrix.getRowDimension(); row++) {
            for (int col = 0; col < outMatrix.getColumnDimension(); col++) {
                outMatrix.setEntry(row, col, a.getEntry(row) * b.getEntry(col));
            }
        }

        // Return
        return outMatrix;

    }

}
