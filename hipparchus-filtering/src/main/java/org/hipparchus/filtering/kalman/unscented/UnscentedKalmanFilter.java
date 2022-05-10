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
import org.hipparchus.util.MerweUnscentedTransform;
import org.hipparchus.util.UnscentedTransformProvider;


/**
 * Unscented Kalman filter for {@link UscentedProcess unscented process}.
 * @param <T> the type of the measurements
 * @see "Woodburn J., and Coppola V., Analysis of Relative Merits of Unscented and Extended Kalman Filter in Orbit
 *       Determination, Reprinted from Astrodynamics 1999, Advances in the Astronautical Sciences, Vol. 171, 2019."
 */

public class UnscentedKalmanFilter<T extends Measurement> implements KalmanFilter<T> {

    /** Unscented process. */
    private UnscentedProcess<T> process;

    /** Predicted states. */
    private ProcessEstimate predicted;

    /** Corrected states. */
    private ProcessEstimate corrected;

    /** Decompose to use for the correction phase. */
    private final MatrixDecomposer decomposer;

    /** Number of estimated parameters. */
    private final int n;

    /** Unscend transform provider. */
    private final UnscentedTransformProvider utProvider;


    /**
     * Default constructor.
     * <p>
     * This constructor uses MerweUnscentedTransform as default transform for <code>utProvider</code>,
     * @param decomposer decomposer to use for the correction phase
     * @param process unscented process to estimate
     * @param initialState initial state
     */
    public UnscentedKalmanFilter(final MatrixDecomposer decomposer, final UnscentedProcess<T> process,
                                 final ProcessEstimate initialState) {
        this(decomposer, process, initialState, new MerweUnscentedTransform(initialState.getState().getDimension()));
    }

    /** Simple constructor.
     * @param decomposer decomposer to use for the correction phase
     * @param initialState initial state
     * @param process unscented process to estimate
     * @param utProvider unscented transform provider
     */

    public UnscentedKalmanFilter(final MatrixDecomposer decomposer, final UnscentedProcess<T> process,
                                 final  ProcessEstimate initialState,
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
        // Unscented transform
        final RealVector[] sigmaPoints = unscentedTransform();

        // Prediction phase
        final UnscentedEvolution evolution = process.getEvolution(getCorrected().getTime(),
                                                                  sigmaPoints,
                                                                  measurement);
        predict(evolution.getCurrentTime(), evolution.getCurrentStates(), evolution.getProcessNoiseMatrix());
        final RealVector predictedMeasurement = getMean(evolution.getCurrentMeasurements());
        final RealMatrix innovationCovarianceMatrix = getInnovationCovarianceMatrix(evolution.getCurrentMeasurements(), predictedMeasurement, measurement.getCovariance());
        final RealMatrix crossCovarianceMatrix = getCrossCovarianceMatrix(evolution.getCurrentStates(), predicted.getState(), evolution.getCurrentMeasurements(), predictedMeasurement);
        final RealVector innovation = (innovationCovarianceMatrix == null) ? null : process.getInnovation(measurement, predictedMeasurement, predicted.getState(), innovationCovarianceMatrix);
        // Correction phase
        correct(measurement, innovationCovarianceMatrix, crossCovarianceMatrix, innovation);
        return getCorrected();
    }

    /** Compute sigma points through unscented transform from previous state.
     * @return sigma points.
     */
    protected RealVector[] unscentedTransform() {

        return utProvider.unscentedTransform(corrected.getState(), corrected.getCovariance());

    }

    /** Perform prediction step.
     * @param time process time
     * @param predictedStates predicted states vector
     * @param noise process noise covariance matrix
     */
    protected void predict(final double time, final RealVector[] predictedStates,  final RealMatrix noise) {

        // Initialize predicted state and covariance
        RealVector predictedState      = new ArrayRealVector(n);
        RealMatrix predictedCovariance = MatrixUtils.createRealMatrix(n, n);

        // Computation of Eq. 23, weighted mean
        predictedState = getMean(predictedStates);

        // Computation of Eq. 24, covariance matrix
        predictedCovariance = getCovariance(predictedStates, predictedState);

        predicted = new ProcessEstimate(time, predictedState, predictedCovariance.add(noise));
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

    protected void correct(final T measurement, final RealMatrix innovationCovarianceMatrix,
                           final RealMatrix crossCovarianceMatrix, final RealVector innovation) throws MathIllegalArgumentException {

        if (innovation == null) {
            // measurement should be ignored
            corrected = predicted;
            return;
        }
        // compute Kalman gain k
        // the following is equivalent to k = P_cross * (R_pred)^-1 (Eq. 31)
        // we don't want to compute the inverse of a matrix,
        // we start by post-multiplying by R_pred and get
        // k.(R_pred) = P_cross
        // then we transpose, knowing that R_pred is a symmetric matrice
        // (R_pred).k^T = P_cross^T
        // then we can use linear system solving instead of matrix inversion
        final RealMatrix k = decomposer.decompose(innovationCovarianceMatrix).solve(crossCovarianceMatrix.transpose()).transpose();
        final RealVector correctedState = predicted.getState().add(k.operate(innovation));
        final RealMatrix correctedCovariance = predicted.getCovariance().subtract(k.multiply(innovationCovarianceMatrix).multiply(k.transpose()));

        corrected = new ProcessEstimate(measurement.getTime(), correctedState, correctedCovariance);

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

    /** Computes innovation covariance matrix. See Eq. 29.
     * @param currentMeasurements current measurements
     * @param predictedMeasurement predicted measurements
     * (may be null if measurement should be ignored)
     * @param r measurement covariance
     * @return innovation covariance matrix (null if predictedMeasurement is null)
     */
    protected RealMatrix getInnovationCovarianceMatrix(final RealVector[] currentMeasurements, final RealVector predictedMeasurement, final RealMatrix r) {
        if (predictedMeasurement == null) {
            return null;
        }
        // Measurement dimension
        final int measDim = predictedMeasurement.getDimension();

        RealMatrix innovationCovarianceMatrix = MatrixUtils.createRealMatrix(measDim, measDim);

        // Computation of Eq. 29
        innovationCovarianceMatrix = getCovariance(currentMeasurements, predictedMeasurement);

        return innovationCovarianceMatrix.add(r);
    }

    /**
     * Computes cross covariance matrix. See Eq. 30.
     * @param predictedStates predicted states
     * @param predictedState predicted state
     * @param predictedMeasurements current measurements
     * @param predictedMeasurement predicted measurement
     * @return cross covariance matrix
     */
    protected RealMatrix getCrossCovarianceMatrix(final RealVector[] predictedStates, final RealVector predictedState, final RealVector[] predictedMeasurements, final RealVector predictedMeasurement) {

        final int stateDim = predictedState.getDimension();
        final int measDim = predictedMeasurement.getDimension();
        final int dim = predictedStates.length;
        final RealVector wc = utProvider.getWc();
        // x*y^T matrix where x stands for the difference between predictedStates and the predictedState
        //                    y stands for the difference between predictedMeasurements and the predictedMeasurement
        final RealMatrix crossDiffSquared = MatrixUtils.createRealMatrix(stateDim, measDim);
        RealMatrix crossCovarianceMatrix = MatrixUtils.createRealMatrix(stateDim, measDim);

        for (int i = 0; i < dim; i++) {
            final RealVector stateDiff = predictedStates[i].subtract(predictedState);
            final RealVector measDiff = predictedMeasurements[i].subtract(predictedMeasurement);

            for (int c = 0; c < measDim; c++) {
                for (int l = 0; l < stateDim; l++) {
                    crossDiffSquared.setEntry(l, c, stateDiff.getEntry(l) * measDiff.getEntry(c));
                }
            }

            crossCovarianceMatrix = crossCovarianceMatrix.add(crossDiffSquared.scalarMultiply(wc.getEntry(i)));
        }
        return crossCovarianceMatrix;
    }

    /**
     * Computes weighted mean from samples. See Eq. 23 and 28.
     * @param samples
     * @return weighted mean
     */

    protected RealVector getMean(final RealVector[] samples) {
        
        final int dim = samples[0].getDimension();
        final int p = samples.length;
        final RealVector wm = utProvider.getWm();
        
        RealVector mean = new ArrayRealVector(dim);
        for (int i = 0; i < p; i++) {
            mean = mean.add(samples[i].mapMultiply(wm.getEntry(i)));
        }

        return mean;
    }
    
    /** Computes covariance from state and samples. See Eq. 24 and 29.
     * @param samples samples
     * @param state state
     * @return covariance matrix
     */
    protected RealMatrix getCovariance(final RealVector[] samples, final RealVector state) {
        // dim can be either state dimension or measurement dimension. It depends on the covariance one wants to compute.
        final int dim = state.getDimension();
        final int p = samples.length;
        final RealVector wc = utProvider.getWc();
        // y*y^T matrix where y stands for the difference between samples and state
        final RealMatrix diffSquared = MatrixUtils.createRealMatrix(dim, dim);

        RealMatrix covarianceMatrix = MatrixUtils.createRealMatrix(dim, dim);

        // Computation of Eq. 18
        for (int i = 0; i < p; i++) {
            final RealVector diff = samples[i].subtract(state);
            for (int c = 0; c < dim; c++) {
                for (int l = 0; l < dim; l++) {
                    diffSquared.setEntry(l, c, diff.getEntry(l) * diff.getEntry(c));
                }
            }
            covarianceMatrix = covarianceMatrix.add(diffSquared.scalarMultiply(wc.getEntry(i)));
        }
        return covarianceMatrix;
    }


}
