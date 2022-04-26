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
import org.hipparchus.linear.CholeskyDecomposition;
import org.hipparchus.linear.MatrixDecomposer;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;


/**
 * Unscented Kalman filter for {@link UscentedProcess unscented process}.
 * @param <T> the type of the measurements
 * @see "Woodburn J., and Coppola V., Analysis of Relative Merits of Unscented and Extended Kalman Filter in Orbit
 *       Determination, Reprinted from Astrodynamics 1999, Advances in the Astronautical Sciences, Vol. 171, 2019."
 */

public class UnscentedKalmanFilter<T extends Measurement> implements KalmanFilter<T> {

    /** Default value for alpha (1.0E-4), see reference. */
    public static final double DEFAULT_ALPHA = 0.0001;

    /** Default value for beta (2.), see reference. */
    public static final double DEFAULT_BETA = 2;

    /** Default value for kappa, (0.), see reference. */
    public static final double DEFAULT_KAPPA = 0;

    /** Unscented process. */
    private UnscentedProcess<T> process;

    /** Predicted states. */
    private ProcessEstimate predicted;

    /** Corrected states. */
    private ProcessEstimate corrected;

    /** Decompose to use for the correction phase. */
    private final MatrixDecomposer decomposer;

    /** Factor applied during Unscented transform for covariance matrix. */
    private final double factor;

    /** Number of estimated parameters. */
    private final int n;

    /** Weights for covariance matrix. */
    private final RealVector wc;

    /** Weights for estimation. */
    private final RealVector wm;


    /**
     * Default constructor.
     * <p>
     * This constructor uses default values for <code>alpha</code>,
     * <code>beta</code>, and <code>kappa</code>.
     * Default alpha is equal to {@link #DEFAULT_ALPHA}, default beta is equal to {@link #DEFAULT_BETA} and default kappa is equal to {@link #DEFAULT_KAPPA}.
     * </p>
     * @param decomposer decomposer to use for the correction phase
     * @param process unscented process to estimate
     * @param initialState initial state
     */
    public UnscentedKalmanFilter(final MatrixDecomposer decomposer, final UnscentedProcess<T> process,
                                 final ProcessEstimate initialState) {
        this(decomposer, process, initialState, DEFAULT_ALPHA, DEFAULT_BETA, DEFAULT_KAPPA);
    }

    /** Simple constructor.
     * @param decomposer decomposer to use for the correction phase
     * @param initialState initial state
     * @param process unscented process to estimate
     * @param alpha scaling control
     * @param beta free parameter
     * @param kappa free parameter
     */

    public UnscentedKalmanFilter(final MatrixDecomposer decomposer, final UnscentedProcess<T> process,
                                 final  ProcessEstimate initialState,
                                 final double alpha, final double beta, final double kappa) {

        this.decomposer = decomposer;
        this.process    = process;
        this.corrected  = initialState;
        this.n          = corrected.getState().getDimension();

        // Check state dimension
        if (n == 0) {
            // State dimension must be different from 0
            throw new MathIllegalArgumentException(LocalizedCoreFormats.ZERO_STATE_SIZE);
        }

        // lambda = alpha² + (n + kappa) - n (see Eq. 10)
        final double lambda = alpha * alpha * (n + kappa) - n;

        // Initialize multiplication factor for covariance matrix
        this.factor = n + lambda;

        // Initialize vectors containing unscented kalman filter weights
        wm = new ArrayRealVector(2 * n + 1);
        wc = new ArrayRealVector(2 * n + 1);

        // Computation of unscented kalman filter weights (See Eq. 15, Eq. 16 and Eq. 17)
        wm.setEntry(0, lambda / (n + lambda));
        wc.setEntry(0, lambda / (n + lambda) + 1 - alpha * alpha + beta);
        for (int i = 1; i <= 2 * n; i++) {
            wm.setEntry(i, 1 / (2 * (n + lambda)));
            wc.setEntry(i, 1 / (2 * (n + lambda)));
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

        // Correction phase
        correct(measurement, evolution.getCurrentMeasurements(), evolution.getCurrentStates(), measurement.getCovariance());
        return getCorrected();
    }

    /** Compute sigma points through unscented transform (cf. Eq.12 and Eq.13) from previous state.
     * @return sigma points.
     */
    protected RealVector[] unscentedTransform() {

        // Initialize array containing sigma points
        final RealVector[] sigmaPoints = new ArrayRealVector[(2 * n) + 1];
        sigmaPoints[0] = corrected.getState();
        final RealMatrix temp = corrected.getCovariance().scalarMultiply(factor);

        // Compute lower triangular matrix of Cholesky decomposition
        final CholeskyDecomposition chdecomposition = new CholeskyDecomposition(temp, 10e-14, 10e-14);
        final RealMatrix L = chdecomposition.getL();

        // Compute sigma points
        for (int i = 1; i <= n; i++) {

            // Computation of Eq. 12
            sigmaPoints[i] = sigmaPoints[0].add(L.getColumnVector(i - 1));

            // Computation of Eq. 13
            sigmaPoints[i + n] = sigmaPoints[0].subtract(L.getColumnVector(i - 1));

        }

        // Return sigma points
        return sigmaPoints;
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

        // x*x^T matrix where x stands for the difference between predicted states and the predicted state
        final RealMatrix statesDiffSquared = MatrixUtils.createRealMatrix(n, n);

        // Computation of Eq. 23, weighted mean
        for (int i = 0; i <= 2 * n; i++) {
            predictedState = predictedState.add(predictedStates[i].mapMultiply(wm.getEntry(i)));
        }

        // Computation of Eq. 24
        for (int i = 0; i <= 2 * n; i++) {
            final RealVector stateDiff = predictedStates[i].subtract(predictedState);
            for (int c = 0; c < n; c++) {
                for (int l = 0; l < n; l++) {
                    statesDiffSquared.setEntry(l, c, stateDiff.getEntry(l) * stateDiff.getEntry(c));
                }
            }

            predictedCovariance = predictedCovariance.add(statesDiffSquared.scalarMultiply(wc.getEntry(i)));
        }

        predicted = new ProcessEstimate(time, predictedState, predictedCovariance.add(noise));
        corrected = null;
    }

    /** Perform correction step.
     * @param measurement single measurement to handle
     * @param currentMeasurements current measurements
     * (may be null if measurement should be ignored)
     * @param predictedStates predicted states
     * @param r innovation covariance matrix
     * (may be null if measurement should be ignored)
     * @exception MathIllegalArgumentException if matrix cannot be decomposed
     */

    protected void correct(final T measurement, final RealVector[] currentMeasurements,
                           final RealVector[] predictedStates, final RealMatrix r) throws MathIllegalArgumentException {

        if (currentMeasurements == null) {
            // measurement should be ignored
            corrected = predicted;
            return;
        }

        // Initialize the predicted measurement
        final int measDim = currentMeasurements[0].getDimension();
        RealVector predictedMeasurement = new ArrayRealVector(measDim);

        // Computation of Eq. 28, weighted mean
        for (int i = 0; i <= 2 * n; i++) {
            predictedMeasurement = predictedMeasurement.add(currentMeasurements[i].mapMultiply(wm.getEntry(i)));
        }

        // y*y^T matrix where y stands for the difference between current measurements and the predicted measurement
        final RealMatrix measurementDiffSquared = MatrixUtils.createRealMatrix(measDim, measDim);
        // x*y^T matrix where x stands for the difference between predicted states and the predicted state
        //                    y is defined above
        final RealMatrix crossDiffSquared = MatrixUtils.createRealMatrix(n, measDim);
        RealMatrix innovationCovarianceMatrix = MatrixUtils.createRealMatrix(measDim, measDim);
        RealMatrix crossCovarianceMatrix = MatrixUtils.createRealMatrix(n, measDim);



        // Computation of Eq. 29 and Eq. 30
        final RealVector innovation = measurement.getValue().subtract(predictedMeasurement);

        for (int i = 0; i <= 2 * n; i++) {
            final RealVector measurementDiff = currentMeasurements[i].subtract(predictedMeasurement);
            final RealVector stateDiff = predictedStates[i].subtract(predicted.getState());

            for (int c = 0; c < measDim; c++) {
                for (int l = 0; l < measDim; l++) {
                    measurementDiffSquared.setEntry(l, c, measurementDiff.getEntry(l) * measurementDiff.getEntry(c));
                }
                for (int l = 0; l < n; l++) {
                    crossDiffSquared.setEntry(l, c, stateDiff.getEntry(l) * measurementDiff.getEntry(c));
                }
            }

            innovationCovarianceMatrix = innovationCovarianceMatrix.add(measurementDiffSquared.scalarMultiply(wc.getEntry(i)));
            crossCovarianceMatrix = crossCovarianceMatrix.add(crossDiffSquared.scalarMultiply(wc.getEntry(i)));
        }
        // compute Kalman gain k
        // the following is equivalent to k = P_cross * (R_pred)^-1 (Eq. 31)
        // we don't want to compute the inverse of a matrix,
        // we start by post-multiplying by R_pred and get
        // k.(R_pred) = P_cross
        // then we transpose, knowing that R_pred is a symmetric matrice
        // (R_pred).k^T = P_cross^T
        // then we can use linear system solving instead of matrix inversion
        final RealMatrix innovationCovarianceMatrixPredicted = innovationCovarianceMatrix.add(r);
        final RealMatrix k = decomposer.decompose(innovationCovarianceMatrixPredicted).solve(crossCovarianceMatrix.transpose()).transpose();
        final RealVector correctedState = predicted.getState().add(k.operate(innovation));
        final RealMatrix correctedCovariance = predicted.getCovariance().subtract(k.multiply(innovationCovarianceMatrixPredicted).multiply(k.transpose()));

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

}
