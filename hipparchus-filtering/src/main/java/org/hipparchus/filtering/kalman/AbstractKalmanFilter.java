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

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.linear.MatrixDecomposer;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;

/**
 * Shared parts between linear and non-linear Kalman filters.
 * @param <T> the type of the measurements
 * @since 1.3
 */
public abstract class AbstractKalmanFilter<T extends Measurement> implements KalmanFilter<T> {

    /** Decomposer decomposer to use for the correction phase. */
    private final MatrixDecomposer decomposer;

    /** Predicted state. */
    private ProcessEstimate predicted;

    /** Corrected state. */
    private ProcessEstimate corrected;

    /** Simple constructor.
     * @param decomposer decomposer to use for the correction phase
     * @param initialState initial state
     */
    protected AbstractKalmanFilter(final MatrixDecomposer decomposer, final ProcessEstimate initialState) {
        this.decomposer = decomposer;
        this.corrected  = initialState;
    }

    /** Perform prediction step.
     * @param time process time
     * @param predictedState predicted state vector
     * @param stm state transition matrix
     * @param noise process noise covariance matrix
     */
    protected void predict(final double time, final RealVector predictedState, final RealMatrix stm, final RealMatrix noise) {
        final RealMatrix predictedCovariance =
                        stm.multiply(corrected.getCovariance().multiplyTransposed(stm)).add(noise);
        predicted = new ProcessEstimate(time, predictedState, predictedCovariance);
        corrected = null;
    }

    /** Compute innovation covariance matrix.
     * @param r measurement covariance
     * @param h Jacobian of the measurement with respect to the state
     * (may be null if measurement should be ignored)
     * @return innovation covariance matrix, defined as \(h.P.h^T + r\), or
     * null if h is null
     */
    protected RealMatrix computeInnovationCovarianceMatrix(final RealMatrix r, final RealMatrix h) {
        if (h == null) {
            return null;
        }
        final RealMatrix phT = predicted.getCovariance().multiplyTransposed(h);
        return h.multiply(phT).add(r);
    }

    /** Perform correction step.
     * @param measurement single measurement to handle
     * @param stm state transition matrix
     * @param innovation innovation vector (i.e. residuals)
     * (may be null if measurement should be ignored)
     * @param h Jacobian of the measurement with respect to the state
     * (may be null if measurement should be ignored)
     * @param s innovation covariance matrix
     * (may be null if measurement should be ignored)
     * @exception MathIllegalArgumentException if matrix cannot be decomposed
     */
    protected void correct(final T measurement, final RealMatrix stm, final RealVector innovation,
                           final RealMatrix h, final RealMatrix s)
        throws MathIllegalArgumentException {

        if (innovation == null) {
            // measurement should be ignored
            corrected = predicted;
            return;
        }

        // compute Kalman gain k
        // the following is equivalent to k = p.h^T * (h.p.h^T + r)^(-1)
        // we don't want to compute the inverse of a matrix,
        // we start by post-multiplying by h.p.h^T + r and get
        // k.(h.p.h^T + r) = p.h^T
        // then we transpose, knowing that both p and r are symmetric matrices
        // (h.p.h^T + r).k^T = h.p
        // then we can use linear system solving instead of matrix inversion
        final RealMatrix k = decomposer.
                             decompose(s).
                             solve(h.multiply(predicted.getCovariance())).
                             transpose();

        // correct state vector
        final RealVector correctedState = predicted.getState().add(k.operate(innovation));

        // here we use the Joseph algorithm (see "Fundamentals of Astrodynamics and Applications,
        // Vallado, Fourth Edition §10.6 eq.10-34) which is equivalent to
        // the traditional Pest = (I - k.h) x Ppred expression but guarantees the output stays symmetric:
        // Pest = (I -k.h) Ppred (I - k.h)^T + k.r.k^T
        final RealMatrix idMkh = k.multiply(h);
        for (int i = 0; i < idMkh.getRowDimension(); ++i) {
            for (int j = 0; j < idMkh.getColumnDimension(); ++j) {
                idMkh.multiplyEntry(i, j, -1);
            }
            idMkh.addToEntry(i, i, 1.0);
        }
        final RealMatrix r = measurement.getCovariance();
        final RealMatrix correctedCovariance =
                        idMkh.multiply(predicted.getCovariance()).multiplyTransposed(idMkh).
                        add(k.multiply(r).multiplyTransposed(k));

        corrected = new ProcessEstimate(measurement.getTime(), correctedState, correctedCovariance,
                                        stm, h, s, k);

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
