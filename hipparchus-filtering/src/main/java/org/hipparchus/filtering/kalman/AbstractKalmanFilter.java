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

package org.hipparchus.filtering.kalman;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.linear.MatrixDecomposer;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;

/**
 * Shared parts between linear and non-linear Kalman filters.
 * @since 1.3
 */
public abstract class AbstractKalmanFilter implements KalmanFilter {

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
    }

    /** Perform correction step.
     * @param measurement single measurement to handle
     * @exception MathIllegalArgumentException if matrix cannot be decomposed
     */
    protected void correct(final Measurement measurement)
        throws MathIllegalArgumentException {

        // correction phase
        final RealMatrix h          = measurement.getJacobian();
        final RealMatrix r          = measurement.getCovariance();
        final RealMatrix phT        = predicted.getCovariance().multiplyTransposed(h);
        final RealVector innovation = measurement.getValue().subtract(h.operate(predicted.getState()));

        // compute Kalman gain k
        // the following is equivalent to k = p.h^T * (h.p.h^T + r)^(-1)
        // we don't want to compute the inverse of a matrix,
        // we start by post-multiplying by h.p.h^T + r and get
        // k.(h.p.h^T + r) = p.h^T
        // then we transpose, knowing that both p and r are symmetric matrices
        // (h.p.h^T + r).k^T = h.p
        // then we can use linear system solving instead of matrix inversion
        final RealMatrix k = decomposer.
                        decompose(h.multiply(phT).add(r)).
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
        final RealMatrix correctedCovariance =
                        idMkh.multiply(predicted.getCovariance()).multiplyTransposed(idMkh).
                        add(k.multiply(r).multiplyTransposed(k));

        corrected = new ProcessEstimate(measurement.getTime(), correctedState, correctedCovariance);

    }

    /** Get the corrected state.
     * @return corrected state
     */
    protected ProcessEstimate getCorrected() {
        return corrected;
    }

}
