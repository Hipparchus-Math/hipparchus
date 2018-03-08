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

package org.hipparchus.filtering;

import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.linear.MatrixDecomposer;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;

/**
 * Kalman estimator for linear process.
 * @since 1.3
 */
public class LinearKalmanEstimator implements KalmanEstimator {

    /** Decomposer decomposer to use for the correction phase. */
    private final MatrixDecomposer decomposer;

    /** Process to be estimated. */
    private final LinearProcess process;

    /** Current estimate. */
    private ProcessEstimate estimate;

    /** Simple constructor.
     * @param decomposer decomposer to use for the correction phase
     * @param process linear process to estimate
     * @param initialState initial state
     */
    public LinearKalmanEstimator(final MatrixDecomposer decomposer,
                                 final LinearProcess process,
                                 final ProcessEstimate initialState) {
        this.decomposer = decomposer;
        this.process    = process;
        this.estimate   = initialState;
    }

    /** {@inheritDoc} */
    @Override
    public ProcessEstimate estimationStep(final Measurement measurement)
        throws MathRuntimeException {

        // prediction phase
        final RealMatrix a  = process.getStateTransitionMatrix(measurement.getTime());
        final RealMatrix b  = process.getControlMatrix(measurement.getTime());
        final RealVector u  = (b == null) ? null : process.getCommand(measurement.getTime());
        final RealMatrix q  = process.getProcessNoiseMatrix(measurement.getTime());

        RealVector predXk = a.operate(estimate.getState());
        if (b != null) {
            predXk = predXk.add(b.operate(u));
        }

        final RealMatrix predCov = a.multiply(estimate.getCovariance().multiply(a.transpose())).add(q);

        // correction phase
        final RealMatrix h          = measurement.getJacobian();
        final RealMatrix r          = measurement.getCovariance();
        final RealMatrix phT        = predCov.multiply(h.transpose());
        final RealVector innovation = measurement.getValue().subtract(h.operate(predXk));

        // the following is equivalent to k = p.h^T * (h.p.h^T + r)^(-1)
        // we don't want to compute the inverse of a matrix,
        // we start by post-multiplying by h.p.h^T +r and get
        // k.(h.p.h^T + r) = p.h^T
        // then we transpose, knowing that both p and r are symmetric matrices
        // (h.p.h^T + r).k^T = h.p
        // and we use linear system solving instead of matrix inversion
        final RealMatrix k = decomposer.
                        decompose(h.multiply(phT).add(r)).
                        solve(h.multiply(predCov)).
                        transpose();
        final RealVector estXk = predXk.add(k.operate(innovation));
        
        // Here we use the Joseph algorithm (see "Fundamentals of Astrodynamics and Applications,
        // Vallado, Fourth Edition §10.6 eq.10-34) which is equivalent to
        // the traditional Pest = (I - k.h) x Ppred expression but guarantees the output stays symmetric:
        // Pest = (I -k.h) Ppred (I - k.h)^T + k.r.k^T
        final RealMatrix idMkh = k.multiply(h);
        for (int i = 0; i < idMkh.getRowDimension(); ++i) {
            for (int j = 0; j < idMkh.getColumnDimension(); ++j) {
                idMkh.setEntry(i, j, -idMkh.getEntry(i, j));
            }
            idMkh.setEntry(i, i, 1.0 + idMkh.getEntry(i, i));
        }
        final RealMatrix estCov =
                        idMkh.multiply(predCov).multiply(idMkh.transpose()).
                        add(k.multiply(r).multiply(k.transpose()));

        estimate = new ProcessEstimate(measurement.getTime(), estXk, estCov);
        return estimate;

    }

}
