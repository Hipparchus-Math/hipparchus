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

package org.hipparchus.filtering.kalman.extended;

import org.hipparchus.filtering.kalman.Measurement;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;

/**
 * Non-linear process that can be estimated by a {@link ExtendedKalmanFilter}.
 * <p>
 * This interface must be implemented by users to represent the behavior
 * of the process to be estimated
 * </p>
 * @param <T> the type of the measurements
 * @see ExtendedKalmanFilter
 * @see org.hipparchus.filtering.kalman.linear.LinearProcess
 * @since 1.3
 */
public interface NonLinearProcess<T extends Measurement> {

    /** Get the state evolution between two times.
     * @param previousTime time of the previous state
     * @param previousState process state at {@code previousTime}
     * @param measurement measurement to process
     * @return state evolution
     */
    NonLinearEvolution getEvolution(double previousTime, RealVector previousState, T measurement);

    /** Get the innovation brought by a measurement.
     * @param measurement measurement to process
     * @param evolution evolution returned by a previous call to {@link #getEvolution(double, RealVector, Measurement)}
     * @param innovationCovarianceMatrix innovation covariance matrix, defined as \(h.P.h^T + r\)
     * where h is the {@link NonLinearEvolution#getMeasurementJacobian() measurement Jacobian},
     * P is the predicted covariance and r is {@link Measurement#getCovariance() measurement covariance}
     * @return innovation brought by a measurement, may be null if measurement should be rejected
     */
    RealVector getInnovation(T measurement, NonLinearEvolution evolution, RealMatrix innovationCovarianceMatrix);

}
