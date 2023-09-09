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

import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;

/**
 * Holder for process state and covariance.
 * <p>
 * The estimate always contains time, state and covariance. These data are
 * the only ones needed to start a Kalman filter. Once a filter has been
 * started and produces new estimates, these new estimates will always
 * contain a state transition matrix and if the measurement has not been
 * ignored, they will also contain measurement Jacobian, innovation covariance
 * and Kalman gain.
 * </p>
 * @since 1.3
 */
public class ProcessEstimate {

    /** Process time (typically the time or index of a measurement). */
    private final double time;

    /** State vector. */
    private final RealVector state;

    /** State covariance. */
    private final RealMatrix covariance;

    /** State transition matrix, may be null.
     * @since 1.4
     */
    private final RealMatrix stateTransitionMatrix;

    /** Jacobian of the measurement with respect to the state (h matrix), may be null.
     * @since 1.4
     */
    private final RealMatrix measurementJacobian;

    /** Innovation covariance matrix, defined as \(h.P.h^T + r\), may be null.
     * @since 1.4
     */
    private final RealMatrix innovationCovarianceMatrix;

    /** Kalman gain (k matrix), may be null.
     * @since 1.4
     */
    private final RealMatrix kalmanGain;

    /** Simple constructor.
     * <p>
     * This constructor sets state transition matrix, covariance matrix H,
     * innovation covariance matrix and Kalman gain k to null.
     * </p>
     * @param time process time (typically the time or index of a measurement)
     * @param state state vector
     * @param covariance state covariance
     */
    public ProcessEstimate(final double time, final RealVector state, final RealMatrix covariance) {
        this(time, state, covariance, null, null, null, null);
    }

    /** Simple constructor.
     * @param time process time (typically the time or index of a measurement)
     * @param state state vector
     * @param covariance state covariance
     * @param stateTransitionMatrix state transition matrix between previous state and estimated (but not yet corrected) state
     * @param measurementJacobian Jacobian of the measurement with respect to the state
     * @param innovationCovariance innovation covariance matrix, defined as \(h.P.h^T + r\), may be null
     * @param kalmanGain Kalman Gain matrix, may be null
     * @since 1.4
     */
    public ProcessEstimate(final double time, final RealVector state, final RealMatrix covariance,
                           final RealMatrix stateTransitionMatrix, final RealMatrix measurementJacobian,
                           final RealMatrix innovationCovariance, final RealMatrix kalmanGain) {
        this.time                       = time;
        this.state                      = state;
        this.covariance                 = covariance;
        this.stateTransitionMatrix      = stateTransitionMatrix;
        this.measurementJacobian        = measurementJacobian;
        this.innovationCovarianceMatrix = innovationCovariance;
        this.kalmanGain                 = kalmanGain;
    }

    /** Get the process time.
     * @return process time (typically the time or index of a measurement)
     */
    public double getTime() {
        return time;
    }

    /** Get the state vector.
     * @return state vector
     */
    public RealVector getState() {
        return state;
    }

    /** Get the state covariance.
     * @return state covariance
     */
    public RealMatrix getCovariance() {
        return covariance;
    }

    /** Get state transition matrix between previous state and estimated (but not yet corrected) state.
     * @return state transition matrix between previous state and estimated state (but not yet corrected)
     * (may be null for initial process estimate)
     * @since 1.4
     */
    public RealMatrix getStateTransitionMatrix() {
        return stateTransitionMatrix;
    }

    /** Get the Jacobian of the measurement with respect to the state (H matrix).
     * @return Jacobian of the measurement with respect to the state (may be null for initial
     * process estimate or if the measurement has been ignored)
     * @since 1.4
     */
    public RealMatrix getMeasurementJacobian() {
        return measurementJacobian;
    }

    /** Get the innovation covariance matrix.
     * @return innovation covariance matrix (may be null for initial
     * process estimate or if the measurement has been ignored)
     * @since 1.4
     */
    public RealMatrix getInnovationCovariance() {
        return innovationCovarianceMatrix;
    }

    /** Get the Kalman gain matrix.
     * @return Kalman gain matrix (may be null for initial
     * process estimate or if the measurement has been ignored)
     * @since 1.4
     */
    public RealMatrix getKalmanGain() {
        return kalmanGain;
    }

}
