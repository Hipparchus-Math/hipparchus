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

import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;

/**
 * Container for {@link NonLinearProcess non-linear process} evolution data.
 * @see NonLinearProcess
 * @since 1.3
 */
public class NonLinearEvolution {

    /** Current time. */
    private final double currentTime;

    /** State vector at current time. */
    private final RealVector currentState;

    /** State transition matrix between previous and current state. */
    private final RealMatrix stateTransitionMatrix;

    /** Process noise matrix. */
    private final RealMatrix processNoiseMatrix;

    /** Jacobian of the measurement with respect to the state (may be null). */
    private final RealMatrix measurementJacobian;

    /** Simple constructor.
     * @param currentTime current time
     * @param currentState state vector at current time
     * @param stateTransitionMatrix state transition matrix between previous and current state
     * @param processNoiseMatrix process noise
     * @param measurementJacobian Jacobian of the measurement with respect to the state
     * (may be null if measurement should be ignored)
     */
    public NonLinearEvolution(final double currentTime, final RealVector currentState,
                              final RealMatrix stateTransitionMatrix, final RealMatrix processNoiseMatrix,
                              final RealMatrix measurementJacobian) {
        this.currentTime           = currentTime;
        this.currentState          = currentState;
        this.stateTransitionMatrix = stateTransitionMatrix;
        this.processNoiseMatrix    = processNoiseMatrix;
        this.measurementJacobian   = measurementJacobian;
    }

    /** Get current time.
     * @return current time
     */
    public double getCurrentTime() {
        return currentTime;
    }

    /** Get current state.
     * @return current state
     */
    public RealVector getCurrentState() {
        return currentState;
    }

    /** Get state transition matrix between previous and current state.
     * @return state transition matrix between previous and current state
     */
    public RealMatrix getStateTransitionMatrix() {
        return stateTransitionMatrix;
    }

    /** Get process noise.
     * @return process noise
     */
    public RealMatrix getProcessNoiseMatrix() {
        return processNoiseMatrix;
    }

    /** Get measurement Jacobian.
     * @return Jacobian of the measurement with respect to the state
     * (may be null if measurement should be ignored)
     */
    public RealMatrix getMeasurementJacobian() {
        return measurementJacobian;
    }

}
