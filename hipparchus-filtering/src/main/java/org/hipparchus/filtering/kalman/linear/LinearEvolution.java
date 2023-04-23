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

package org.hipparchus.filtering.kalman.linear;

import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;

/**
 * Container for {@link LinearProcess linear process} evolution data.
 * @see LinearProcess
 * @since 1.3
 */
public class LinearEvolution {

    /** State transition matrix A<sub>k-1</sub>. */
    private final RealMatrix stateTransitionMatrix;

    /** Control matrix B<sub>k-1</sub> (can be null if the process is not controlled). */
    private final RealMatrix controlMatrix;

    /** Command u<sub>k-1</sub>. (can be null if the process is not controlled). */
    private final RealVector command;

    /** Process noise matrix Q<sub>k-1</sub>. */
    private final RealMatrix processNoiseMatrix;

    /** Jacobian of the measurement with respect to the state (may be null). */
    private final RealMatrix measurementJacobian;

    /** Simple constructor.
     * @param stateTransitionMatrix state transition matrix A<sub>k-1</sub>
     * @param controlMatrix control matrix B<sub>k-1</sub> (can be null if the process is not controlled)
     * @param command u<sub>k-1</sub>. (can be null if the process is not controlled)
     * @param processNoiseMatrix process noise matrix Q<sub>k-1</sub>
     * @param measurementJacobian Jacobian of the measurement with respect to the state
     * (may be null if measurement should be ignored)
     */
    public LinearEvolution(final RealMatrix stateTransitionMatrix,
                           final RealMatrix controlMatrix, final RealVector command,
                           final RealMatrix processNoiseMatrix,
                           final RealMatrix measurementJacobian) {
        this.stateTransitionMatrix = stateTransitionMatrix;
        this.controlMatrix         = controlMatrix;
        this.command               = command;
        this.processNoiseMatrix    = processNoiseMatrix;
        this.measurementJacobian   = measurementJacobian;
    }

    /** Get the state transition matrix A<sub>k-1</sub>.
     * @return state transition matrix A<sub>k-1</sub>
     */
    public RealMatrix getStateTransitionMatrix() {
        return stateTransitionMatrix;
    }

    /** Get the control matrix B<sub>k-1</sub>.
     * @return control matrix B<sub>k-1</sub> (can be null if there is no control)
     */
    public RealMatrix getControlMatrix() {
        return controlMatrix;
    }

    /** Get the command u<sub>k-1</sub>.
     * @return command vector u<sub>k-1</sub> (can be null if there is no control)
     */
    public RealVector getCommand() {
        return command;
    }

    /** Get the process noise matrix Q<sub>k-1</sub>.
     * @return process noise matrix<sub>k-1</sub>
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
