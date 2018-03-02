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

import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;

/**
 * Linear process that can be estimated by a {@link LinearKalmanEstimator}.
 * <p>
 * This interface must be implemented by users to represent the behaviour
 * of the process to be estimated
 * </p>
 * <p>
 * A linear process is governed by the equation:<br>
 *  x<sub>k</sub> = A<sub>k-1</sub> x<sub>k-1</sub> + B<sub>k-1</sub> u<sub>k-1</sub> + w<sub>k-1</sub><br>
 * where <ul>
 *   <li>A<sub>k-1</sub> is the state transition matrix in the absence of control,</li>
 *   <li>B<sub>k-1</sub> is the control matrix,</li>
 *   <li>u<sub>k-1</sub> is the command</li>
 *   <li>w<sub>k-1</sub> is the process noise, which as covariance matrix Q<sub>k-1</sub></li>
 * </ul>
 * </p>
 * @since 1.3
 */
public interface LinearProcess {

    /** Get the state transition matrix A<sub>k-1</sub>.
     * @param time time k-1 of the <em>previous</em> state
     * @return state transition matrix A<sub>k-1</sub> at specified time
     */
    RealMatrix getStateTransitionMatrix(double time);

    /** Get the control matrix B<sub>k-1</sub>.
     * @param time time k-1 of the <em>previous</em> state
     * @return control matrix B<sub>k-1</sub> at specified time
     * (can be null if there is no control)
     */
    RealMatrix getControlMatrix(double time);

    /** Get the command u<sub>k-1</sub>.
     * @param time time k-1 of the <em>previous</em> state
     * @return command vector u<sub>k-1</sub> at specified time
     * (not called if {@link #getControlMatrix(double) getControlMatrix}
     * returns null)
     */
    RealVector getCommand(double time);

    /** Get the process noise matrix Q<sub>k-1</sub>.
     * @param time time k-1 of the <em>previous</em> state
     * @returnprocess noise matrix<sub>k-1</sub> at specified time
     * (can be null if there is no control)
     */
    RealMatrix getProcessNoiseMatrix(double time);

}
