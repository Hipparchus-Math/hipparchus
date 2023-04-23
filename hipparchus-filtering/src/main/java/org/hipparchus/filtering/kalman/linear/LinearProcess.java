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

import org.hipparchus.filtering.kalman.Measurement;

/**
 * Linear process that can be estimated by a {@link LinearKalmanFilter}.
 * <p>
 * This interface must be implemented by users to represent the behavior
 * of the process to be estimated
 * </p>
 * <p>
 * A linear process is governed by the equation:
 * \(
 *  x_k = A_{k-1} x_{k-1} + B_{k-1} u_{k-1} + w_{k-1}
 * \)
 * where</p>
 * <ul>
 *   <li>A<sub>k-1</sub> is the state transition matrix in the absence of control,</li>
 *   <li>B<sub>k-1</sub> is the control matrix,</li>
 *   <li>u<sub>k-1</sub> is the command</li>
 *   <li>w<sub>k-1</sub> is the process noise, which has covariance matrix Q<sub>k-1</sub></li>
 * </ul>
 * @param <T> the type of the measurements
 * @see LinearKalmanFilter
 * @see org.hipparchus.filtering.kalman.extended.NonLinearProcess
 * @since 1.3
 */
public interface LinearProcess<T extends Measurement> {

    /** Get the state evolution between two times.
     * @param measurement measurement to process
     * @return state evolution
     */
    LinearEvolution getEvolution(T measurement);

}
