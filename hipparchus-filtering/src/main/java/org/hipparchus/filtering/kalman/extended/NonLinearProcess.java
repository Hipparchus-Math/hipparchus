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

package org.hipparchus.filtering.kalman.extended;

import org.hipparchus.filtering.kalman.linear.LinearProcess;
import org.hipparchus.linear.RealVector;

/**
 * Non-linear process that can be estimated by a {@link ExtendedKalmanEstimator}.
 * <p>
 * This interface must be implemented by users to represent the behavior
 * of the process to be estimated
 * </p>
 * <p>
 * A linear process is governed by the equation:<br>
 *  x<sub>k</sub> = f(x<sub>k-1</sub>) + w<sub>k-1</sub><br>
 * where <ul>
 *   <li>f is the non-linear state evolution function,</li>
 *   <li>w<sub>k-1</sub> is the process noise, which as covariance matrix Q<sub>k-1</sub></li>
 * </ul>
 * </p>
 * @see ExtendedKalmanEstimator
 * @see LinearProcess
 * @since 1.3
 */
public interface NonLinearProcess {

    /** Get the state evolution between two times.
     * @param previousTime time of the previous state
     * @param previousState process state at {@code previousTime}
     * @param currentTime time at which 
     * @return state evolution
     */
    NonLinearEvolution getEvolution(double previousTime, RealVector previousState, double currentTime);

}
