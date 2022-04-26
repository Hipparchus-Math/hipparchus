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

package org.hipparchus.filtering.kalman.unscented;

import org.hipparchus.filtering.kalman.Measurement;
import org.hipparchus.linear.RealVector;

/**
 * Unscented process that can be estimated by a {@link UnscentedKalmanFilter}.
 * <p>
 * This interface must be implemented by users to represent the behavior
 * of the process to be estimated
 * </p>
 * @param <T> the type of the measurements
 * @see UnscentedKalmanFilter
 * @see org.hipparchus.filtering.kalman.unscented.UnscentedProcess
 */

public interface UnscentedProcess<T extends Measurement>  {

    /** Get the state evolution between two times.
     * @param previousTime time of the previous state
     * @param sigmaPoints sigma points at {@code previousTime}
     * @param measurement measurement to process
     * @return states evolution
     */
    UnscentedEvolution getEvolution(double previousTime, RealVector[] sigmaPoints, T measurement);


}
