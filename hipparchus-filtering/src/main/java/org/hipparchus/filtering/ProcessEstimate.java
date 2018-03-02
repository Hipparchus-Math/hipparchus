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
 * Holder for process state and covariance.
 * @since 1.3
 */
public class ProcessEstimate {

    /** Process time (typically the time or index of a measurement). */
    private final double time;

    /** State vector. */
    private final RealVector state;

    /** State covariance. */
    private final RealMatrix covariance;

    /** Simple constructor.
     * @param time process time (typically the time or index of a measurement)
     * @param state state vector
     * @param covariance state covariance
     */
    public ProcessEstimate(final double time, final RealVector state, final RealMatrix covariance) {
        this.time       = time;
        this.state      = state;
        this.covariance = covariance;
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

}
