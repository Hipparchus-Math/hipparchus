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
 * Holder for a measurement on process.
 * @since 1.3
 */
public class Measurement {

    /** Measurement time or index. */
    private final double time;

    /** Measurement vector. */
    private final RealVector value;

    /** Measurement Jacobian matrix (partial derivatives of measurement with respect to process state). */
    private final RealMatrix jacobian;

    /** Measurement covariance. */
    private final RealMatrix covariance;

    /** Simple constructor.
     * @param time measurement time or index
     * @param value measurement vector
     * @param jacobian measurement Jacobian matrix (partial derivatives of measurement with respect to process state)
     * @param covariance measurement covariance
     */
    public Measurement(final double time, final RealVector value, final RealMatrix jacobian,
                       final RealMatrix covariance) {
        this.time       = time;
        this.value      = value;
        this.jacobian   = jacobian;
        this.covariance = covariance;
    }

    /** Get the process time.
     * @return process time (typically the time or index of a measurement)
     */
    public double getTime() {
        return time;
    }

    /** Get the measurement vector.
     * @return measurement vector
     */
    public RealVector getValue() {
        return value;
    }

    /** Get the measurement Jacobian matrix.
     * <p>
     * The jacobian matrix contains the partial derivatives of measurement
     * with respect to process state.
     * </p>
     * @return measurement Jacobian matrix
     */
    public RealMatrix getJacobian() {
        return jacobian;
    }

    /** Get the measurement covariance.
     * @return measurement covariance
     */
    public RealMatrix getCovariance() {
        return covariance;
    }

}
