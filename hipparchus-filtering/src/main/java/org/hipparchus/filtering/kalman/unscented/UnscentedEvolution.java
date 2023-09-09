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
package org.hipparchus.filtering.kalman.unscented;

import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;

/**
 * Container for {@link UnscentedProcess unscented process} evolution data.
 * @see UnscentedProcess
 * @since 2.2
 */
public class UnscentedEvolution {

    /** Current time. */
    private final double currentTime;

    /** State vectors at current time. */
    private final RealVector[] currentStates;

    /** Process noise matrix. */
    private final RealMatrix processNoiseMatrix;

    /**
     * Constructor.
     * @param currentTime current time
     * @param currentStates state vectors at current time
     * @param processNoiseMatrix process noise matrix
     */
    public UnscentedEvolution(final double currentTime, final RealVector[] currentStates,
                              final RealMatrix processNoiseMatrix) {
        this.currentTime         = currentTime;
        this.currentStates       = currentStates.clone();
        this.processNoiseMatrix  = processNoiseMatrix;
    }

    /** Get current time.
     * @return current time
     */
    public double getCurrentTime() {
        return currentTime;
    }

    /** Get current states.
     * @return current states
     */
    public RealVector[] getCurrentStates() {
        return currentStates.clone();
    }

    /** Get process noise.
     * @return process noise
     */
    public RealMatrix getProcessNoiseMatrix() {
        return processNoiseMatrix;
    }

}
