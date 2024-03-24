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
package org.hipparchus.random;

import org.hipparchus.util.FastMath;

/** This class is a Gauss-Markov order 1 autoregressive process generator for scalars.
 * @since 3.1
 */

public class GaussMarkovGenerator {

    /** Correlation time. */
    private final double tau;

    /** Standard deviation of the stationary process. */
    private final double stationarySigma;

    /** Underlying generator. */
    private final RandomGenerator generator;

    /** Last generated value. */
    private double last;

    /** Create a new generator.
     * @param tau correlation time
     * @param stationarySigma standard deviation of the stationary process
     * @param generator underlying random generator to use
     */
    public GaussMarkovGenerator(final double tau, final double stationarySigma,
                                final RandomGenerator generator) {
        this.tau             = tau;
        this.stationarySigma = stationarySigma;
        this.generator       = generator;
        this.last            = Double.NaN;
    }

    /** Get the correlation time.
     * @return correlation time
     */
    public double getTau() {
        return tau;
    }

    /** Get the standard deviation of the stationary process.
     * @return standard deviation of the stationary process
     */
    public double getStationarySigma() {
        return stationarySigma;
    }

    /** Generate next step in the autoregressive process.
     * @param deltaT time step since previous estimate (unused at first call)
     * @return a random scalar obeying autoregressive model
     */
    public double next(final double deltaT) {

        if (Double.isNaN(last)) {
            // first generation: use the stationary process
            last = stationarySigma * generator.nextGaussian();
        } else {
            // regular generation: use the autoregressive process
            final double phi    = FastMath.exp(-deltaT / tau);
            final double sigmaE = FastMath.sqrt(1 - phi * phi) * stationarySigma;
            last = phi * last + sigmaE * generator.nextGaussian();
        }

        return last;

    }

}
