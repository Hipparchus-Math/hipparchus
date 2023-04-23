/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */

package org.hipparchus.ode;

import org.hipparchus.ode.events.ODEEventDetector;

/**
 * This class is used as the base class of the problems that are
 * integrated during the junit tests for the ODE integrators.
 */
public abstract class TestProblemAbstract
    implements OrdinaryDifferentialEquation {

    /** Number of functions calls. */
    private int calls;

    /** Initial state */
    private final ODEState s0;

    /** Final time */
    private final double t1;

    /** Error scale */
    private final double[] errorScale;

    /**
     * Simple constructor.
     * @param t0 initial time
     * @param y0 initial state vector
     * @param t1 final time
     * @param errorScale error scale
     */
    protected TestProblemAbstract(double t0, double[] y0, double t1, double[] errorScale) {
        calls           = 0;
        s0              = new ODEState(t0, y0);
        this.t1         = t1;
        this.errorScale = errorScale.clone();
    }

    public int getDimension() {
        return s0.getPrimaryState().length;
    }

    /**
     * Get the initial time.
     * @return initial time
     */
    public double getInitialTime() {
        return s0.getTime();
    }

    /**
     * Get the initial state vector.
     * @return initial state vector
     */
    public ODEState getInitialState() {
        return s0;
    }

    /**
     * Get the final time.
     * @return final time
     */
    public double getFinalTime() {
        return t1;
    }

    /**
     * Get the error scale.
     * @return error scale
     */
    public double[] getErrorScale() {
        return errorScale;
    }

    /**
     * Get the event detectors.
     * @param maxCheck maximum checking interval, must be strictly positive (s)
     * @param threshold convergence threshold (s)
     * @param maxIter maximum number of iterations in the event time search
     * @return events detectors   */
    public ODEEventDetector[] getEventDetectors(final double maxCheck, final double threshold, final int maxIter) {
        return new ODEEventDetector[0];
    }

    /**
     * Get the theoretical events times.
     * @return theoretical events times
     */
    public double[] getTheoreticalEventsTimes() {
        return new double[0];
    }

    /**
     * Get the number of calls.
     * @return nuber of calls
     */
    public int getCalls() {
        return calls;
    }

    public double[] computeDerivatives(double t, double[] y) {
        ++calls;
        return doComputeDerivatives(t, y);
    }

    abstract public double[] doComputeDerivatives(double t, double[] y);

    /**
     * Compute the theoretical state at the specified time.
     * @param t time at which the state is required
     * @return state vector at time t
     */
    abstract public double[] computeTheoreticalState(double t);

}
