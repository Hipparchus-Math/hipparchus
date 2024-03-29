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

package org.hipparchus.ode.sampling;

import org.hipparchus.ode.ODEStateAndDerivative;


/**
 * This interface represents a handler that should be called after
 * each successful step.
 *
 * <p>The ODE integrators compute the evolution of the state vector at
 * some grid points that depend on their own internal algorithm. Once
 * they have found a new grid point (possibly after having computed
 * several evaluation of the derivative at intermediate points), they
 * provide it to objects implementing this interface. These objects
 * typically either ignore the intermediate steps and wait for the
 * last one, store the points in an ephemeris, or forward them to
 * specialized processing or output methods.</p>
 *
 * @see org.hipparchus.ode.ODEIntegrator
 * @see ODEStateInterpolator
 */

public interface ODEStepHandler {

    /** Initialize step handler at the start of an ODE integration.
     * <p>
     * This method is called once at the start of the integration. It
     * may be used by the step handler to initialize some internal data
     * if needed.
     * </p>
     * <p>
     * The default implementation does nothing
     * </p>
     * @param initialState initial time, state vector and derivative
     * @param finalTime target time for the integration
     */
    default void init(ODEStateAndDerivative initialState, double finalTime) {
        // nothing by default
    }

    /**
     * Handle the last accepted step.
     * @param interpolator interpolator for the last accepted step
     */
    void handleStep(ODEStateInterpolator interpolator);

    /**
     * Finalize integration.
     * @param finalState state at integration end
     * @since 2.0
     */
    default void finish(ODEStateAndDerivative finalState) {
        // nothing by default
    }

}
