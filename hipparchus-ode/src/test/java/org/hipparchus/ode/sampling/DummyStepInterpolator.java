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

package org.hipparchus.ode.sampling;

import org.hipparchus.ode.EquationsMapper;
import org.hipparchus.ode.ODEStateAndDerivative;

/** This class is a step interpolator that does nothing.
 *
 * <p>This class is used when the {@link ODEStepHandler "step handler"}
 * set up by the user does not need step interpolation. It does not
 * recompute the state when {@link AbstractODEStateInterpolator#getInterpolatedState(double)
 * getInterpolatedState} is called, only updating time and copying current state
 * vector.</p>
 *
 * @see ODEStepHandler
 *
 */

public class DummyStepInterpolator extends AbstractODEStateInterpolator {

    /** Serializable version identifier. */
    private static final long serialVersionUID = 20160402L;

    /** Simple constructor.
     * @param isForward integration direction indicator
     * @param globalPreviousState start of the global step
     * @param globalCurrentState end of the global step
     * @param softPreviousState start of the restricted step
     * @param softCurrentState end of the restricted step
     * @param equationsMapper mapper for ODE equations primary and secondary components
     */
    public DummyStepInterpolator(final boolean isForward,
                                 final ODEStateAndDerivative globalPreviousState,
                                 final ODEStateAndDerivative globalCurrentState,
                                 final ODEStateAndDerivative softPreviousState,
                                 final ODEStateAndDerivative softCurrentState,
                                 final EquationsMapper equationsMapper) {
        super(isForward, globalPreviousState, globalCurrentState,
              softPreviousState, softCurrentState, equationsMapper);
    }

    /** {@inheritDoc} */
    @Override
    protected DummyStepInterpolator create(boolean newForward,
                                           ODEStateAndDerivative newGlobalPreviousState,
                                           ODEStateAndDerivative newGlobalCurrentState,
                                           ODEStateAndDerivative newSoftPreviousState,
                                           ODEStateAndDerivative newSoftCurrentState,
                                           EquationsMapper newMapper) {
        return new DummyStepInterpolator(newForward,
                                         newGlobalPreviousState, newGlobalCurrentState,
                                         newSoftPreviousState, newSoftCurrentState,
                                         newMapper);
    }

    /** {@inheritDoc}.
     * In this class, this method does nothing: the interpolated state
     * is always the state at the end of the current step.
     */
    @Override
    protected ODEStateAndDerivative computeInterpolatedStateAndDerivatives(EquationsMapper equationsMapper,
                                                                           double time, double theta,
                                                                           double thetaH, double oneMinusThetaH) {
        return new ODEStateAndDerivative(time,
                                         getCurrentState().getCompleteState(),
                                         getCurrentState().getCompleteDerivative());
    }

}
