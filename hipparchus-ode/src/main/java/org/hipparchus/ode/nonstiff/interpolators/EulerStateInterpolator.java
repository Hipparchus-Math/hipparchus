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

package org.hipparchus.ode.nonstiff.interpolators;

import org.hipparchus.ode.EquationsMapper;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.nonstiff.EulerIntegrator;

/**
 * This class implements a linear interpolator for step.
 *
 * <p>This interpolator computes dense output inside the last
 * step computed. The interpolation equation is consistent with the
 * integration scheme :</p>
 * <ul>
 *   <li>Using reference point at step start:<br>
 *     y(t<sub>n</sub> + &theta; h) = y (t<sub>n</sub>) + &theta; h y'
 *   </li>
 *   <li>Using reference point at step end:<br>
 *     y(t<sub>n</sub> + &theta; h) = y (t<sub>n</sub> + h) - (1-&theta;) h y'
 *   </li>
 * </ul>
 *
 * <p>where &theta; belongs to [0 ; 1] and where y' is the evaluation of
 * the derivatives already computed during the step.</p>
 *
 * @see EulerIntegrator
 */

public class EulerStateInterpolator extends RungeKuttaStateInterpolator {

    /** Serializable version identifier. */
    private static final long serialVersionUID = 20160328L;

    /** Simple constructor.
     * @param forward integration direction indicator
     * @param yDotK slopes at the intermediate points
     * @param globalPreviousState start of the global step
     * @param globalCurrentState end of the global step
     * @param softPreviousState start of the restricted step
     * @param softCurrentState end of the restricted step
     * @param mapper equations mapper for the all equations
     */
    public EulerStateInterpolator(final boolean forward,
                                  final double[][] yDotK,
                                  final ODEStateAndDerivative globalPreviousState,
                                  final ODEStateAndDerivative globalCurrentState,
                                  final ODEStateAndDerivative softPreviousState,
                                  final ODEStateAndDerivative softCurrentState,
                                  final EquationsMapper mapper) {
        super(forward, yDotK, globalPreviousState, globalCurrentState, softPreviousState, softCurrentState, mapper);
    }

    /** {@inheritDoc} */
    @Override
    protected EulerStateInterpolator create(final boolean newForward, final double[][] newYDotK,
                                            final ODEStateAndDerivative newGlobalPreviousState,
                                            final ODEStateAndDerivative newGlobalCurrentState,
                                            final ODEStateAndDerivative newSoftPreviousState,
                                            final ODEStateAndDerivative newSoftCurrentState,
                                            final EquationsMapper newMapper) {
        return new EulerStateInterpolator(newForward, newYDotK,
                                          newGlobalPreviousState, newGlobalCurrentState,
                                          newSoftPreviousState, newSoftCurrentState,
                                          newMapper);
    }

    /** {@inheritDoc} */
    @Override
    protected ODEStateAndDerivative computeInterpolatedStateAndDerivatives(final EquationsMapper mapper,
                                                                           final double time, final double theta,
                                                                           final double thetaH, final double oneMinusThetaH) {
        final double[] interpolatedState;
        if (getGlobalPreviousState() != null && theta <= 0.5) {
            interpolatedState       = previousStateLinearCombination(thetaH);
        } else {
            interpolatedState       = currentStateLinearCombination(-oneMinusThetaH);
        }
        final double[] interpolatedDerivatives = derivativeLinearCombination(1.0);

        return mapper.mapStateAndDerivative(time, interpolatedState, interpolatedDerivatives);

    }

}
