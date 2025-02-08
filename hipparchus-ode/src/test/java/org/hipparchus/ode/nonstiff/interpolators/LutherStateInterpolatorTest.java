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
import org.hipparchus.ode.nonstiff.ButcherArrayProvider;
import org.hipparchus.ode.nonstiff.LutherIntegrator;
import org.junit.jupiter.api.Test;

class LutherStateInterpolatorTest extends RungeKuttaStateInterpolatorAbstractTest {

    @Override
    protected RungeKuttaStateInterpolator
        createInterpolator(boolean forward, double[][] yDotK,
                           ODEStateAndDerivative globalPreviousState,
                           ODEStateAndDerivative globalCurrentState,
                           ODEStateAndDerivative softPreviousState,
                           ODEStateAndDerivative softCurrentState,
                           EquationsMapper mapper) {
        return new LutherStateInterpolator(forward, yDotK,
                                          globalPreviousState, globalCurrentState,
                                          softPreviousState, softCurrentState,
                                          mapper);
    }

    @Override
    protected ButcherArrayProvider createButcherArrayProvider() {
        return new LutherIntegrator(1.0);
    }

    @Test
    public void interpolationAtBounds() {
        doInterpolationAtBounds(1.0e-15);
    }

    @Test
    public void interpolationInside() {
        doInterpolationInside(1.1e-7, 9.6e-9);
    }

    @Override
    @Test
    public void restrictPrevious() {
        doRestrictPrevious(1.0e-15, 1.0e-15);
    }

    @Override
    @Test
    public void restrictCurrent() {
        doRestrictCurrent(1.0e-15, 1.0e-15);
    }

    @Override
    @Test
    public void restrictBothEnds() {
        doRestrictBothEnds(1.0e-15, 1.0e-15);
    }

    @Override
    @Test
    public void degenerateInterpolation() {
        doDegenerateInterpolation();
    }

}
