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

package org.hipparchus.ode.nonstiff;


import org.hipparchus.analysis.differentiation.DSFactory;
import org.hipparchus.analysis.differentiation.DerivativeStructure;
import org.hipparchus.ode.EquationsMapper;
import org.hipparchus.ode.ExpandableODE;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.sampling.AbstractODEStateInterpolator;

public class GraggBulirschStoerStateInterpolatorTest extends ODEStateInterpolatorAbstractTest {

    protected AbstractODEStateInterpolator setUpInterpolator(final ReferenceODE eqn,
                                                             final double t0, final double[] y0,
                                                             final double t1) {

        // evaluate scaled derivatives at mid-step
        final int derivationOrder = 7;
        DerivativeStructure middleT = new DSFactory(1, derivationOrder).variable(0, 0.5 * (t0 + t1));
        DerivativeStructure[] derivatives =  eqn.theoreticalState(middleT);

        final double[][] yMidDots = new double[derivationOrder + 1][eqn.getDimension()];
        final double     h        = t1 - t0;
        double           hK       = 1.0;
        for (int k = 0; k < yMidDots.length; ++k) {
            for (int i = 0; i < derivatives.length; ++i) {
                yMidDots[k][i] = hK * derivatives[i].getPartialDerivative(k);
            }
            hK *= h;
        }

        EquationsMapper mapper = new ExpandableODE(eqn).getMapper();
        ODEStateAndDerivative s0 = mapper.mapStateAndDerivative(t0, y0, eqn.computeDerivatives(t0, y0));
        double[] y1 = eqn.theoreticalState(t1);
        ODEStateAndDerivative s1 = mapper.mapStateAndDerivative(t1, y1, eqn.computeDerivatives(t1, y1));

        GraggBulirschStoerStateInterpolator interpolator =
                        new GraggBulirschStoerStateInterpolator(t1 >= t0, s0, s1, s0, s1,
                                                                mapper, yMidDots, derivationOrder);
        return interpolator;

    }

    @Override
    public void interpolationAtBounds() {
        doInterpolationAtBounds(1.0e-15);
    }

    @Override
    public void interpolationInside() {
        doInterpolationInside(3.5e-18, 1.2e-16);
    }

    @Override
    public void restrictPrevious() {
        doRestrictPrevious(1.0e-15, 1.0e-15);
    }

    @Override
    public void restrictCurrent() {
        doRestrictCurrent(1.0e-15, 1.0e-15);
    }

    @Override
    public void restrictBothEnds() {
        doRestrictBothEnds(1.0e-15, 1.0e-15);
    }

    @Override
    public void degenerateInterpolation() {
        doDegenerateInterpolation();
    }

}
