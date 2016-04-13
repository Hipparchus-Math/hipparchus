/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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

package org.hipparchus.ode.nonstiff;


import org.hipparchus.analysis.differentiation.DerivativeStructure;
import org.hipparchus.ode.EquationsMapper;
import org.hipparchus.ode.ExpandableODE;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.sampling.ODEStateInterpolator;
import org.junit.Test;

public class GraggBulirschStoerStateInterpolatorTest extends ODEStateInterpolatorAbstractTest {

    @Test
    public void interpolationAtBounds() {
        doInterpolationAtBounds(1.0e-15);
    }

    @Test
    public void interpolationInside() {
        doInterpolationInside(1.0e-50, 1.0e-50);
    }

    protected ODEStateInterpolator setUpInterpolator(final ReferenceODE eqn,
                                                     final double t0, final double[] y0,
                                                     final double t1) {

        // evaluate derivatives at mid-step
        final int derivationOrder = 7;
        DerivativeStructure middleT = new DerivativeStructure(1, derivationOrder, 0, 0.5 * (t0 + t1));
        DerivativeStructure[] derivatives =  eqn.theoreticalState(middleT);

        double[][] yMidDots = new double[derivationOrder + 1][eqn.getDimension()];
        for (int k = 0; k < yMidDots.length; ++k) {
            for (int i = 0; i < derivatives.length; ++i) {
                yMidDots[k][i] = derivatives[i].getPartialDerivative(k);
            }
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

}
