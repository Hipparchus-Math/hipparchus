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


import org.hipparchus.ode.EquationsMapper;
import org.hipparchus.ode.ExpandableODE;
import org.hipparchus.ode.ODEStateAndDerivative;

public abstract class RungeKuttaStateInterpolatorAbstractTest extends ODEStateInterpolatorAbstractTest {

    protected abstract RungeKuttaStateInterpolator
        createInterpolator(boolean forward, double[][] yDotK,
                           ODEStateAndDerivative globalPreviousState,
                           ODEStateAndDerivative globalCurrentState,
                           ODEStateAndDerivative softPreviousState,
                           ODEStateAndDerivative softCurrentState,
                           EquationsMapper mapper);

    protected abstract ButcherArrayProvider createButcherArrayProvider();

    protected RungeKuttaStateInterpolator setUpInterpolator(final ReferenceODE eqn,
                                                            final double t0, final double[] y0,
                                                            final double t1) {

        // get the Butcher arrays from the field integrator
        ButcherArrayProvider provider = createButcherArrayProvider();
        double[][] a = provider.getA();
        double[]   b = provider.getB();
        double[]   c = provider.getC();

        // store initial state
        EquationsMapper mapper = new ExpandableODE(eqn).getMapper();
        double[][] yDotK = new double[b.length][];
        yDotK[0] = eqn.computeDerivatives(t0, y0);
        ODEStateAndDerivative s0 = mapper.mapStateAndDerivative(t0, y0, yDotK[0]);

        // perform one integration step, in order to get consistent derivatives
        double h = t1 - t0;
        for (int k = 0; k < a.length; ++k) {
            double[] y = y0.clone();
            for (int i = 0; i < y0.length; ++i) {
                for (int s = 0; s <= k; ++s) {
                    y[i] += h * a[k][s] * yDotK[s][i];
                }
            }
            yDotK[k + 1] = eqn.computeDerivatives(t0 + h * c[k], y);
        }

        // store state at step end
        double[] y = y0.clone();
        for (int i = 0; i < y0.length; ++i) {
            for (int s = 0; s < b.length; ++s) {
                y[i] += h * b[s] * yDotK[s][i];
            }
        }
        ODEStateAndDerivative s1 = mapper.mapStateAndDerivative(t1, y, eqn.computeDerivatives(t1, y));

        return createInterpolator(t1 > t0, yDotK, s0, s1, s0, s1, mapper);

    }

}
