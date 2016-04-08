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


import org.hipparchus.ode.EquationsMapper;
import org.hipparchus.ode.ExpandableODE;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

public abstract class RungeKuttaStepInterpolatorAbstractTest {

    protected abstract RungeKuttaStepInterpolator
        createInterpolator(boolean forward, double[][] yDotK,
                           ODEStateAndDerivative globalPreviousState,
                           ODEStateAndDerivative globalCurrentState,
                           ODEStateAndDerivative softPreviousState,
                           ODEStateAndDerivative softCurrentState,
                           EquationsMapper mapper);

    protected abstract ButcherArrayProvider createButcherArrayProvider();

    @Test
    public abstract void interpolationAtBounds();

    protected void doInterpolationAtBounds(double epsilon) {
        RungeKuttaStepInterpolator interpolator = setUpInterpolator(new SinCos(),
                                                                    0.0, new double[] { 0.0, 1.0 }, 0.125);

        Assert.assertEquals(0.0, interpolator.getPreviousState().getTime(), 1.0e-15);
        for (int i = 0; i < 2; ++i) {
            Assert.assertEquals(interpolator.getPreviousState().getState()[i],
                                interpolator.getInterpolatedState(interpolator.getPreviousState().getTime()).getState()[i],
                                epsilon);
        }
        Assert.assertEquals(0.125, interpolator.getCurrentState().getTime(), 1.0e-15);
        for (int i = 0; i < 2; ++i) {
            Assert.assertEquals(interpolator.getCurrentState().getState()[i],
                                interpolator.getInterpolatedState(interpolator.getCurrentState().getTime()).getState()[i],
                                epsilon);
        }
    }

    @Test
    public abstract void interpolationInside();

    protected void doInterpolationInside(double epsilonSin, double epsilonCos) {

        SinCos sinCos = new SinCos();
        RungeKuttaStepInterpolator interpolator = setUpInterpolator(sinCos, 0.0, new double[] { 0.0, 1.0 }, 0.0125);

        int n = 100;
        double maxErrorSin = 0;
        double maxErrorCos = 0;
        for (int i = 0; i <= n; ++i) {
            double t =     ((n - i) * interpolator.getPreviousState().getTime() +
                                 i  * interpolator.getCurrentState().getTime()) / n;
            final double[] interpolated = interpolator.getInterpolatedState(t).getState();
            double[] reference = sinCos.theoreticalState(t);
            maxErrorSin = FastMath.max(maxErrorSin, FastMath.abs(interpolated[0] - reference[0]));
            maxErrorCos = FastMath.max(maxErrorCos, FastMath.abs(interpolated[1] - reference[1]));
        }

        Assert.assertEquals(0.0, maxErrorSin, epsilonSin);
        Assert.assertEquals(0.0, maxErrorCos, epsilonCos);

    }

    private RungeKuttaStepInterpolator setUpInterpolator(final OrdinaryDifferentialEquation eqn,
                                                         final double t0, final double[] y0,
                                                         final double t1) {

        // get the Butcher arrays from the field integrator
        ButcherArrayProvider provider = createButcherArrayProvider();
        double[][] a = provider.getA();
        double[]   b = provider.getB();
        double[]   c = provider.getC();

        // store initial state
        double[][] yDotK = new double[b.length][];
        yDotK[0] = eqn.computeDerivatives(t0, y0);
        ODEStateAndDerivative s0 = new ODEStateAndDerivative(t0, y0, yDotK[0]);

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
        ODEStateAndDerivative s1 = new ODEStateAndDerivative(t1, y, eqn.computeDerivatives(t1, y));

        return createInterpolator(t1 > t0, yDotK, s0, s1, s0, s1,
                                  new ExpandableODE(eqn).getMapper());

    }

    private static class SinCos implements OrdinaryDifferentialEquation {
        public int getDimension() {
            return 2;
        }
        public void init(final double t0, final double[] y0, final double finalTime) {
        }
        public double[] computeDerivatives(final double t, final double[] y) {
            return new double[] {
                y[1], -y[0]
            };
        }
        public double[] theoreticalState(final double t) {
            return new double[] {
                FastMath.sin(t), FastMath.cos(t)
            };
        }
    }

}
