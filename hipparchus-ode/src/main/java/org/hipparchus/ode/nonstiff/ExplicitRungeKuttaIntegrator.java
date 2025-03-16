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

import org.hipparchus.ode.ExpandableODE;
import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.OrdinaryDifferentialEquation;


/**
 * This interface implements the part of Runge-Kutta
 * integrators for Ordinary Differential Equations
 * common to fixed- and adaptive steps.
 *
 * <p>These methods are explicit Runge-Kutta methods, their Butcher
 * arrays are as follows :</p>
 * <pre>
 *    0  |
 *   c2  | a21
 *   c3  | a31  a32
 *   ... |        ...
 *   cs  | as1  as2  ...  ass-1
 *       |--------------------------
 *       |  b1   b2  ...   bs-1  bs
 * </pre>
 *
 * @see ButcherArrayProvider
 * @see FixedStepRungeKuttaIntegrator
 * @see EmbeddedRungeKuttaIntegrator
 * @since 3.1
 */

public interface ExplicitRungeKuttaIntegrator extends ButcherArrayProvider, ODEIntegrator {

    /**
     * Getter for the number of stages corresponding to the Butcher array.
     *
     * @return number of stages
     */
    default int getNumberOfStages() {
        return getB().length;
    }

    /** Fast computation of a single step of ODE integration.
     * <p>This method is intended for the limited use case of
     * very fast computation of only one step without using any of the
     * rich features of general integrators that may take some time
     * to set up (i.e. no step handlers, no events handlers, no additional
     * states, no interpolators, no error control, no evaluations count,
     * no sanity checks ...). It handles the strict minimum of computation,
     * so it can be embedded in outer loops.</p>
     * <p>
     * This method is <em>not</em> used at all by the {@link #integrate(ExpandableODE, ODEState, double)}
     * method. It also completely ignores the step set at construction time, and
     * uses only a single step to go from {@code t0} to {@code t}.
     * </p>
     * <p>
     * As this method does not use any of the state-dependent features of the integrator,
     * it should be reasonably thread-safe <em>if and only if</em> the provided differential
     * equations are themselves thread-safe.
     * </p>
     * @param equations differential equations to integrate
     * @param t0 initial time
     * @param y0 initial value of the state vector at t0
     * @param t target time for the integration
     * (can be set to a value smaller than {@code t0} for backward integration)
     * @return state vector at {@code t}
     */
    default double[] singleStep(final OrdinaryDifferentialEquation equations, final double t0, final double[] y0,
                                final double t) {

        // create some internal working arrays
        final int stages       = getNumberOfStages();
        final double[][] yDotK = new double[stages][];

        // first stage
        final double h = t - t0;
        final ExpandableODE expandableODE = new ExpandableODE(equations);
        yDotK[0] = expandableODE.computeDerivatives(t0, y0);

        // next stages
        applyInternalButcherWeights(expandableODE, t0, y0, h, getA(), getC(), yDotK);

        // estimate the state at the end of the step
        return applyExternalButcherWeights(y0, yDotK, h, getB());

    }

    /**
     * Apply internal weights of Butcher array, with corresponding times.
     * @param equations differential equations to integrate
     * @param t0        initial time
     * @param y0        initial value of the state vector at t0
     * @param h         step size
     * @param a         internal weights of Butcher array
     * @param c         times of Butcher array
     * @param yDotK     array where to store result
     */
    static void applyInternalButcherWeights(final ExpandableODE equations, final double t0, final double[] y0,
                                            final double h, final double[][] a, final double[] c,
                                            final double[][] yDotK) {
        // create some internal working arrays
        final int stages = c.length + 1;
        final double[] yTmp = y0.clone();

        for (int k = 1; k < stages; ++k) {

            for (int j = 0; j < y0.length; ++j) {
                double sum = yDotK[0][j] * a[k - 1][0];
                for (int l = 1; l < k; ++l) {
                    sum += yDotK[l][j] * a[k - 1][l];
                }
                yTmp[j] = y0[j] + h * sum;
            }

            yDotK[k] = equations.computeDerivatives(t0 + h * c[k - 1], yTmp);
        }
    }

    /** Apply external weights of Butcher array, assuming internal ones have been applied.
     * @param yDotK output of stages
     * @param y0 initial value of the state vector at t0
     * @param h step size
     * @param b external weights of Butcher array
     * @return state vector
     */
    static double[] applyExternalButcherWeights(final double[] y0, final double[][] yDotK, final double h,
                                                final double[] b) {
        final double[] y = y0.clone();
        final int stages = b.length;
        for (int j = 0; j < y0.length; ++j) {
            double sum = yDotK[0][j] * b[0];
            for (int l = 1; l < stages; ++l) {
                sum += yDotK[l][j] * b[l];
            }
            y[j] += h * sum;
        }
        return y;
    }
}
