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
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.nonstiff.interpolators.HighamHall54StateInterpolator;
import org.hipparchus.util.FastMath;


/**
 * This class implements the 5(4) Higham and Hall integrator for
 * Ordinary Differential Equations.
 *
 * <p>This integrator is an embedded Runge-Kutta integrator
 * of order 5(4) used in local extrapolation mode (i.e. the solution
 * is computed using the high order formula) with stepsize control
 * (and automatic step initialization) and continuous output. This
 * method uses 7 functions evaluations per step.</p>
 *
 */

public class HighamHall54Integrator extends EmbeddedRungeKuttaIntegrator {

    /** Name of integration scheme. */
    public static final String METHOD_NAME = "Higham-Hall 5(4)";

    /** Error weights Butcher array. */
    static final double[] STATIC_E = {
        -1.0/20.0, 0.0, 81.0/160.0, -6.0/5.0, 25.0/32.0, 1.0/16.0, -1.0/10.0
    };

    /** Simple constructor.
     * Build a fifth order Higham and Hall integrator with the given step bounds
     * @param minStep minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param maxStep maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param scalAbsoluteTolerance allowed absolute error
     * @param scalRelativeTolerance allowed relative error
     */
    public HighamHall54Integrator(final double minStep, final double maxStep,
                                  final double scalAbsoluteTolerance,
                                  final double scalRelativeTolerance) {
        super(METHOD_NAME, -1,
              minStep, maxStep, scalAbsoluteTolerance, scalRelativeTolerance);
    }

    /** Simple constructor.
     * Build a fifth order Higham and Hall integrator with the given step bounds
     * @param minStep minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param maxStep maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param vecAbsoluteTolerance allowed absolute error
     * @param vecRelativeTolerance allowed relative error
     */
    public HighamHall54Integrator(final double minStep, final double maxStep,
                                  final double[] vecAbsoluteTolerance,
                                  final double[] vecRelativeTolerance) {
        super(METHOD_NAME, -1,
              minStep, maxStep, vecAbsoluteTolerance, vecRelativeTolerance);
    }

    /** {@inheritDoc} */
    @Override
    public double[] getC() {
        return new double[] {
            2.0/9.0, 1.0/3.0, 1.0/2.0, 3.0/5.0, 1.0, 1.0
        };
    }

    /** {@inheritDoc} */
    @Override
    public double[][] getA() {
        return new double[][] {
            {2.0/9.0},
            {1.0/12.0, 1.0/4.0},
            {1.0/8.0, 0.0, 3.0/8.0},
            {91.0/500.0, -27.0/100.0, 78.0/125.0, 8.0/125.0},
            {-11.0/20.0, 27.0/20.0, 12.0/5.0, -36.0/5.0, 5.0},
            {1.0/12.0, 0.0, 27.0/32.0, -4.0/3.0, 125.0/96.0, 5.0/48.0}
        };
    }

    /** {@inheritDoc} */
    @Override
    public double[] getB() {
        return new double[] {
            1.0/12.0, 0.0, 27.0/32.0, -4.0/3.0, 125.0/96.0, 5.0/48.0, 0.0
        };
    }

    /** {@inheritDoc} */
    @Override
    protected HighamHall54StateInterpolator
    createInterpolator(final boolean forward, double[][] yDotK,
                       final ODEStateAndDerivative globalPreviousState,
                       final ODEStateAndDerivative globalCurrentState,
                       final EquationsMapper mapper) {
        return new HighamHall54StateInterpolator(forward, yDotK,
                                                globalPreviousState, globalCurrentState,
                                                globalPreviousState, globalCurrentState,
                                                mapper);
    }

    /** {@inheritDoc} */
    @Override
    public int getOrder() {
        return 5;
    }

    /** {@inheritDoc} */
    @Override
    protected double estimateError(final double[][] yDotK,
                                   final double[] y0, final double[] y1,
                                   final double h) {

        final StepsizeHelper helper = getStepSizeHelper();
        double error = 0;

        for (int j = 0; j < helper.getMainSetDimension(); ++j) {
            double errSum = STATIC_E[0] * yDotK[0][j];
            for (int l = 1; l < STATIC_E.length; ++l) {
                errSum += STATIC_E[l] * yDotK[l][j];
            }

            final double tol = helper.getTolerance(j, FastMath.max(FastMath.abs(y0[j]), FastMath.abs(y1[j])));
            final double ratio  = h * errSum / tol;
            error += ratio * ratio;

        }

        return FastMath.sqrt(error / helper.getMainSetDimension());

    }

}
