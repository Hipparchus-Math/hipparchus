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


import org.hipparchus.CalculusFieldElement;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.FieldExpandableODE;
import org.hipparchus.ode.FieldODEIntegrator;
import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.TestFieldProblemAbstract;
import org.hipparchus.ode.TestProblemAbstract;
import org.hipparchus.util.FastMath;
import org.junit.Assert;

public class StepInterpolatorTestUtils {

    public static void checkDerivativesConsistency(final ODEIntegrator integrator,
                                                   final TestProblemAbstract problem,
                                                   final double finiteDifferencesRatio,
                                                   final double threshold)
        throws MathIllegalArgumentException, MathIllegalStateException {
        integrator.addStepHandler(new ODEStepHandler() {

            public void handleStep(ODEStateInterpolator interpolator) {

                final double dt = interpolator.getCurrentState().getTime() - interpolator.getPreviousState().getTime();
                final double h  = finiteDifferencesRatio * dt;
                final double t  = interpolator.getCurrentState().getTime() - 0.3 * dt;

                if (FastMath.abs(h) < 10 * FastMath.ulp(t)) {
                    return;
                }

                final double[] yM4h = interpolator.getInterpolatedState(t - 4 * h).getPrimaryState();
                final double[] yM3h = interpolator.getInterpolatedState(t - 3 * h).getPrimaryState();
                final double[] yM2h = interpolator.getInterpolatedState(t - 2 * h).getPrimaryState();
                final double[] yM1h = interpolator.getInterpolatedState(t - h).getPrimaryState();
                final double[] yP1h = interpolator.getInterpolatedState(t + h).getPrimaryState();
                final double[] yP2h = interpolator.getInterpolatedState(t + 2 * h).getPrimaryState();
                final double[] yP3h = interpolator.getInterpolatedState(t + 3 * h).getPrimaryState();
                final double[] yP4h = interpolator.getInterpolatedState(t + 4 * h).getPrimaryState();

                final double[] yDot = interpolator.getInterpolatedState(t).getPrimaryDerivative();

                for (int i = 0; i < yDot.length; ++i) {
                    final double approYDot = ( -3 * (yP4h[i] - yM4h[i]) +
                                               32 * (yP3h[i] - yM3h[i]) +
                                             -168 * (yP2h[i] - yM2h[i]) +
                                              672 * (yP1h[i] - yM1h[i])) / (840 * h);
                    Assert.assertEquals("" + (approYDot - yDot[i]), approYDot, yDot[i], threshold);
                }

            }

        });

        integrator.integrate(problem, problem.getInitialState(), problem.getFinalTime());

    }

    public static <T extends CalculusFieldElement<T>> void checkDerivativesConsistency(final FieldODEIntegrator<T> integrator,
                                                                                       final TestFieldProblemAbstract<T> problem,
                                                                                       final double threshold) {
        integrator.addStepHandler(new FieldODEStepHandler<T>() {

            public void handleStep(FieldODEStateInterpolator<T> interpolator) {

                final T h = interpolator.getCurrentState().getTime().subtract(interpolator.getPreviousState().getTime()).multiply(0.001);
                final T t = interpolator.getCurrentState().getTime().subtract(h.multiply(300));

                if (h.abs().subtract(FastMath.ulp(t.getReal()) * 10).getReal() < 0) {
                    return;
                }

                final T[] yM4h = interpolator.getInterpolatedState(t.add(h.multiply(-4))).getPrimaryState();
                final T[] yM3h = interpolator.getInterpolatedState(t.add(h.multiply(-3))).getPrimaryState();
                final T[] yM2h = interpolator.getInterpolatedState(t.add(h.multiply(-2))).getPrimaryState();
                final T[] yM1h = interpolator.getInterpolatedState(t.add(h.multiply(-1))).getPrimaryState();
                final T[] yP1h = interpolator.getInterpolatedState(t.add(h.multiply( 1))).getPrimaryState();
                final T[] yP2h = interpolator.getInterpolatedState(t.add(h.multiply( 2))).getPrimaryState();
                final T[] yP3h = interpolator.getInterpolatedState(t.add(h.multiply( 3))).getPrimaryState();
                final T[] yP4h = interpolator.getInterpolatedState(t.add(h.multiply( 4))).getPrimaryState();

                final T[] yDot = interpolator.getInterpolatedState(t).getPrimaryDerivative();

                for (int i = 0; i < yDot.length; ++i) {
                    final T approYDot =     yP4h[i].subtract(yM4h[i]).multiply(  -3).
                                        add(yP3h[i].subtract(yM3h[i]).multiply(  32)).
                                        add(yP2h[i].subtract(yM2h[i]).multiply(-168)).
                                        add(yP1h[i].subtract(yM1h[i]).multiply( 672)).
                                        divide(h.multiply(840));
                    Assert.assertEquals(approYDot.getReal(), yDot[i].getReal(), threshold);
                }

            }

        });

        integrator.integrate(new FieldExpandableODE<T>(problem), problem.getInitialState(), problem.getFinalTime());

    }
}

