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


import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.ode.TestProblem1;
import org.hipparchus.ode.TestProblem2;
import org.hipparchus.ode.TestProblem3;
import org.hipparchus.ode.TestProblem4;
import org.hipparchus.ode.TestProblem5;
import org.hipparchus.ode.TestProblem6;
import org.hipparchus.ode.TestProblemAbstract;
import org.hipparchus.ode.TestProblemHandler;
import org.hipparchus.ode.events.ODEEventHandler;
import org.hipparchus.ode.sampling.ODEStateInterpolator;
import org.hipparchus.ode.sampling.ODEStepHandler;
import org.hipparchus.ode.sampling.StepInterpolatorTestUtils;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

public class EulerIntegratorTest {

    @Test(expected=MathIllegalArgumentException.class)
    public void testDimensionCheck()
                    throws MathIllegalArgumentException, MathIllegalStateException {
        TestProblem1 pb = new TestProblem1();
        new EulerIntegrator(0.01).integrate(pb,
                                            new ODEState(0.0, new double[pb.getDimension()+10]),
                                            1.0);
    }

    @Test
    public void testDecreasingSteps()
                    throws MathIllegalArgumentException, MathIllegalStateException {

        for (TestProblemAbstract pb : new TestProblemAbstract[] {
            new TestProblem1(), new TestProblem2(), new TestProblem3(),
            new TestProblem4(), new TestProblem5(), new TestProblem6()
        }) {

            double previousValueError = Double.NaN;
            double previousTimeError = Double.NaN;
            for (int i = 4; i < 8; ++i) {

                double step = (pb.getFinalTime() - pb.getInitialTime()) * FastMath.pow(2.0, -i);

                ODEIntegrator integ = new EulerIntegrator(step);
                TestProblemHandler handler = new TestProblemHandler(pb, integ);
                integ.addStepHandler(handler);
                ODEEventHandler[] functions = pb.getEventsHandlers();
                for (int l = 0; l < functions.length; ++l) {
                    integ.addEventHandler(functions[l],
                                          Double.POSITIVE_INFINITY, 1.0e-6 * step, 1000);
                }
                double stopTime = integ.integrate(pb, pb.getInitialState(), pb.getFinalTime()).getTime();
                if (functions.length == 0) {
                    Assert.assertEquals(pb.getFinalTime(), stopTime, 1.0e-10);
                }

                double valueError = handler.getMaximalValueError();
                if (i > 4) {
                    Assert.assertTrue(valueError < FastMath.abs(previousValueError));
                }
                previousValueError = valueError;

                double timeError = handler.getMaximalTimeError();
                if (i > 4) {
                    Assert.assertTrue(timeError <= FastMath.abs(previousTimeError));
                }
                previousTimeError = timeError;

            }

        }

    }

    @Test
    public void testSmallStep()
                    throws MathIllegalArgumentException, MathIllegalStateException {

        TestProblem1 pb  = new TestProblem1();
        double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.001;

        ODEIntegrator integ = new EulerIntegrator(step);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());

        Assert.assertTrue(handler.getLastError() < 2.0e-4);
        Assert.assertTrue(handler.getMaximalValueError() < 1.0e-3);
        Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);
        Assert.assertEquals("Euler", integ.getName());

    }

    @Test
    public void testBigStep()
                    throws MathIllegalArgumentException, MathIllegalStateException {

        TestProblem1 pb  = new TestProblem1();
        double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.2;

        ODEIntegrator integ = new EulerIntegrator(step);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());

        Assert.assertTrue(handler.getLastError() > 0.01);
        Assert.assertTrue(handler.getMaximalValueError() > 0.2);
        Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);

    }

    @Test
    public void testBackward()
                    throws MathIllegalArgumentException, MathIllegalStateException {

        TestProblem5 pb = new TestProblem5();
        double step = FastMath.abs(pb.getFinalTime() - pb.getInitialTime()) * 0.001;

        ODEIntegrator integ = new EulerIntegrator(step);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());

        Assert.assertTrue(handler.getLastError() < 0.45);
        Assert.assertTrue(handler.getMaximalValueError() < 0.45);
        Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);
        Assert.assertEquals("Euler", integ.getName());
    }

    @Test
    public void testStepSize()
                    throws MathIllegalArgumentException, MathIllegalStateException {
        final double step = 1.23456;
        ODEIntegrator integ = new EulerIntegrator(step);
        integ.addStepHandler(new ODEStepHandler() {
            public void handleStep(ODEStateInterpolator interpolator, boolean isLast) {
                if (! isLast) {
                    Assert.assertEquals(step,
                                        interpolator.getCurrentState().getTime() -
                                        interpolator.getPreviousState().getTime(),
                                        1.0e-12);
                }
            }
        });
        integ.integrate(new OrdinaryDifferentialEquation() {
            public double[] computeDerivatives(double t, double[] y) {
                return new double[] { 1.0 };
            }
            public int getDimension() {
                return 1;
            }
        }, new ODEState(0.0, new double[] { 0.0 }), 5.0);
    }

    @Test
    public void derivativesConsistency()
        throws MathIllegalArgumentException, MathIllegalStateException {
        TestProblem3 pb = new TestProblem3();
        double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.001;
        EulerIntegrator integ = new EulerIntegrator(step);
        StepInterpolatorTestUtils.checkDerivativesConsistency(integ, pb, 0.01, 5.1e-12);
    }

}
