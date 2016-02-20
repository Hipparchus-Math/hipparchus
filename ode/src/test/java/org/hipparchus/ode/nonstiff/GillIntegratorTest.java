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
import org.hipparchus.ode.FirstOrderDifferentialEquations;
import org.hipparchus.ode.FirstOrderIntegrator;
import org.hipparchus.ode.TestProblem1;
import org.hipparchus.ode.TestProblem2;
import org.hipparchus.ode.TestProblem3;
import org.hipparchus.ode.TestProblem4;
import org.hipparchus.ode.TestProblem5;
import org.hipparchus.ode.TestProblem6;
import org.hipparchus.ode.TestProblemAbstract;
import org.hipparchus.ode.TestProblemHandler;
import org.hipparchus.ode.events.EventHandler;
import org.hipparchus.ode.sampling.StepHandler;
import org.hipparchus.ode.sampling.StepInterpolator;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

public class GillIntegratorTest {

  @Test(expected=MathIllegalArgumentException.class)
  public void testDimensionCheck()
      throws MathIllegalArgumentException,
             MathIllegalArgumentException, MathIllegalStateException {
      TestProblem1 pb = new TestProblem1();
      new GillIntegrator(0.01).integrate(pb,
                                         0.0, new double[pb.getDimension()+10],
                                         1.0, new double[pb.getDimension()+10]);
        Assert.fail("an exception should have been thrown");
  }

  @Test
  public void testDecreasingSteps()
      throws MathIllegalArgumentException,
             MathIllegalArgumentException, MathIllegalStateException {

      for (TestProblemAbstract pb : new TestProblemAbstract[] {
          new TestProblem1(), new TestProblem2(), new TestProblem3(),
          new TestProblem4(), new TestProblem5(), new TestProblem6()
      }) {

      double previousValueError = Double.NaN;
      double previousTimeError = Double.NaN;
      for (int i = 5; i < 10; ++i) {

        double step = (pb.getFinalTime() - pb.getInitialTime()) * FastMath.pow(2.0, -i);

        FirstOrderIntegrator integ = new GillIntegrator(step);
        TestProblemHandler handler = new TestProblemHandler(pb, integ);
        integ.addStepHandler(handler);
        EventHandler[] functions = pb.getEventsHandlers();
        for (int l = 0; l < functions.length; ++l) {
          integ.addEventHandler(functions[l],
                                     Double.POSITIVE_INFINITY, 1.0e-6 * step, 1000);
        }
        double stopTime = integ.integrate(pb, pb.getInitialTime(), pb.getInitialState(),
                                          pb.getFinalTime(), new double[pb.getDimension()]);
        if (functions.length == 0) {
            Assert.assertEquals(pb.getFinalTime(), stopTime, 1.0e-10);
        }

        double valueError = handler.getMaximalValueError();
        if (i > 5) {
          Assert.assertTrue(valueError < 1.01 * FastMath.abs(previousValueError));
        }
        previousValueError = valueError;

        double timeError = handler.getMaximalTimeError();
        if (i > 5) {
          Assert.assertTrue(timeError <= FastMath.abs(previousTimeError));
        }
        previousTimeError = timeError;

      }

    }

  }

  @Test
  public void testSmallStep()
      throws MathIllegalArgumentException,
             MathIllegalArgumentException, MathIllegalStateException {

    TestProblem1 pb = new TestProblem1();
    double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.001;

    FirstOrderIntegrator integ = new GillIntegrator(step);
    TestProblemHandler handler = new TestProblemHandler(pb, integ);
    integ.addStepHandler(handler);
    integ.integrate(pb, pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);

    Assert.assertTrue(handler.getLastError() < 2.0e-13);
    Assert.assertTrue(handler.getMaximalValueError() < 4.0e-12);
    Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);
    Assert.assertEquals("Gill", integ.getName());

  }

  @Test
  public void testBigStep()
      throws MathIllegalArgumentException,
             MathIllegalArgumentException, MathIllegalStateException {

    TestProblem1 pb = new TestProblem1();
    double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.2;

    FirstOrderIntegrator integ = new GillIntegrator(step);
    TestProblemHandler handler = new TestProblemHandler(pb, integ);
    integ.addStepHandler(handler);
    integ.integrate(pb, pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);

    Assert.assertTrue(handler.getLastError() > 0.0004);
    Assert.assertTrue(handler.getMaximalValueError() > 0.005);
    Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);

  }

  @Test
  public void testBackward()
      throws MathIllegalArgumentException,
             MathIllegalArgumentException, MathIllegalStateException {

      TestProblem5 pb = new TestProblem5();
      double step = FastMath.abs(pb.getFinalTime() - pb.getInitialTime()) * 0.001;

      FirstOrderIntegrator integ = new GillIntegrator(step);
      TestProblemHandler handler = new TestProblemHandler(pb, integ);
      integ.addStepHandler(handler);
      integ.integrate(pb, pb.getInitialTime(), pb.getInitialState(),
                      pb.getFinalTime(), new double[pb.getDimension()]);

      Assert.assertTrue(handler.getLastError() < 5.0e-10);
      Assert.assertTrue(handler.getMaximalValueError() < 7.0e-10);
      Assert.assertEquals(0, handler.getMaximalTimeError(), 1.0e-12);
      Assert.assertEquals("Gill", integ.getName());
  }

  @Test
  public void testKepler()
      throws MathIllegalArgumentException,
             MathIllegalArgumentException, MathIllegalStateException {

    final TestProblem3 pb  = new TestProblem3(0.9);
    double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.0003;

    FirstOrderIntegrator integ = new GillIntegrator(step);
    integ.addStepHandler(new KeplerStepHandler(pb));
    integ.integrate(pb,
                    pb.getInitialTime(), pb.getInitialState(),
                    pb.getFinalTime(), new double[pb.getDimension()]);
  }

  @Test
  public void testUnstableDerivative()
      throws MathIllegalArgumentException,
             MathIllegalArgumentException, MathIllegalStateException {
    final StepProblem stepProblem = new StepProblem(0.0, 1.0, 2.0);
    FirstOrderIntegrator integ = new GillIntegrator(0.3);
    integ.addEventHandler(stepProblem, 1.0, 1.0e-12, 1000);
    double[] y = { Double.NaN };
    integ.integrate(stepProblem, 0.0, new double[] { 0.0 }, 10.0, y);
    Assert.assertEquals(8.0, y[0], 1.0e-12);
  }

  private static class KeplerStepHandler implements StepHandler {
    public KeplerStepHandler(TestProblem3 pb) {
      this.pb = pb;
    }
    public void init(double t0, double[] y0, double t) {
      maxError = 0;
    }
    public void handleStep(StepInterpolator interpolator, boolean isLast)
        throws MathIllegalStateException {

      double[] interpolatedY = interpolator.getInterpolatedState();
      double[] theoreticalY  = pb.computeTheoreticalState(interpolator.getCurrentTime());
      double dx = interpolatedY[0] - theoreticalY[0];
      double dy = interpolatedY[1] - theoreticalY[1];
      double error = dx * dx + dy * dy;
      if (error > maxError) {
        maxError = error;
      }
      if (isLast) {
        // even with more than 1000 evaluations per period,
        // Gill is not able to integrate such an eccentric
        // orbit with a good accuracy
        Assert.assertTrue(maxError > 0.001);
      }
    }
    private double maxError;
    private TestProblem3 pb;
  }

  @Test
  public void testStepSize()
      throws MathIllegalArgumentException,
             MathIllegalArgumentException, MathIllegalStateException {
      final double step = 1.23456;
      FirstOrderIntegrator integ = new GillIntegrator(step);
      integ.addStepHandler(new StepHandler() {
          public void handleStep(StepInterpolator interpolator, boolean isLast) {
              if (! isLast) {
                  Assert.assertEquals(step,
                               interpolator.getCurrentTime() - interpolator.getPreviousTime(),
                               1.0e-12);
              }
          }
          public void init(double t0, double[] y0, double t) {
          }
      });
      integ.integrate(new FirstOrderDifferentialEquations() {
          public void computeDerivatives(double t, double[] y, double[] dot) {
              dot[0] = 1.0;
          }
          public int getDimension() {
              return 1;
          }
      }, 0.0, new double[] { 0.0 }, 5.0, new double[1]);
  }

}
