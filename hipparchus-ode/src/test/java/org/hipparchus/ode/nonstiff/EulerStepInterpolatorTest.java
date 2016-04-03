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


import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Random;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.DenseOutputModel;
import org.hipparchus.ode.ExpandableODE;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.ode.TestProblem1;
import org.hipparchus.ode.TestProblem3;
import org.hipparchus.ode.sampling.ODEStepHandler;
import org.hipparchus.ode.sampling.StepInterpolatorTestUtils;
import org.junit.Assert;
import org.junit.Test;

public class EulerStepInterpolatorTest {

    @Test
    public void noReset() throws MathIllegalStateException {

        OrdinaryDifferentialEquation equation = new OrdinaryDifferentialEquation() {
            public int getDimension() {
                return 3;
            }
            public double[] computeDerivatives(double t, double[] y) {
                return new double[] { 1.0, 2.0, -2.0 };
            }
        };
        double     t0   = 0;
        double[]   y0   =   { 0.0, 1.0, -2.0 };
        double[]   yDot = equation.computeDerivatives(t0, y0);
        double     t1   = 0.1;
        double[]   y1   =  new double[y0.length];
        for (int i = 0; i < y1.length; ++i) {
            y1[i] = y0[i] + (t1 - t0) * yDot[i];
        }
        EulerStepInterpolator interpolator =
                        new EulerStepInterpolator(true,
                                                  new double[0][],
                                                  new ODEStateAndDerivative(0, y0, yDot),
                                                  new ODEStateAndDerivative(0, y1, yDot),
                                                  new ODEStateAndDerivative(0, y0, yDot),
                                                  new ODEStateAndDerivative(0, y1, yDot),
                                                  new ExpandableODE(equation).getMapper());
        double[] result = interpolator.getInterpolatedState(0.5 * (t0 + t1)).getState();
        for (int i = 0; i < equation.getDimension(); ++i) {
            Assert.assertEquals(0.5 * (y0[i] + y1[i]), result[i], 1.0e-10);
        }

    }

    @Test
    public void interpolationAtBounds() throws MathIllegalStateException {

        OrdinaryDifferentialEquation equation = new OrdinaryDifferentialEquation() {
            public int getDimension() {
                return 3;
            }
            public double[] computeDerivatives(double t, double[] y) {
                return new double[] { 1.0, 2.0, -2.0 };
            }
        };
        double     t0   = 0;
        double[]   y0   =   { 0.0, 1.0, -2.0 };
        double[]   yDot = equation.computeDerivatives(t0, y0);
        double     t1   = 0.1;
        double[]   y1   =  new double[y0.length];
        for (int i = 0; i < y1.length; ++i) {
            y1[i] = y0[i] + (t1 - t0) * yDot[i];
        }
        EulerStepInterpolator interpolator =
                        new EulerStepInterpolator(true,
                                                  new double[0][],
                                                  new ODEStateAndDerivative(0, y0, yDot),
                                                  new ODEStateAndDerivative(0, y1, yDot),
                                                  new ODEStateAndDerivative(0, y0, yDot),
                                                  new ODEStateAndDerivative(0, y1, yDot),
                                                  new ExpandableODE(equation).getMapper());
        double[] previous = interpolator.getPreviousState().getState();
        double[] interpolated = interpolator.getInterpolatedState(interpolator.getPreviousState().getTime()).getState();
        for (int i = 0; i < equation.getDimension(); ++i) {
            Assert.assertEquals(previous[i], interpolated[i], 1.0e-10);
        }

        double[] current = interpolator.getCurrentState().getState();
        interpolated = interpolator.getInterpolatedState(interpolator.getCurrentState().getTime()).getState();
        for (int i = 0; i < equation.getDimension(); ++i) {
            Assert.assertEquals(current[i], interpolated[i], 1.0e-10);
        }

    }

    @Test
    public void derivativesConsistency()
        throws MathIllegalArgumentException, MathIllegalStateException {
        TestProblem3 pb = new TestProblem3();
        double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.001;
        EulerIntegrator integ = new EulerIntegrator(step);
        StepInterpolatorTestUtils.checkDerivativesConsistency(integ, pb, 0.01, 5.1e-12);
    }

    @Test
    public void serialization()
        throws IOException, ClassNotFoundException,
               MathIllegalArgumentException, MathIllegalStateException {

        TestProblem1 pb = new TestProblem1();
        double step = (pb.getFinalTime() - pb.getInitialTime()) * 0.001;
        EulerIntegrator integ = new EulerIntegrator(step);
        integ.addStepHandler(new DenseOutputModel());
        integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());

        ByteArrayOutputStream bos = new ByteArrayOutputStream();
        ObjectOutputStream    oos = new ObjectOutputStream(bos);
        for (ODEStepHandler handler : integ.getStepHandlers()) {
            oos.writeObject(handler);
        }

        ByteArrayInputStream  bis = new ByteArrayInputStream(bos.toByteArray());
        ObjectInputStream     ois = new ObjectInputStream(bis);
        DenseOutputModel cm  = (DenseOutputModel) ois.readObject();

        Random random = new Random(347588535632l);
        double maxError = 0.0;
        for (int i = 0; i < 1000; ++i) {
            double r = random.nextDouble();
            double time = r * pb.getInitialTime() + (1.0 - r) * pb.getFinalTime();
            double[] interpolatedY = cm.getInterpolatedState(time).getState();
            double[] theoreticalY  = pb.computeTheoreticalState(time);
            double dx = interpolatedY[0] - theoreticalY[0];
            double dy = interpolatedY[1] - theoreticalY[1];
            double error = dx * dx + dy * dy;
            if (error > maxError) {
                maxError = error;
            }
        }
        Assert.assertTrue(maxError < 0.001);

    }

}
