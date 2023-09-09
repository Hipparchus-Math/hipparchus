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

package org.hipparchus.ode;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Random;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.nonstiff.DormandPrince54Integrator;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.ode.nonstiff.EulerIntegrator;
import org.hipparchus.ode.sampling.DummyStepInterpolator;
import org.hipparchus.ode.sampling.ODEStateInterpolator;
import org.hipparchus.ode.sampling.ODEStepHandler;
import org.hipparchus.util.FastMath;
import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

public class DenseOutputModelTest {

    TestProblem3 pb;
    ODEIntegrator integ;

    @Test
    public void testBoundaries() throws MathIllegalArgumentException, MathIllegalStateException {
        integ.addStepHandler(new DenseOutputModel());
        integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());
        DenseOutputModel dom = (DenseOutputModel) integ.getStepHandlers().iterator().next();
        double tBefore = 2.0 * pb.getInitialTime() - pb.getFinalTime();
        Assert.assertEquals(tBefore, dom.getInterpolatedState(tBefore).getTime(), 1.0e-10);
        double tAfter = 2.0 * pb.getFinalTime() - pb.getInitialTime();
        Assert.assertEquals(tAfter, dom.getInterpolatedState(tAfter).getTime(), 1.0e-10);
        double tMiddle = 2.0 * pb.getFinalTime() - pb.getInitialTime();
        Assert.assertEquals(tMiddle, dom.getInterpolatedState(tMiddle).getTime(), 1.0e-10);
    }

    @Test
    public void testRandomAccess() throws MathIllegalArgumentException, MathIllegalStateException {

        DenseOutputModel dom = new DenseOutputModel();
        integ.addStepHandler(dom);
        integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());

        Random random = new Random(347588535632l);
        double maxError    = 0.0;
        double maxErrorDot = 0.0;
        for (int i = 0; i < 1000; ++i) {
            double r = random.nextDouble();
            double time = r * pb.getInitialTime() + (1.0 - r) * pb.getFinalTime();
            ODEStateAndDerivative sd = dom.getInterpolatedState(time);
            double[] interpolatedY    = sd.getPrimaryState();
            double[] interpolatedYDot = sd.getPrimaryDerivative();
            double[] theoreticalY     = pb.computeTheoreticalState(time);
            double[] theoreticalYDot  = pb.doComputeDerivatives(time, theoreticalY);
            double dx = interpolatedY[0] - theoreticalY[0];
            double dy = interpolatedY[1] - theoreticalY[1];
            double error = dx * dx + dy * dy;
            maxError = FastMath.max(maxError, error);
            double dxDot = interpolatedYDot[0] - theoreticalYDot[0];
            double dyDot = interpolatedYDot[1] - theoreticalYDot[1];
            double errorDot = dxDot * dxDot + dyDot * dyDot;
            maxErrorDot = FastMath.max(maxErrorDot, errorDot);
        }

        Assert.assertEquals(0.0, maxError,    1.0e-9);
        Assert.assertEquals(0.0, maxErrorDot, 4.0e-7);

    }

    @Test
    public void testModelsMerging() throws MathIllegalArgumentException, MathIllegalStateException {

        // theoretical solution: y[0] = cos(t), y[1] = sin(t)
        OrdinaryDifferentialEquation problem =
                        new OrdinaryDifferentialEquation() {
            @Override
            public double[] computeDerivatives(double t, double[] y) {
                return new double[] { -y[1], y[0] };
            }
            @Override
            public int getDimension() {
                return 2;
            }
        };

        // integrate backward from &pi; to 0;
        DenseOutputModel dom1 = new DenseOutputModel();
        ODEIntegrator integ1 = new DormandPrince853Integrator(0, 1.0, 1.0e-8, 1.0e-8);
        integ1.addStepHandler(dom1);
        integ1.integrate(problem, new ODEState(FastMath.PI, new double[] { -1.0, 0.0 }), 0);

        // integrate backward from 2&pi; to &pi;
        DenseOutputModel dom2 = new DenseOutputModel();
        ODEIntegrator integ2 = new DormandPrince853Integrator(0, 0.1, 1.0e-12, 1.0e-12);
        integ2.addStepHandler(dom2);
        integ2.integrate(problem, new ODEState(2.0 * FastMath.PI, new double[] { 1.0, 0.0 }), FastMath.PI);

        // merge the two half circles
        DenseOutputModel dom = new DenseOutputModel();
        dom.append(dom2);
        dom.append(new DenseOutputModel());
        dom.append(dom1);

        // check circle
        Assert.assertEquals(2.0 * FastMath.PI, dom.getInitialTime(), 1.0e-12);
        Assert.assertEquals(0, dom.getFinalTime(), 1.0e-12);
        for (double t = 0; t < 2.0 * FastMath.PI; t += 0.1) {
            final double[] y = dom.getInterpolatedState(t).getPrimaryState();
            Assert.assertEquals(FastMath.cos(t), y[0], 1.0e-7);
            Assert.assertEquals(FastMath.sin(t), y[1], 1.0e-7);
        }

    }

    @Test
    public void testErrorConditions() throws MathIllegalArgumentException, MathIllegalStateException {

        DenseOutputModel cm = new DenseOutputModel();
        cm.handleStep(buildInterpolator(0, new double[] { 0.0, 1.0, -2.0 }, 1));

        // dimension mismatch
        Assert.assertTrue(checkAppendError(cm, 1.0, new double[] { 0.0, 1.0 }, 2.0));

        // hole between time ranges
        Assert.assertTrue(checkAppendError(cm, 10.0, new double[] { 0.0, 1.0, -2.0 }, 20.0));

        // propagation direction mismatch
        Assert.assertTrue(checkAppendError(cm, 1.0, new double[] { 0.0, 1.0, -2.0 }, 0.0));

        // no errors
        Assert.assertFalse(checkAppendError(cm, 1.0, new double[] { 0.0, 1.0, -2.0 }, 2.0));

    }

    @Test
    public void testSerialization() {
        try {
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

            int expectedSize = 131976;
            Assert.assertTrue("size = " + bos.size(), bos.size () >  9 * expectedSize / 10);
            Assert.assertTrue("size = " + bos.size(), bos.size () < 11 * expectedSize / 10);

            ByteArrayInputStream  bis = new ByteArrayInputStream(bos.toByteArray());
            ObjectInputStream     ois = new ObjectInputStream(bis);
            DenseOutputModel cm  = (DenseOutputModel) ois.readObject();

            Random random = new Random(347588535632l);
            double maxError = 0.0;
            for (int i = 0; i < 1000; ++i) {
                double r = random.nextDouble();
                double time = r * pb.getInitialTime() + (1.0 - r) * pb.getFinalTime();
                double[] interpolatedY = cm.getInterpolatedState(time).getPrimaryState();
                double[] theoreticalY  = pb.computeTheoreticalState(time);
                double dx = interpolatedY[0] - theoreticalY[0];
                double dy = interpolatedY[1] - theoreticalY[1];
                double error = dx * dx + dy * dy;
                if (error > maxError) {
                    maxError = error;
                }
            }
            Assert.assertEquals(0.0, maxError, 5.5e-7);
        } catch (ClassNotFoundException | IOException e) {
            Assert.fail(e.getLocalizedMessage());
        }

    }
    private boolean checkAppendError(DenseOutputModel cm,
                                     double t0, double[] y0, double t1)
                                                     throws MathIllegalArgumentException, MathIllegalStateException {
        try {
            DenseOutputModel otherCm = new DenseOutputModel();
            otherCm.handleStep(buildInterpolator(t0, y0, t1));
            cm.append(otherCm);
        } catch(MathIllegalArgumentException iae) {
            return true; // there was an allowable error
        }
        return false; // no allowable error
    }

    private ODEStateInterpolator buildInterpolator(double t0, double[] y0, double t1) {
        return new DummyStepInterpolator(t1 >= t0,
                                         new ODEStateAndDerivative(t0, y0,  new double[y0.length]),
                                         new ODEStateAndDerivative(t1, y0,  new double[y0.length]),
                                         new ODEStateAndDerivative(t0, y0,  new double[y0.length]),
                                         new ODEStateAndDerivative(t1, y0,  new double[y0.length]),
                                         new EquationsMapper(null, y0.length));
    }

    public void checkValue(double value, double reference) {
        Assert.assertTrue(FastMath.abs(value - reference) < 1.0e-10);
    }

    @Before
    public void setUp() {
        pb = new TestProblem3(0.9);
        double minStep = 0;
        double maxStep = pb.getFinalTime() - pb.getInitialTime();
        integ = new DormandPrince54Integrator(minStep, maxStep, 1.0e-8, 1.0e-8);
    }

    @After
    public void tearDown() {
        pb    = null;
        integ = null;
    }

}
