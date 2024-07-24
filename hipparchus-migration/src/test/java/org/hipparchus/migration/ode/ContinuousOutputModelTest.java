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

package org.hipparchus.migration.ode;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.EquationsMapper;
import org.hipparchus.ode.ExpandableODE;
import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.TestProblem3;
import org.hipparchus.ode.nonstiff.DormandPrince54Integrator;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.ode.sampling.DummyStepInterpolator;
import org.hipparchus.ode.sampling.ODEStateInterpolator;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

@Deprecated
public class ContinuousOutputModelTest {

    TestProblem3 pb;
    ODEIntegrator integ;

    @Test
    void testBoundaries() throws MathIllegalArgumentException, MathIllegalStateException {
        integ.addStepHandler(new ContinuousOutputModel());
        integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());
        ContinuousOutputModel cm = (ContinuousOutputModel) integ.getStepHandlers().iterator().next();
        cm.setInterpolatedTime(2.0 * pb.getInitialTime() - pb.getFinalTime());
        cm.setInterpolatedTime(2.0 * pb.getFinalTime() - pb.getInitialTime());
        cm.setInterpolatedTime(0.5 * (pb.getFinalTime() + pb.getInitialTime()));
    }

    @Test
    void testRandomAccess() throws MathIllegalArgumentException, MathIllegalStateException {

        ContinuousOutputModel cm = new ContinuousOutputModel();
        integ.addStepHandler(cm);
        integ.integrate(pb, pb.getInitialState(), pb.getFinalTime());

        Random random = new Random(347588535632l);
        double maxError    = 0.0;
        double maxErrorDot = 0.0;
        for (int i = 0; i < 1000; ++i) {
            double r = random.nextDouble();
            double time = r * pb.getInitialTime() + (1.0 - r) * pb.getFinalTime();
            cm.setInterpolatedTime(time);
            double[] interpolatedY    = cm.getInterpolatedState();
            double[] interpolatedYDot = cm.getInterpolatedDerivatives();
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

        assertEquals(0.0, maxError,    1.0e-9);
        assertEquals(0.0, maxErrorDot, 4.0e-7);

    }

    @Test
    void testModelsMerging() throws MathIllegalArgumentException, MathIllegalStateException {

        // theoretical solution: y[0] = cos(t), y[1] = sin(t)
        FirstOrderDifferentialEquations problem =
                        new FirstOrderDifferentialEquations() {
            @Override
            public void computeDerivatives(double t, double[] y, double[] dot) {
                dot[0] = -y[1];
                dot[1] =  y[0];
            }
            @Override
            public int getDimension() {
                return 2;
            }
        };

        // integrate backward from &pi; to 0;
        ContinuousOutputModel cm1 = new ContinuousOutputModel();
        ODEIntegrator integ1 = new DormandPrince853Integrator(0, 1.0, 1.0e-8, 1.0e-8);
        integ1.addStepHandler(cm1);
        integ1.integrate(new ExpandableODE(problem),
                         new ODEState(FastMath.PI, new double[] { -1.0, 0.0 }),
                         0);

        // integrate backward from 2&pi; to &pi;
        ContinuousOutputModel cm2 = new ContinuousOutputModel();
        ODEIntegrator integ2 = new DormandPrince853Integrator(0, 0.1, 1.0e-12, 1.0e-12);
        integ2.addStepHandler(cm2);
        integ2.integrate(new ExpandableODE(problem),
                         new ODEState(2.0 * FastMath.PI, new double[] { 1.0, 0.0 }),
                         FastMath.PI);

        // merge the two half circles
        ContinuousOutputModel cm = new ContinuousOutputModel();
        cm.append(cm2);
        cm.append(new ContinuousOutputModel());
        cm.append(cm1);

        // check circle
        assertEquals(2.0 * FastMath.PI, cm.getInitialTime(), 1.0e-12);
        assertEquals(0, cm.getFinalTime(), 1.0e-12);
        assertEquals(cm.getFinalTime(), cm.getInterpolatedTime(), 1.0e-12);
        for (double t = 0; t < 2.0 * FastMath.PI; t += 0.1) {
            cm.setInterpolatedTime(t);
            double[] y = cm.getInterpolatedState();
            assertEquals(FastMath.cos(t), y[0], 1.0e-7);
            assertEquals(FastMath.sin(t), y[1], 1.0e-7);
        }

    }

    @Test
    void testErrorConditions() throws MathIllegalArgumentException, MathIllegalStateException {

        ContinuousOutputModel cm = new ContinuousOutputModel();
        cm.handleStep(buildInterpolator(0, new double[] { 0.0, 1.0, -2.0 }, 1));

        // dimension mismatch
        assertTrue(checkAppendError(cm, 1.0, new double[] { 0.0, 1.0 }, 2.0));

        // hole between time ranges
        assertTrue(checkAppendError(cm, 10.0, new double[] { 0.0, 1.0, -2.0 }, 20.0));

        // propagation direction mismatch
        assertTrue(checkAppendError(cm, 1.0, new double[] { 0.0, 1.0, -2.0 }, 0.0));

        // no errors
        assertFalse(checkAppendError(cm, 1.0, new double[] { 0.0, 1.0, -2.0 }, 2.0));

    }

    private boolean checkAppendError(ContinuousOutputModel cm,
                                     double t0, double[] y0, double t1)
                                                     throws MathIllegalArgumentException, MathIllegalStateException {
        try {
            ContinuousOutputModel otherCm = new ContinuousOutputModel();
            otherCm.handleStep(buildInterpolator(t0, y0, t1));
            cm.append(otherCm);
        } catch(MathIllegalArgumentException iae) {
            return true; // there was an allowable error
        }
        return false; // no allowable error
    }

    private ODEStateInterpolator buildInterpolator(double t0, double[] y0, double t1) {
        EquationsMapper mapper = null;
        try {
            Constructor<EquationsMapper> ctr;
            ctr = EquationsMapper.class.getDeclaredConstructor(EquationsMapper.class, Integer.TYPE);
            ctr.setAccessible(true);
            mapper = ctr.newInstance(null, y0.length);
        } catch (NoSuchMethodException | SecurityException | InstantiationException |
                 IllegalAccessException | IllegalArgumentException | InvocationTargetException e) {
            fail(e.getLocalizedMessage());
        }
        return new DummyStepInterpolator(t1 >= t0,
                                         new ODEStateAndDerivative(t0, y0,  new double[y0.length]),
                                         new ODEStateAndDerivative(t1, y0,  new double[y0.length]),
                                         new ODEStateAndDerivative(t0, y0,  new double[y0.length]),
                                         new ODEStateAndDerivative(t1, y0,  new double[y0.length]),
                                         mapper);
    }

    public void checkValue(double value, double reference) {
        assertTrue(FastMath.abs(value - reference) < 1.0e-10);
    }

    @BeforeEach
    void setUp() {
        pb = new TestProblem3(0.9);
        double minStep = 0;
        double maxStep = pb.getFinalTime() - pb.getInitialTime();
        integ = new DormandPrince54Integrator(minStep, maxStep, 1.0e-8, 1.0e-8);
    }

    @AfterEach
    void tearDown() {
        pb    = null;
        integ = null;
    }

}
