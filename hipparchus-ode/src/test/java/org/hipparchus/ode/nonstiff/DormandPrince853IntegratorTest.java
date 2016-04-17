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
import org.junit.Assert;
import org.junit.Test;


public class DormandPrince853IntegratorTest extends EmbeddedRungeKuttaIntegratorAbstractTest {

    protected EmbeddedRungeKuttaIntegrator
    createIntegrator(final double minStep, final double maxStep,
                     final double scalAbsoluteTolerance, final double scalRelativeTolerance) {
        return new DormandPrince853Integrator(minStep, maxStep, scalAbsoluteTolerance, scalRelativeTolerance);
    }

    protected EmbeddedRungeKuttaIntegrator
    createIntegrator(final double minStep, final double maxStep,
                     final double[] vecAbsoluteTolerance, final double[] vecRelativeTolerance) {
        return new DormandPrince853Integrator(minStep, maxStep, vecAbsoluteTolerance, vecRelativeTolerance);
    }

    @Override
    public void testBackward() {
        doTestBackward(8.1e-8, 1.1e-7, 1.0e-12, "Dormand-Prince 8 (5, 3)");
    }

    @Override
    public void testKepler() {
        doTestKepler(4.4e-11);
    }

    @Override
    public void testForwardBackwardExceptions() {
        doTestForwardBackwardExceptions();
    }

    @Override
    public void testIncreasingTolerance() {
        // the 1.3 factor is only valid for this test
        // and has been obtained from trial and error
        // there is no general relation between local and global errors
        doTestIncreasingTolerance(1.3, 1.0e-12);
    }

    @Override
    public void testEvents() {
        doTestEvents(2.1e-7, "Dormand-Prince 8 (5, 3)");
    }

    @Test
    public void testMissedEndEvent() {
        doTestMissedEndEvent(1.0e-15, 1.0e-15);
    }

    @Test
    public void testVariableSteps() {
        doTestVariableSteps(0.00763, 0.836);
    }

    @Test
    public void testUnstableDerivative() {
     doTestUnstableDerivative(1.0e-12);
    }

    @Override
    public void testPartialDerivatives() {
        doTestPartialDerivatives(2.6e-12, 2.0e-11);
    }

    @Deprecated
    @Test
    public void testDeprecatedInterfaces()
                    throws MathIllegalArgumentException, MathIllegalStateException {

        final DeprecatedODE pb = new DeprecatedODE();
        double minStep = 0;
        double maxStep = 1.0;
        double scalAbsoluteTolerance = 1.0e-8;
        double scalRelativeTolerance = scalAbsoluteTolerance;

        ODEIntegrator integ = new DormandPrince853Integrator(minStep, maxStep,
                                                             scalAbsoluteTolerance,
                                                             scalRelativeTolerance);
        final DeprecatedStepHandler stepHandler = new DeprecatedStepHandler(integ);
        integ.addStepHandler(stepHandler);
        final DeprecatedEventHandler eventHandler = new DeprecatedEventHandler();
        integ.addEventHandler(eventHandler, 1.0, 1.0e-10, 100);
        double[] y = new double[1];
        double finalT = integ.integrate(pb, 0.0, y, 10.0, y);

        Assert.assertTrue(pb.initCalled);
        Assert.assertTrue(stepHandler.initCalled);
        Assert.assertTrue(stepHandler.lastStepSeen);
        Assert.assertTrue(eventHandler.initCalled);
        Assert.assertTrue(eventHandler.resetCalled);

        Assert.assertEquals(4.0, finalT, 1.0e-10);

    }

    @Deprecated
    private static class DeprecatedODE
    implements org.hipparchus.ode.FirstOrderDifferentialEquations {

        private boolean initCalled = false;

        public void init(double t0, double[] y0, double finalTime) {
            initCalled = true;
        }

        public int getDimension() {
            return 1;
        }

        public void computeDerivatives(double t, double[] y, double[] yDot) {
            yDot[0] = 1.0;
        }

    }

    @Deprecated
    private static class DeprecatedStepHandler
    implements org.hipparchus.ode.sampling.StepHandler {

        private final ODEIntegrator integrator;
        private boolean initCalled   = false;
        private boolean lastStepSeen = false;

        public DeprecatedStepHandler(final ODEIntegrator integrator) {
            this.integrator = integrator;
        }

        public void init(double t0, double[] y0, double t) {
            initCalled = true;
        }

        public void handleStep(org.hipparchus.ode.sampling.StepInterpolator interpolator, boolean isLast) {
            Assert.assertEquals(interpolator.getPreviousTime(),
                                integrator.getCurrentStepStart(),
                                1.0e-10);
            if (isLast) {
                lastStepSeen = true;
            }
        }

    }

    @Deprecated
    private static class DeprecatedEventHandler
    implements org.hipparchus.ode.events.EventHandler {

        private boolean initCalled = false;
        private boolean resetCalled = false;

        public void init(double t0, double[] y0, double t) {
            initCalled = true;
        }

        public double g(double t, double[] y) {
            return (t - 2.0) * (t - 4.0);
        }

        public Action eventOccurred(double t, double[] y, boolean increasing) {
            return t < 3 ? Action.RESET_STATE : Action.STOP;
        }

        public void resetState(double t, double[] y) {
            resetCalled = true;
        }

    }

}
