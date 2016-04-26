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

package org.hipparchus.migration.ode;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.ode.nonstiff.DormandPrince853IntegratorTest;
import org.junit.Assert;
import org.junit.Test;


@Deprecated
public class IntegratorDeprecatedMethodTest extends DormandPrince853IntegratorTest {

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

    private static class DeprecatedODE
    implements org.hipparchus.migration.ode.FirstOrderDifferentialEquations {

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

    private static class DeprecatedStepHandler
    implements org.hipparchus.migration.ode.StepHandler {

        private final ODEIntegrator integrator;
        private boolean initCalled   = false;
        private boolean lastStepSeen = false;

        public DeprecatedStepHandler(final ODEIntegrator integrator) {
            this.integrator = integrator;
        }

        public void init(double t0, double[] y0, double t) {
            initCalled = true;
        }

        public void handleStep(MigrationStepInterpolator interpolator, boolean isLast) {
            Assert.assertEquals(interpolator.getPreviousTime(),
                                integrator.getCurrentStepStart(),
                                1.0e-10);
            if (isLast) {
                lastStepSeen = true;
            }
        }

    }

    private static class DeprecatedEventHandler
    implements org.hipparchus.migration.ode.events.EventHandler {

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
