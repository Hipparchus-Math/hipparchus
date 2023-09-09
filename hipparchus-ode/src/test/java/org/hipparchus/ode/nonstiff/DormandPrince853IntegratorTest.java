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
    public void testTorqueFreeMotionOmegaOnly() {
        doTestTorqueFreeMotionOmegaOnly(4.0e-16);
    }

    @Override
    public void testTorqueFreeMotion() {
        doTestTorqueFreeMotion(1.3e-12, 9.0e-12);
    }

    @Override
    public void testTorqueFreeMotionIssue230() {
        doTestTorqueFreeMotionIssue230(2.9e-14, 7.8e-14);
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

    @Override
    public void testStepEnd() {
        doTestStepEnd(20, "Dormand-Prince 8 (5, 3)");
    }

    @Override
    public void testStopAfterStep() {
        doTestStopAfterStep(12, 6.842171);
    }

    @Override
    public void testResetAfterStep() {
        doTestResetAfterStep(12, 13);
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

    @Test
    public void testSecondaryEquations() {
        doTestSecondaryEquations(3.3e-12, 8.9e-15);
    }

}