/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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


import org.junit.Test;

public class HighamHall54IntegratorTest extends EmbeddedRungeKuttaIntegratorAbstractTest {

    protected EmbeddedRungeKuttaIntegrator
    createIntegrator(final double minStep, final double maxStep,
                     final double scalAbsoluteTolerance, final double scalRelativeTolerance) {
        return new HighamHall54Integrator(minStep, maxStep, scalAbsoluteTolerance, scalRelativeTolerance);
    }

    protected EmbeddedRungeKuttaIntegrator
    createIntegrator(final double minStep, final double maxStep,
                     final double[] vecAbsoluteTolerance, final double[] vecRelativeTolerance) {
        return new HighamHall54Integrator(minStep, maxStep, vecAbsoluteTolerance, vecRelativeTolerance);
    }

    @Override
    public void testBackward() {
        doTestBackward(5.0e-7, 5.0e-7, 1.0e-12, "Higham-Hall 5(4)");
    }

    @Override
    public void testKepler() {
        doTestKepler(1.5e-4);
    }

    @Override
    public void testTorqueFreeMotion() {
        doTestTorqueFreeMotion(1.0e-15, 1.0e-15);
    }

    @Override
    public void testTorqueFreeMotionDebug() {
        doTestTorqueFreeMotionDebug(1.0e-15, 1.0e-15);
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
        doTestEvents(1.0e-7, "Higham-Hall 5(4)");
    }

    @Test
    public void testMissedEndEvent() {
        doTestMissedEndEvent(1.0e-15, 1.0e-15);
    }

    @Test
    public void testVariableSteps() {
        doTestVariableSteps(0.00169, 0.267);
    }

    @Test
    public void testUnstableDerivative() {
     doTestUnstableDerivative(1.0e-12);
    }

    @Override
    public void testPartialDerivatives() {
        doTestPartialDerivatives(1.2e-11, 8.0e-11);
    }

    @Test
    public void testSecondaryEquations() {
        doTestSecondaryEquations(9.7e-12, 1.3e-14);
    }

	@Override
	public void testTorqueFreeMotionOmegaOnly() {
		// TODO Auto-generated method stub
		
	}

}
