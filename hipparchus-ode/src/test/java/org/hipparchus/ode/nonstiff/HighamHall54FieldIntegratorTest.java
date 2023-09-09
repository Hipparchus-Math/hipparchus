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


import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.util.Binary64Field;
import org.junit.Test;

public class HighamHall54FieldIntegratorTest extends EmbeddedRungeKuttaFieldIntegratorAbstractTest {

    protected <T extends CalculusFieldElement<T>> EmbeddedRungeKuttaFieldIntegrator<T>
    createIntegrator(Field<T> field, final double minStep, final double maxStep,
                     final double scalAbsoluteTolerance, final double scalRelativeTolerance) {
        return new HighamHall54FieldIntegrator<T>(field, minStep, maxStep, scalAbsoluteTolerance, scalRelativeTolerance);
    }

    protected <T extends CalculusFieldElement<T>> EmbeddedRungeKuttaFieldIntegrator<T>
    createIntegrator(Field<T> field, final double minStep, final double maxStep,
                     final double[] vecAbsoluteTolerance, final double[] vecRelativeTolerance) {
        return new HighamHall54FieldIntegrator<T>(field, minStep, maxStep, vecAbsoluteTolerance, vecRelativeTolerance);
    }

    @Override
    public void testNonFieldIntegratorConsistency() {
        doTestNonFieldIntegratorConsistency(Binary64Field.getInstance());
    }

    @Override
    public void testSanityChecks() {
        doTestSanityChecks(Binary64Field.getInstance());
    }

    @Override
    public void testBackward() {
        doTestBackward(Binary64Field.getInstance(), 5.0e-7, 5.0e-7, 1.0e-12, "Higham-Hall 5(4)");
    }

    @Override
    public void testKepler() {
        doTestKepler(Binary64Field.getInstance(), 1.5e-4);
    }

    @Override
    public void testTorqueFreeMotionOmegaOnly() {
        doTestTorqueFreeMotionOmegaOnly(Binary64Field.getInstance(), 3.0e-15);
    }

    @Override
    public void testTorqueFreeMotion() {
        doTestTorqueFreeMotion(Binary64Field.getInstance(), 1.6e-15, 5.6e-16);
    }

    @Override
    public void testTorqueFreeMotionIssue230() {
        doTestTorqueFreeMotionIssue230(Binary64Field.getInstance(), 5.4e-15, 1.5e-15);
    }

    @Override
    public void testForwardBackwardExceptions() {
        doTestForwardBackwardExceptions(Binary64Field.getInstance());
    }

    @Override
    public void testMinStep() {
        doTestMinStep(Binary64Field.getInstance());
    }

    @Override
    public void testIncreasingTolerance() {
        // the 1.3 factor is only valid for this test
        // and has been obtained from trial and error
        // there is no general relation between local and global errors
        doTestIncreasingTolerance(Binary64Field.getInstance(), 1.3, 1.0e-12);
    }

    @Override
    public void testEvents() {
        doTestEvents(Binary64Field.getInstance(), 1.0e-7, "Higham-Hall 5(4)");
    }

    @Override
    public void testStepEnd() {
        doTestStepEnd(Binary64Field.getInstance(), 119, "Higham-Hall 5(4)");
    }

    @Override
    public void testStopAfterStep() {
        doTestStopAfterStep(Binary64Field.getInstance(), 12, 1.110791);
    }

    @Override
    public void testResetAfterStep() {
        doTestResetAfterStep(Binary64Field.getInstance(), 12, 14);
    }

    @Override
    public void testEventsErrors() {
        doTestEventsErrors(Binary64Field.getInstance());
    }

    @Override
    public void testEventsNoConvergence() {
        doTestEventsNoConvergence(Binary64Field.getInstance());
    }

    @Override
    public void testPartialDerivatives() {
        doTestPartialDerivatives(1.2e-11, new double[] { 8.0e-11, 1.8e-11, 2.6e-12, 2.6e-12, 1.8e-11 });
    }

    @Test
    public void testSecondaryEquations() {
        doTestSecondaryEquations(Binary64Field.getInstance(), 9.7e-12, 1.3e-14);
    }

}
