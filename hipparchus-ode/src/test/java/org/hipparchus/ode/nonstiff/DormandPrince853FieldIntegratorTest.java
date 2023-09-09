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

public class DormandPrince853FieldIntegratorTest extends EmbeddedRungeKuttaFieldIntegratorAbstractTest {

    protected <T extends CalculusFieldElement<T>> EmbeddedRungeKuttaFieldIntegrator<T>
    createIntegrator(Field<T> field, final double minStep, final double maxStep,
                     final double scalAbsoluteTolerance, final double scalRelativeTolerance) {
        return new DormandPrince853FieldIntegrator<T>(field, minStep, maxStep, scalAbsoluteTolerance, scalRelativeTolerance);
    }

    protected <T extends CalculusFieldElement<T>> EmbeddedRungeKuttaFieldIntegrator<T>
    createIntegrator(Field<T> field, final double minStep, final double maxStep,
                     final double[] vecAbsoluteTolerance, final double[] vecRelativeTolerance) {
        return new DormandPrince853FieldIntegrator<T>(field, minStep, maxStep, vecAbsoluteTolerance, vecRelativeTolerance);
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
        doTestBackward(Binary64Field.getInstance(), 8.1e-8, 1.1e-7, 1.0e-12, "Dormand-Prince 8 (5, 3)");
    }

    @Override
    public void testKepler() {
        doTestKepler(Binary64Field.getInstance(), 4.4e-11);
    }

    @Override
    public void testTorqueFreeMotionOmegaOnly() {
        doTestTorqueFreeMotionOmegaOnly(Binary64Field.getInstance(), 4.0e-16);
    }

    @Override
    public void testTorqueFreeMotion() {
        doTestTorqueFreeMotion(Binary64Field.getInstance(), 1.3e-12, 9.0e-12);
    }

    @Override
    public void testTorqueFreeMotionIssue230() {
        doTestTorqueFreeMotionIssue230(Binary64Field.getInstance(), 2.9e-14, 7.8e-14);
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
        doTestEvents(Binary64Field.getInstance(), 2.1e-7, "Dormand-Prince 8 (5, 3)");
    }

    @Override
    public void testStepEnd() {
        doTestStepEnd(Binary64Field.getInstance(), 20, "Dormand-Prince 8 (5, 3)");
    }

    @Override
    public void testStopAfterStep() {
        doTestStopAfterStep(Binary64Field.getInstance(), 12, 6.842171);
    }

    @Override
    public void testResetAfterStep() {
        doTestResetAfterStep(Binary64Field.getInstance(), 12, 13);
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
        doTestPartialDerivatives(2.6e-12, new double[] { 1.9e-11, 3.6e-12, 5.3e-13, 5.3e-13, 3.6e-12 });
    }

    @Test
    public void testSecondaryEquations() {
        doTestSecondaryEquations(Binary64Field.getInstance(), 3.3e-12, 8.9e-15);
    }

}
