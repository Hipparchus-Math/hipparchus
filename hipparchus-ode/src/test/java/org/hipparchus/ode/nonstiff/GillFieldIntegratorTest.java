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

public class GillFieldIntegratorTest extends RungeKuttaFieldIntegratorAbstractTest {

    protected <T extends CalculusFieldElement<T>> RungeKuttaFieldIntegrator<T>
    createIntegrator(Field<T> field, T step) {
        return new GillFieldIntegrator<T>(field, step);
    }

    @Override
    public void testNonFieldIntegratorConsistency() {
        doTestNonFieldIntegratorConsistency(Binary64Field.getInstance());
    }

    @Override
    public void testMissedEndEvent() {
        doTestMissedEndEvent(Binary64Field.getInstance(), 1.0e-15, 6.0e-5);
    }

    @Override
    public void testSanityChecks() {
        doTestSanityChecks(Binary64Field.getInstance());
    }

    @Override
    public void testDecreasingSteps() {
        doTestDecreasingSteps(Binary64Field.getInstance(), 1.0, 1.0, 1.0e-10);
    }

    @Override
    public void testSmallStep() {
        doTestSmallStep(Binary64Field.getInstance(), 2.0e-13, 4.0e-12, 1.0e-12, "Gill");
    }

    @Override
    public void testBigStep() {
        doTestBigStep(Binary64Field.getInstance(), 0.0004, 0.005, 1.0e-12, "Gill");

    }

    @Override
    public void testBackward() {
        doTestBackward(Binary64Field.getInstance(), 5.0e-10, 7.0e-10, 1.0e-12, "Gill");
    }

    @Override
    public void testKepler() {
        doTestKepler(Binary64Field.getInstance(), 1.72e-3, 1.0e-5);
    }

    @Override
    public void testStepSize() {
        doTestStepSize(Binary64Field.getInstance(), 1.0e-12);
    }

    @Override
    public void testSingleStep() {
        doTestSingleStep(Binary64Field.getInstance(), 0.21);
    }

    @Override
    public void testTooLargeFirstStep() {
        doTestTooLargeFirstStep(Binary64Field.getInstance());
    }

    @Override
    public void testUnstableDerivative() {
        doTestUnstableDerivative(Binary64Field.getInstance(), 1.0e-12);
    }

    @Override
    public void testDerivativesConsistency() {
        doTestDerivativesConsistency(Binary64Field.getInstance(), 1.0e-10);
    }

    @Override
    public void testPartialDerivatives() {
        doTestPartialDerivatives(3.2e-10, new double[] { 2.1e-9, 5.9e-10, 7.0e-11, 7.0e-11, 5.9e-10 });
    }

    @Test
    public void testSecondaryEquations() {
        doTestSecondaryEquations(Binary64Field.getInstance(), 1.1e-12, 5.6e-13);
    }

}
