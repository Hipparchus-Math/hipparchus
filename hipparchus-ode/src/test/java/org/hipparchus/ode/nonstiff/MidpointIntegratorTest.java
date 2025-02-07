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

import org.junit.jupiter.api.Test;

class MidpointIntegratorTest extends RungeKuttaIntegratorAbstractTest {

    protected FixedStepRungeKuttaIntegrator createIntegrator(double step) {
        return new MidpointIntegrator(step);
    }

    @Override
    @Test
    public void testMissedEndEvent() {
        doTestMissedEndEvent(1.0e-15, 6.0e-5);
    }

    @Override
    @Test
    public void testSanityChecks() {
        doTestSanityChecks();
    }

    @Override
    @Test
    public void testDecreasingSteps() {
        doTestDecreasingSteps(1.0, 1.0, 1.0e-10);
    }

    @Override
    @Test
    public void testSmallStep() {
        doTestSmallStep(2.0e-7, 1.0e-6, 1.0e-12, "midpoint");
    }

    @Override
    @Test
    public void testBigStep() {
        doTestBigStep(0.01, 0.05, 1.0e-12, "midpoint");

    }

    @Override
    @Test
    public void testBackward() {
        doTestBackward(6.0e-4, 6.0e-4, 1.0e-12, "midpoint");
    }

    @Override
    @Test
    public void testKepler() {
        doTestKepler(1.19, 0.01);
    }

    @Override
    @Test
    public void testStepSize() {
        doTestStepSize(1.0e-12);
    }

    @Override
    @Test
    public void testSingleStep() {
        doTestSingleStep(0.21);
    }

    @Override
    @Test
    public void testTooLargeFirstStep() {
        doTestTooLargeFirstStep();
    }

    @Override
    @Test
    public void testUnstableDerivative() {
        doTestUnstableDerivative(1.0e-12);
    }

    @Override
    @Test
    public void testDerivativesConsistency() {
        doTestDerivativesConsistency(1.0e-10);
    }

    @Override
    @Test
    public void testSerialization() {
        doTestSerialization(737825, 1.19);
    }

    @Test
    public void testSecondaryEquations() {
        doTestSecondaryEquations(1.6e-6, 5.6e-13);
    }

}
