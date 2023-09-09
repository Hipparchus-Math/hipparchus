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

public class EulerIntegratorTest extends RungeKuttaIntegratorAbstractTest {

    protected RungeKuttaIntegrator createIntegrator(double step) {
        return new EulerIntegrator(step);
    }

    @Override
    public void testMissedEndEvent() {
        doTestMissedEndEvent(1.0e-15, 6.0e-5);
    }

    @Override
    public void testSanityChecks() {
        doTestSanityChecks();
    }

    @Override
    public void testDecreasingSteps() {
        doTestDecreasingSteps(1.0, 1.5, 1.0e-10);
    }

    @Override
    public void testSmallStep() {
        doTestSmallStep(2.0e-4, 1.0e-3, 1.0e-12, "Euler");
    }

    @Override
    public void testBigStep() {
        doTestBigStep(0.01, 0.2, 1.0e-12, "Euler");

    }

    @Override
    public void testBackward() {
        doTestBackward(0.45, 0.45, 1.0e-12, "Euler");
    }

    @Override
    public void testKepler() {
        // Euler integrator is clearly not able to solve this problem
        doTestKepler(881.176, 0.001);
    }

    @Override
    public void testStepSize() {
        doTestStepSize(1.0e-12);
    }

    @Override
    public void testSingleStep() {
        doTestSingleStep(0.21);
    }

    @Override
    public void testTooLargeFirstStep() {
        doTestTooLargeFirstStep();
    }

    @Override
    public void testUnstableDerivative() {
        doTestUnstableDerivative(1.0e-12);
    }

    @Override
    public void testDerivativesConsistency() {
        doTestDerivativesConsistency(1.0e-10);
    }

    @Override
    public void testSerialization() {
        doTestSerialization(597794, 881);
    }

    @Test
    public void testSecondaryEquations() {
        doTestSecondaryEquations(4.8e-3, 5.6e-13);
    }

}
