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

public class ClassicalRungeKuttaIntegratorTest extends RungeKuttaIntegratorAbstractTest {

    protected RungeKuttaIntegrator createIntegrator(double step) {
        return new ClassicalRungeKuttaIntegrator(step);
    }

    @Override
    public void testMissedEndEvent() {
        doTestMissedEndEvent(5.0e-6, 1.0e-9);
    }

    @Override
    public void testSanityChecks() {
        doTestSanityChecks();
    }

    @Override
    public void testDecreasingSteps() {
        doTestDecreasingSteps(1.0, 1.0, 1.0e-10);
    }

    @Override
    public void testSmallStep() {
        doTestSmallStep(2.0e-13, 4.0e-12, 1.0e-12, "classical Runge-Kutta");
    }

    @Override
    public void testBigStep() {
        doTestBigStep(0.0004, 0.005, 1.0e-12, "classical Runge-Kutta");

    }

    @Override
    public void testBackward() {
        doTestBackward(5.0e-10, 7.0e-10, 1.0e-12, "classical Runge-Kutta");
    }

    @Override
    public void testKepler() {
        doTestKepler(5.82e-3, 1.0e-5);
    }

    @Override
    public void testStepSize() {
        doTestStepSize(1.0e-12);
    }

    @Override
    public void testSingleStep() {
        doTestSingleStep(9.3e-9);
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
        doTestSerialization(1017892, 5.5e-3);
    }

    @Test
    public void testSecondaryEquations() {
        doTestSecondaryEquations(1.1e-12, 5.6e-13);
    }

}
