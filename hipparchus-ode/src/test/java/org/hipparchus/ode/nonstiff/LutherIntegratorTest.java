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

public class LutherIntegratorTest extends RungeKuttaIntegratorAbstractTest {

    protected RungeKuttaIntegrator createIntegrator(double step) {
        return new LutherIntegrator(step);
    }

    @Override
    public void testMissedEndEvent() {
        doTestMissedEndEvent(1.0e-15, 1.0e-15);
    }

    @Override
    public void testSanityChecks() {
        doTestSanityChecks();
    }

    @Override
    public void testDecreasingSteps() {
        doTestDecreasingSteps(3.6, 1.0, 1.0e-10);
    }

    @Override
    public void testSmallStep() {
        doTestSmallStep(8.7e-17, 3.6e-15, 1.0e-12, "Luther");
    }

    @Override
    public void testBigStep() {
        doTestBigStep(2.7e-5, 1.7e-3, 1.0e-12, "Luther");
    }

    @Override
    public void testBackward() {
        doTestBackward(2.4e-13, 4.3e-13, 1.0e-12, "Luther");
    }

    @Override
    public void testKepler() {
        doTestKepler(2.18e-7, 4.0e-10);
    }

    @Override
    public void testStepSize() {
        doTestStepSize(1.0e-22);
    }

    @Override
    public void testSingleStep() {
        doTestSingleStep(6.0e-12);
    }

    @Override
    public void testTooLargeFirstStep() {
        doTestTooLargeFirstStep();
    }

    @Override
    public void testUnstableDerivative() {
        doTestUnstableDerivative(4.0e-15);
    }

    @Override
    public void testDerivativesConsistency() {
        doTestDerivativesConsistency(1.0e-20);
    }

    @Override
    public void testSerialization() {
        doTestSerialization(1437963, 2.18e-7);
    }

    @Test
    public void testSecondaryEquations() {
        doTestSecondaryEquations(1.0e-12, 5.6e-13);
    }

}
