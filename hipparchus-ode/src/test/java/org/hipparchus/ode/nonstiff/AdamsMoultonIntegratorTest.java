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


import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.junit.Test;

public class AdamsMoultonIntegratorTest extends AdamsIntegratorAbstractTest {

    protected AdamsIntegrator
    createIntegrator(final int nSteps, final double minStep, final double maxStep,
                     final double scalAbsoluteTolerance, final double scalRelativeTolerance) {
        return new AdamsMoultonIntegrator(nSteps, minStep, maxStep,
                                          scalAbsoluteTolerance, scalRelativeTolerance);
    }

    protected AdamsIntegrator
    createIntegrator(final int nSteps, final double minStep, final double maxStep,
                     final double[] vecAbsoluteTolerance, final double[] vecRelativeTolerance) {
        return new AdamsMoultonIntegrator(nSteps, minStep, maxStep,
                                          vecAbsoluteTolerance, vecRelativeTolerance);
    }

    @Test
    public void testNbPoints() {
        doNbPointsTest();
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testMinStep() {
        doDimensionCheck();
    }

    @Test
    public void testIncreasingTolerance() {
        // the 0.45 and 8.69 factors are only valid for this test
        // and has been obtained from trial and error
        // there are no general relationship between local and global errors
        doTestIncreasingTolerance(0.45, 8.69);
    }

    @Test(expected = MathIllegalStateException.class)
    public void exceedMaxEvaluations() {
        doExceedMaxEvaluations(650);
    }

    @Test
    public void backward() {
        doBackward(3.0e-9, 3.0e-9, 1.0e-16, "Adams-Moulton");
    }

    @Test
    public void polynomial() {
        doPolynomial(5, 2.2e-05, 2.0e-11);
    }

    @Test
    public void testSecondaryEquations() {
        doTestSecondaryEquations(1.9e-11, 7.2e-15);
    }

    @Test(expected=MathIllegalStateException.class)
    public void testStartFailure() {
        doTestStartFailure();
    }

}
