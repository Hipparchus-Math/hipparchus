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


import org.hipparchus.Field;
import org.hipparchus.RealFieldElement;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.Decimal64Field;
import org.junit.Test;

public class AdamsMoultonFieldIntegratorTest extends AdamsFieldIntegratorAbstractTest {

    protected <T extends RealFieldElement<T>> AdamsFieldIntegrator<T>
    createIntegrator(Field<T> field, final int nSteps, final double minStep, final double maxStep,
                     final double scalAbsoluteTolerance, final double scalRelativeTolerance) {
        return new AdamsMoultonFieldIntegrator<T>(field, nSteps, minStep, maxStep,
                        scalAbsoluteTolerance, scalRelativeTolerance);
    }

    protected <T extends RealFieldElement<T>> AdamsFieldIntegrator<T>
    createIntegrator(Field<T> field, final int nSteps, final double minStep, final double maxStep,
                     final double[] vecAbsoluteTolerance, final double[] vecRelativeTolerance) {
        return new AdamsMoultonFieldIntegrator<T>(field, nSteps, minStep, maxStep,
                        vecAbsoluteTolerance, vecRelativeTolerance);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testMinStep() {
        doDimensionCheck(Decimal64Field.getInstance());
    }

    @Test
    public void testIncreasingTolerance() {
        // the 0.45 and 8.69 factors are only valid for this test
        // and has been obtained from trial and error
        // there are no general relationship between local and global errors
        doTestIncreasingTolerance(Decimal64Field.getInstance(), 0.45, 8.69);
    }

    @Test(expected = MathIllegalStateException.class)
    public void exceedMaxEvaluations() {
        doExceedMaxEvaluations(Decimal64Field.getInstance(), 650);
    }

    @Test
    public void backward() {
        doBackward(Decimal64Field.getInstance(), 3.0e-9, 3.0e-9, 1.0e-16, "Adams-Moulton");
    }

    @Test
    public void polynomial() {
        doPolynomial(Decimal64Field.getInstance(), 5, 2.2e-05, 1.1e-11);
    }

    @Test
    public void testSecondaryEquations() {
        doTestSecondaryEquations(Decimal64Field.getInstance(), 1.9e-11, 7.2e-15);
    }

    @Test(expected=MathIllegalStateException.class)
    public void testStartFailure() {
        doTestStartFailure(Decimal64Field.getInstance());
    }

}
