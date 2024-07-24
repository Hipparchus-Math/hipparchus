/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */

package org.hipparchus.ode;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.MathArrays;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertThrows;

public class FieldExpandableODETest {

    @Test
    public void testOnlyMainEquation() {
        doTestOnlyMainEquation(Binary64Field.getInstance());
    }

    private <T extends CalculusFieldElement<T>> void doTestOnlyMainEquation(final Field<T> field) {
        FieldOrdinaryDifferentialEquation<T> main = new Linear<T>(field, 3, 0);
        FieldExpandableODE<T> equation = new FieldExpandableODE<T>(main);
        Assertions.assertEquals(main.getDimension(), equation.getMapper().getTotalDimension());
        Assertions.assertEquals(1, equation.getMapper().getNumberOfEquations());
        T t0 = field.getZero().add(10);
        T t  = field.getZero().add(100);
        T[] complete    = MathArrays.buildArray(field, equation.getMapper().getTotalDimension());
        for (int i = 0; i < complete.length; ++i) {
            complete[i] = field.getZero().add(i);
        }
        T[] completeDot = equation.computeDerivatives(t0, complete);
        equation.init(equation.getMapper().mapStateAndDerivative(t0, complete, completeDot), t);
        FieldODEStateAndDerivative<T> state = equation.getMapper().mapStateAndDerivative(t0, complete, completeDot);
        Assertions.assertEquals(0, state.getNumberOfSecondaryStates());
        T[] mainState    = state.getPrimaryState();
        T[] mainStateDot = state.getPrimaryDerivative();
        Assertions.assertEquals(main.getDimension(), mainState.length);
        for (int i = 0; i < main.getDimension(); ++i) {
            Assertions.assertEquals(i, mainState[i].getReal(),   1.0e-15);
            Assertions.assertEquals(i, mainStateDot[i].getReal(), 1.0e-15);
            Assertions.assertEquals(i, completeDot[i].getReal(),  1.0e-15);
        }
    }

    @Test
    public void testPrimaryAndSecondary() {
        doTestPrimaryAndSecondary(Binary64Field.getInstance());
    }

    private <T extends CalculusFieldElement<T>> void doTestPrimaryAndSecondary(final Field<T> field) {

        FieldOrdinaryDifferentialEquation<T> main = new Linear<T>(field, 3, 0);
        FieldExpandableODE<T> equation = new FieldExpandableODE<T>(main);
        FieldSecondaryODE<T> secondary1 = new Linear<T>(field, 3, main.getDimension());
        int i1 = equation.addSecondaryEquations(secondary1);
        FieldSecondaryODE<T> secondary2 = new Linear<T>(field, 5, main.getDimension() + secondary1.getDimension());
        int i2 = equation.addSecondaryEquations(secondary2);
        Assertions.assertEquals(main.getDimension() + secondary1.getDimension() + secondary2.getDimension(),
                            equation.getMapper().getTotalDimension());
        Assertions.assertEquals(3, equation.getMapper().getNumberOfEquations());
        Assertions.assertEquals(1, i1);
        Assertions.assertEquals(2, i2);

        T t0 = field.getZero().add(10);
        T t  = field.getZero().add(100);
        T[] complete    = MathArrays.buildArray(field, equation.getMapper().getTotalDimension());
        for (int i = 0; i < complete.length; ++i) {
            complete[i] = field.getZero().add(i);
        }
        T[] completeDot = equation.computeDerivatives(t0, complete);
        equation.init(equation.getMapper().mapStateAndDerivative(t0, complete, completeDot), t);

        T[] mainState    = equation.getMapper().extractEquationData(0,  complete);
        T[] mainStateDot = equation.getMapper().extractEquationData(0,  completeDot);
        Assertions.assertEquals(main.getDimension(), mainState.length);
        for (int i = 0; i < main.getDimension(); ++i) {
            Assertions.assertEquals(i, mainState[i].getReal(),   1.0e-15);
            Assertions.assertEquals(i, mainStateDot[i].getReal(), 1.0e-15);
            Assertions.assertEquals(i, completeDot[i].getReal(),  1.0e-15);
        }

        T[] secondaryState1    = equation.getMapper().extractEquationData(i1,  complete);
        T[] secondaryState1Dot = equation.getMapper().extractEquationData(i1,  completeDot);
        Assertions.assertEquals(secondary1.getDimension(), secondaryState1.length);
        for (int i = 0; i < secondary1.getDimension(); ++i) {
            Assertions.assertEquals(i + main.getDimension(), secondaryState1[i].getReal(),   1.0e-15);
            Assertions.assertEquals(-i, secondaryState1Dot[i].getReal(), 1.0e-15);
            Assertions.assertEquals(-i, completeDot[i + main.getDimension()].getReal(),  1.0e-15);
        }

        T[] secondaryState2    = equation.getMapper().extractEquationData(i2,  complete);
        T[] secondaryState2Dot = equation.getMapper().extractEquationData(i2,  completeDot);
        Assertions.assertEquals(secondary2.getDimension(), secondaryState2.length);
        for (int i = 0; i < secondary2.getDimension(); ++i) {
            Assertions.assertEquals(i + main.getDimension() + secondary1.getDimension(), secondaryState2[i].getReal(),   1.0e-15);
            Assertions.assertEquals(-i, secondaryState2Dot[i].getReal(), 1.0e-15);
            Assertions.assertEquals(-i, completeDot[i + main.getDimension() + secondary1.getDimension()].getReal(),  1.0e-15);
        }

    }

    @Test
    public void testMap() {
        doTestMap(Binary64Field.getInstance());
    }

    private <T extends CalculusFieldElement<T>> void doTestMap(final Field<T> field) {

        FieldOrdinaryDifferentialEquation<T> main = new Linear<T>(field, 3, 0);
        FieldExpandableODE<T> equation = new FieldExpandableODE<T>(main);
        FieldSecondaryODE<T> secondary1 = new Linear<T>(field, 3, main.getDimension());
        int i1 = equation.addSecondaryEquations(secondary1);
        FieldSecondaryODE<T> secondary2 = new Linear<T>(field, 5, main.getDimension() + secondary1.getDimension());
        int i2 = equation.addSecondaryEquations(secondary2);
        Assertions.assertEquals(main.getDimension() + secondary1.getDimension() + secondary2.getDimension(),
                            equation.getMapper().getTotalDimension());
        Assertions.assertEquals(3, equation.getMapper().getNumberOfEquations());
        Assertions.assertEquals(1, i1);
        Assertions.assertEquals(2, i2);

        T t0 = field.getZero().add(10);
        T t  = field.getZero().add(100);
        T[] complete    = MathArrays.buildArray(field, equation.getMapper().getTotalDimension());
        for (int i = 0; i < complete.length; ++i) {
            complete[i]    = field.getZero().add(i);
        }
        T[] completeDot = equation.computeDerivatives(t0, complete);
        equation.init(equation.getMapper().mapStateAndDerivative(t0, complete, completeDot), t);

        try {
            equation.getMapper().mapStateAndDerivative(t0, MathArrays.buildArray(field, complete.length + 1), completeDot);
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException dme) {
            // expected
        }
        try {
            equation.getMapper().mapStateAndDerivative(t0, complete, MathArrays.buildArray(field, completeDot.length + 1));
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException dme) {
            // expected
        }
        FieldODEStateAndDerivative<T> state = equation.getMapper().mapStateAndDerivative(t0, complete, completeDot);
        Assertions.assertEquals(2, state.getNumberOfSecondaryStates());
        Assertions.assertEquals(main.getDimension(),       state.getSecondaryStateDimension(0));
        Assertions.assertEquals(secondary1.getDimension(), state.getSecondaryStateDimension(i1));
        Assertions.assertEquals(secondary2.getDimension(), state.getSecondaryStateDimension(i2));

        T[] mainState             = state.getPrimaryState();
        T[] mainStateDot          = state.getPrimaryDerivative();
        T[] mainStateAlternate    = state.getSecondaryState(0);
        T[] mainStateDotAlternate = state.getSecondaryDerivative(0);
        Assertions.assertEquals(main.getDimension(), mainState.length);
        for (int i = 0; i < main.getDimension(); ++i) {
            Assertions.assertEquals(i, mainState[i].getReal(),             1.0e-15);
            Assertions.assertEquals(i, mainStateDot[i].getReal(),          1.0e-15);
            Assertions.assertEquals(i, mainStateAlternate[i].getReal(),    1.0e-15);
            Assertions.assertEquals(i, mainStateDotAlternate[i].getReal(), 1.0e-15);
            Assertions.assertEquals(i, completeDot[i].getReal(),           1.0e-15);
        }

        T[] secondaryState1    = state.getSecondaryState(i1);
        T[] secondaryState1Dot = state.getSecondaryDerivative(i1);
        Assertions.assertEquals(secondary1.getDimension(), secondaryState1.length);
        for (int i = 0; i < secondary1.getDimension(); ++i) {
            Assertions.assertEquals(i + main.getDimension(), secondaryState1[i].getReal(),   1.0e-15);
            Assertions.assertEquals(-i, secondaryState1Dot[i].getReal(), 1.0e-15);
            Assertions.assertEquals(-i, completeDot[i + main.getDimension()].getReal(),  1.0e-15);
        }

        T[] secondaryState2    = state.getSecondaryState(i2);
        T[] secondaryState2Dot = state.getSecondaryDerivative(i2);
        Assertions.assertEquals(secondary2.getDimension(), secondaryState2.length);
        for (int i = 0; i < secondary2.getDimension(); ++i) {
            Assertions.assertEquals(i + main.getDimension() + secondary1.getDimension(), secondaryState2[i].getReal(),   1.0e-15);
            Assertions.assertEquals(-i, secondaryState2Dot[i].getReal(), 1.0e-15);
            Assertions.assertEquals(-i, completeDot[i + main.getDimension() + secondary1.getDimension()].getReal(),  1.0e-15);
        }

        T[] remappedState = state.getCompleteState();
        T[] remappedDerivative = state.getCompleteDerivative();
        Assertions.assertEquals(equation.getMapper().getTotalDimension(), remappedState.length);
        Assertions.assertEquals(equation.getMapper().getTotalDimension(), remappedDerivative.length);
        for (int i = 0; i < remappedState.length; ++i) {
            Assertions.assertEquals(complete[i].getReal(),    remappedState[i].getReal(),      1.0e-15);
            Assertions.assertEquals(completeDot[i].getReal(), remappedDerivative[i].getReal(), 1.0e-15);
        }
    }

    @Test
    public void testExtractDimensionMismatch() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            doTestExtractDimensionMismatch(Binary64Field.getInstance());
        });
    }

    private <T extends CalculusFieldElement<T>> void doTestExtractDimensionMismatch(final Field<T> field)
        throws MathIllegalArgumentException {

        FieldOrdinaryDifferentialEquation<T> main = new Linear<T>(field, 3, 0);
        FieldExpandableODE<T> equation = new FieldExpandableODE<T>(main);
        FieldSecondaryODE<T> secondary1 = new Linear<T>(field, 3, main.getDimension());
        int i1 = equation.addSecondaryEquations(secondary1);
        T[] tooShort    = MathArrays.buildArray(field, main.getDimension());
        equation.getMapper().extractEquationData(i1, tooShort);
    }

    @Test
    public void testInsertTooShortComplete() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            doTestInsertTooShortComplete(Binary64Field.getInstance());
        });
    }

    private <T extends CalculusFieldElement<T>> void doTestInsertTooShortComplete(final Field<T> field)
        throws MathIllegalArgumentException {

        FieldOrdinaryDifferentialEquation<T> main = new Linear<T>(field, 3, 0);
        FieldExpandableODE<T> equation = new FieldExpandableODE<T>(main);
        FieldSecondaryODE<T> secondary1 = new Linear<T>(field, 3, main.getDimension());
        int i1 = equation.addSecondaryEquations(secondary1);
        T[] equationData = MathArrays.buildArray(field, secondary1.getDimension());
        T[] tooShort     = MathArrays.buildArray(field, main.getDimension());
        equation.getMapper().insertEquationData(i1, equationData, tooShort);
    }

    @Test
    public void testInsertWrongEquationData() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            doTestInsertWrongEquationData(Binary64Field.getInstance());
        });
    }

    private <T extends CalculusFieldElement<T>> void doTestInsertWrongEquationData(final Field<T> field)
        throws MathIllegalArgumentException {

        FieldOrdinaryDifferentialEquation<T> main = new Linear<T>(field, 3, 0);
        FieldExpandableODE<T> equation = new FieldExpandableODE<T>(main);
        FieldSecondaryODE<T> secondary1 = new Linear<T>(field, 3, main.getDimension());
        int i1 = equation.addSecondaryEquations(secondary1);
        T[] wrongEquationData = MathArrays.buildArray(field, secondary1.getDimension() + 1);
        T[] complete          = MathArrays.buildArray(field, equation.getMapper().getTotalDimension());
        equation.getMapper().insertEquationData(i1, wrongEquationData, complete);
    }

    @Test
    public void testNegativeIndex() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            doTestNegativeIndex(Binary64Field.getInstance());
        });
    }

    private <T extends CalculusFieldElement<T>> void doTestNegativeIndex(final Field<T> field)
        throws MathIllegalArgumentException {

        FieldOrdinaryDifferentialEquation<T> main = new Linear<T>(field, 3, 0);
        FieldExpandableODE<T> equation = new FieldExpandableODE<T>(main);
        T[] complete = MathArrays.buildArray(field, equation.getMapper().getTotalDimension());
        equation.getMapper().extractEquationData(-1, complete);
    }

    @Test
    public void testTooLargeIndex() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            doTestTooLargeIndex(Binary64Field.getInstance());
        });
    }

    private <T extends CalculusFieldElement<T>> void doTestTooLargeIndex(final Field<T> field)
        throws MathIllegalArgumentException {

        FieldOrdinaryDifferentialEquation<T> main = new Linear<T>(field, 3, 0);
        FieldExpandableODE<T> equation = new FieldExpandableODE<T>(main);
        T[] complete = MathArrays.buildArray(field, equation.getMapper().getTotalDimension());
        equation.getMapper().extractEquationData(+1, complete);
    }

    private static class  Linear<T extends CalculusFieldElement<T>>
        implements  FieldOrdinaryDifferentialEquation<T>, FieldSecondaryODE<T> {

        private final Field<T> field;
        private final int dimension;
        private final int start;

        private Linear(final Field<T> field, final int dimension, final int start) {
            this.field     = field;
            this.dimension = dimension;
            this.start     = start;
        }

        public int getDimension() {
            return dimension;
        }

        public void init(final T t0, final T[] y0, final T finalTime) {
            Assertions.assertEquals(dimension, y0.length);
            Assertions.assertEquals(10.0,  t0.getReal(), 1.0e-15);
            Assertions.assertEquals(100.0, finalTime.getReal(), 1.0e-15);
            for (int i = 0; i < y0.length; ++i) {
                Assertions.assertEquals(i, y0[i].getReal(), 1.0e-15);
            }
        }

        public T[] computeDerivatives(final T t, final T[] y) {
            final T[] yDot = MathArrays.buildArray(field, dimension);
            for (int i = 0; i < dimension; ++i) {
                yDot[i] = field.getZero().add(i);
            }
            return yDot;
        }

        public void init(final T t0, final T[] primary0, final T[] secondary0, final T finalTime) {
            Assertions.assertEquals(dimension, secondary0.length);
            Assertions.assertEquals(10.0,  t0.getReal(), 1.0e-15);
            Assertions.assertEquals(100.0, finalTime.getReal(), 1.0e-15);
            for (int i = 0; i < primary0.length; ++i) {
                Assertions.assertEquals(i, primary0[i].getReal(), 1.0e-15);
            }
            for (int i = 0; i < secondary0.length; ++i) {
                Assertions.assertEquals(start + i, secondary0[i].getReal(), 1.0e-15);
            }
        }

        public T[] computeDerivatives(final T t, final T[] primary, final T[] primaryDot, final T[] secondary) {
            final T[] secondaryDot = MathArrays.buildArray(field, dimension);
            for (int i = 0; i < dimension; ++i) {
                secondaryDot[i] = field.getZero().subtract(i);
            }
            return secondaryDot;
        }

    }

}
