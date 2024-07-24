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

import org.hipparchus.exception.MathIllegalArgumentException;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.fail;

class ExpandableODETest {

    @Test
    void testOnlyMainEquation() {
        OrdinaryDifferentialEquation main = new Linear(3, 0);
        ExpandableODE equation = new ExpandableODE(main);
        assertEquals(main.getDimension(), equation.getMapper().getTotalDimension());
        assertEquals(1, equation.getMapper().getNumberOfEquations());
        double t0 = 10;
        double t  = 100;
        double[] complete    = new double[equation.getMapper().getTotalDimension()];
        for (int i = 0; i < complete.length; ++i) {
            complete[i] = i;
        }
        double[] completeDot = equation.computeDerivatives(t0, complete);
        equation.init(equation.getMapper().mapStateAndDerivative(t0, complete, completeDot), t);
        ODEStateAndDerivative state = equation.getMapper().mapStateAndDerivative(t0, complete, completeDot);
        assertEquals(0, state.getNumberOfSecondaryStates());
        double[] mainState    = state.getPrimaryState();
        double[] mainStateDot = state.getPrimaryDerivative();
        assertEquals(main.getDimension(), mainState.length);
        for (int i = 0; i < main.getDimension(); ++i) {
            assertEquals(i, mainState[i],   1.0e-15);
            assertEquals(i, mainStateDot[i], 1.0e-15);
            assertEquals(i, completeDot[i],  1.0e-15);
        }
    }

    @Test
    void testPrimaryAndSecondary() {
        OrdinaryDifferentialEquation main = new Linear(3, 0);
        ExpandableODE equation = new ExpandableODE(main);
        SecondaryODE secondary1 = new Linear(3, main.getDimension());
        int i1 = equation.addSecondaryEquations(secondary1);
        SecondaryODE secondary2 = new Linear(5, main.getDimension() + secondary1.getDimension());
        int i2 = equation.addSecondaryEquations(secondary2);
        assertEquals(main.getDimension() + secondary1.getDimension() + secondary2.getDimension(),
                            equation.getMapper().getTotalDimension());
        assertEquals(3, equation.getMapper().getNumberOfEquations());
        assertEquals(1, i1);
        assertEquals(2, i2);

        double t0 = 10;
        double t  = 100;
        double[] complete    = new double[equation.getMapper().getTotalDimension()];
        for (int i = 0; i < complete.length; ++i) {
            complete[i] = i;
        }
        double[] completeDot = equation.computeDerivatives(t0, complete);
        equation.init(equation.getMapper().mapStateAndDerivative(t0, complete, completeDot), t);

        double[] mainState    = equation.getMapper().extractEquationData(0,  complete);
        double[] mainStateDot = equation.getMapper().extractEquationData(0,  completeDot);
        assertEquals(main.getDimension(), mainState.length);
        for (int i = 0; i < main.getDimension(); ++i) {
            assertEquals(i, mainState[i],   1.0e-15);
            assertEquals(i, mainStateDot[i], 1.0e-15);
            assertEquals(i, completeDot[i],  1.0e-15);
        }

        double[] secondaryState1    = equation.getMapper().extractEquationData(i1,  complete);
        double[] secondaryState1Dot = equation.getMapper().extractEquationData(i1,  completeDot);
        assertEquals(secondary1.getDimension(), secondaryState1.length);
        for (int i = 0; i < secondary1.getDimension(); ++i) {
            assertEquals(i + main.getDimension(), secondaryState1[i],   1.0e-15);
            assertEquals(-i, secondaryState1Dot[i], 1.0e-15);
            assertEquals(-i, completeDot[i + main.getDimension()],  1.0e-15);
        }

        double[] secondaryState2    = equation.getMapper().extractEquationData(i2,  complete);
        double[] secondaryState2Dot = equation.getMapper().extractEquationData(i2,  completeDot);
        assertEquals(secondary2.getDimension(), secondaryState2.length);
        for (int i = 0; i < secondary2.getDimension(); ++i) {
            assertEquals(i + main.getDimension() + secondary1.getDimension(), secondaryState2[i],   1.0e-15);
            assertEquals(-i, secondaryState2Dot[i], 1.0e-15);
            assertEquals(-i, completeDot[i + main.getDimension() + secondary1.getDimension()],  1.0e-15);
        }

    }

    @Test
    void testMap() {
        OrdinaryDifferentialEquation main = new Linear(3, 0);
        ExpandableODE equation = new ExpandableODE(main);
        SecondaryODE secondary1 = new Linear(3, main.getDimension());
        int i1 = equation.addSecondaryEquations(secondary1);
        SecondaryODE secondary2 = new Linear(5, main.getDimension() + secondary1.getDimension());
        int i2 = equation.addSecondaryEquations(secondary2);
        assertEquals(main.getDimension() + secondary1.getDimension() + secondary2.getDimension(),
                            equation.getMapper().getTotalDimension());
        assertEquals(3, equation.getMapper().getNumberOfEquations());
        assertEquals(1, i1);
        assertEquals(2, i2);

        double t0 = 10;
        double t  = 100;
        double[] complete    = new double[equation.getMapper().getTotalDimension()];
        for (int i = 0; i < complete.length; ++i) {
            complete[i]    = i;
        }
        double[] completeDot = equation.computeDerivatives(t0, complete);
        equation.init(equation.getMapper().mapStateAndDerivative(t0, complete, completeDot), t);

        try {
            equation.getMapper().mapStateAndDerivative(t0, new double[complete.length + 1], completeDot);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException dme) {
            // expected
        }
        try {
            equation.getMapper().mapStateAndDerivative(t0, complete, new double[completeDot.length + 1]);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException dme) {
            // expected
        }
        ODEStateAndDerivative state = equation.getMapper().mapStateAndDerivative(t0, complete, completeDot);
        assertEquals(2, state.getNumberOfSecondaryStates());
        assertEquals(main.getDimension(),       state.getSecondaryStateDimension(0));
        assertEquals(secondary1.getDimension(), state.getSecondaryStateDimension(i1));
        assertEquals(secondary2.getDimension(), state.getSecondaryStateDimension(i2));

        double[] mainState             = state.getPrimaryState();
        double[] mainStateDot          = state.getPrimaryDerivative();
        double[] mainStateAlternate    = state.getSecondaryState(0);
        double[] mainStateDotAlternate = state.getSecondaryDerivative(0);
        assertEquals(main.getDimension(), mainState.length);
        for (int i = 0; i < main.getDimension(); ++i) {
            assertEquals(i, mainState[i],             1.0e-15);
            assertEquals(i, mainStateDot[i],          1.0e-15);
            assertEquals(i, mainStateAlternate[i],    1.0e-15);
            assertEquals(i, mainStateDotAlternate[i], 1.0e-15);
            assertEquals(i, completeDot[i],           1.0e-15);
        }

        double[] secondaryState1    = state.getSecondaryState(i1);
        double[] secondaryState1Dot = state.getSecondaryDerivative(i1);
        assertEquals(secondary1.getDimension(), secondaryState1.length);
        for (int i = 0; i < secondary1.getDimension(); ++i) {
            assertEquals(i + main.getDimension(), secondaryState1[i],   1.0e-15);
            assertEquals(-i, secondaryState1Dot[i], 1.0e-15);
            assertEquals(-i, completeDot[i + main.getDimension()],  1.0e-15);
        }

        double[] secondaryState2    = state.getSecondaryState(i2);
        double[] secondaryState2Dot = state.getSecondaryDerivative(i2);
        assertEquals(secondary2.getDimension(), secondaryState2.length);
        for (int i = 0; i < secondary2.getDimension(); ++i) {
            assertEquals(i + main.getDimension() + secondary1.getDimension(), secondaryState2[i],   1.0e-15);
            assertEquals(-i, secondaryState2Dot[i], 1.0e-15);
            assertEquals(-i, completeDot[i + main.getDimension() + secondary1.getDimension()],  1.0e-15);
        }

        double[] remappedState = state.getCompleteState();
        double[] remappedDerivative = state.getCompleteDerivative();
        assertEquals(equation.getMapper().getTotalDimension(), remappedState.length);
        assertEquals(equation.getMapper().getTotalDimension(), remappedDerivative.length);
        for (int i = 0; i < remappedState.length; ++i) {
            assertEquals(complete[i],    remappedState[i],      1.0e-15);
            assertEquals(completeDot[i], remappedDerivative[i], 1.0e-15);
        }
    }

    @Test
    void testExtractDimensionMismatch() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            OrdinaryDifferentialEquation main = new Linear(3, 0);
            ExpandableODE equation = new ExpandableODE(main);
            SecondaryODE secondary1 = new Linear(3, main.getDimension());
            int i1 = equation.addSecondaryEquations(secondary1);
            double[] tooShort = new double[main.getDimension()];
            equation.getMapper().extractEquationData(i1, tooShort);
        });
    }

    @Test
    void testInsertTooShortComplete() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            OrdinaryDifferentialEquation main = new Linear(3, 0);
            ExpandableODE equation = new ExpandableODE(main);
            SecondaryODE secondary1 = new Linear(3, main.getDimension());
            int i1 = equation.addSecondaryEquations(secondary1);
            double[] equationData = new double[secondary1.getDimension()];
            double[] tooShort = new double[main.getDimension()];
            equation.getMapper().insertEquationData(i1, equationData, tooShort);
        });
    }

    @Test
    void testInsertWrongEquationData() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            OrdinaryDifferentialEquation main = new Linear(3, 0);
            ExpandableODE equation = new ExpandableODE(main);
            SecondaryODE secondary1 = new Linear(3, main.getDimension());
            int i1 = equation.addSecondaryEquations(secondary1);
            double[] wrongEquationData = new double[secondary1.getDimension() + 1];
            double[] complete = new double[equation.getMapper().getTotalDimension()];
            equation.getMapper().insertEquationData(i1, wrongEquationData, complete);
        });
    }

    @Test
    void testNegativeIndex() {
        assertThrows(MathIllegalArgumentException.class, () -> {

            OrdinaryDifferentialEquation main = new Linear(3, 0);
            ExpandableODE equation = new ExpandableODE(main);
            double[] complete = new double[equation.getMapper().getTotalDimension()];
            equation.getMapper().extractEquationData(-1, complete);
        });
    }

    @Test
    void testTooLargeIndex() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            OrdinaryDifferentialEquation main = new Linear(3, 0);
            ExpandableODE equation = new ExpandableODE(main);
            double[] complete = new double[equation.getMapper().getTotalDimension()];
            equation.getMapper().extractEquationData(+1, complete);
        });
    }

    private static class  Linear implements  OrdinaryDifferentialEquation, SecondaryODE {

        private final int dimension;
        private final int start;

        private Linear(final int dimension, final int start) {
            this.dimension = dimension;
            this.start     = start;
        }

        public int getDimension() {
            return dimension;
        }

        public void init(final double t0, final double[] y0, final double finalTime) {
            assertEquals(dimension, y0.length);
            assertEquals(10.0,  t0, 1.0e-15);
            assertEquals(100.0, finalTime, 1.0e-15);
            for (int i = 0; i < y0.length; ++i) {
                assertEquals(i, y0[i], 1.0e-15);
            }
        }

        public double[] computeDerivatives(final double t, final double[] y) {
            final double[] yDot = new double[dimension];
            for (int i = 0; i < dimension; ++i) {
                yDot[i] = i;
            }
            return yDot;
        }

        public void init(final double t0, final double[] primary0, final double[] secondary0, final double finalTime) {
            assertEquals(dimension, secondary0.length);
            assertEquals(10.0,  t0, 1.0e-15);
            assertEquals(100.0, finalTime, 1.0e-15);
            for (int i = 0; i < primary0.length; ++i) {
                assertEquals(i, primary0[i], 1.0e-15);
            }
            for (int i = 0; i < secondary0.length; ++i) {
                assertEquals(start + i, secondary0[i], 1.0e-15);
            }
        }

        public double[] computeDerivatives(final double t, final double[] primary, final double[] primaryDot, final double[] secondary) {
            final double[] secondaryDot = new double[dimension];
            for (int i = 0; i < dimension; ++i) {
                secondaryDot[i] = -i;
            }
            return secondaryDot;
        }

    }

}
