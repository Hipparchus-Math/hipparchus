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
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertThrows;

public class ExpandableODETest {

    @Test
    public void testOnlyMainEquation() {
        OrdinaryDifferentialEquation main = new Linear(3, 0);
        ExpandableODE equation = new ExpandableODE(main);
        Assertions.assertEquals(main.getDimension(), equation.getMapper().getTotalDimension());
        Assertions.assertEquals(1, equation.getMapper().getNumberOfEquations());
        double t0 = 10;
        double t  = 100;
        double[] complete    = new double[equation.getMapper().getTotalDimension()];
        for (int i = 0; i < complete.length; ++i) {
            complete[i] = i;
        }
        double[] completeDot = equation.computeDerivatives(t0, complete);
        equation.init(equation.getMapper().mapStateAndDerivative(t0, complete, completeDot), t);
        ODEStateAndDerivative state = equation.getMapper().mapStateAndDerivative(t0, complete, completeDot);
        Assertions.assertEquals(0, state.getNumberOfSecondaryStates());
        double[] mainState    = state.getPrimaryState();
        double[] mainStateDot = state.getPrimaryDerivative();
        Assertions.assertEquals(main.getDimension(), mainState.length);
        for (int i = 0; i < main.getDimension(); ++i) {
            Assertions.assertEquals(i, mainState[i],   1.0e-15);
            Assertions.assertEquals(i, mainStateDot[i], 1.0e-15);
            Assertions.assertEquals(i, completeDot[i],  1.0e-15);
        }
    }

    @Test
    public void testPrimaryAndSecondary() {
        OrdinaryDifferentialEquation main = new Linear(3, 0);
        ExpandableODE equation = new ExpandableODE(main);
        SecondaryODE secondary1 = new Linear(3, main.getDimension());
        int i1 = equation.addSecondaryEquations(secondary1);
        SecondaryODE secondary2 = new Linear(5, main.getDimension() + secondary1.getDimension());
        int i2 = equation.addSecondaryEquations(secondary2);
        Assertions.assertEquals(main.getDimension() + secondary1.getDimension() + secondary2.getDimension(),
                            equation.getMapper().getTotalDimension());
        Assertions.assertEquals(3, equation.getMapper().getNumberOfEquations());
        Assertions.assertEquals(1, i1);
        Assertions.assertEquals(2, i2);

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
        Assertions.assertEquals(main.getDimension(), mainState.length);
        for (int i = 0; i < main.getDimension(); ++i) {
            Assertions.assertEquals(i, mainState[i],   1.0e-15);
            Assertions.assertEquals(i, mainStateDot[i], 1.0e-15);
            Assertions.assertEquals(i, completeDot[i],  1.0e-15);
        }

        double[] secondaryState1    = equation.getMapper().extractEquationData(i1,  complete);
        double[] secondaryState1Dot = equation.getMapper().extractEquationData(i1,  completeDot);
        Assertions.assertEquals(secondary1.getDimension(), secondaryState1.length);
        for (int i = 0; i < secondary1.getDimension(); ++i) {
            Assertions.assertEquals(i + main.getDimension(), secondaryState1[i],   1.0e-15);
            Assertions.assertEquals(-i, secondaryState1Dot[i], 1.0e-15);
            Assertions.assertEquals(-i, completeDot[i + main.getDimension()],  1.0e-15);
        }

        double[] secondaryState2    = equation.getMapper().extractEquationData(i2,  complete);
        double[] secondaryState2Dot = equation.getMapper().extractEquationData(i2,  completeDot);
        Assertions.assertEquals(secondary2.getDimension(), secondaryState2.length);
        for (int i = 0; i < secondary2.getDimension(); ++i) {
            Assertions.assertEquals(i + main.getDimension() + secondary1.getDimension(), secondaryState2[i],   1.0e-15);
            Assertions.assertEquals(-i, secondaryState2Dot[i], 1.0e-15);
            Assertions.assertEquals(-i, completeDot[i + main.getDimension() + secondary1.getDimension()],  1.0e-15);
        }

    }

    @Test
    public void testMap() {
        OrdinaryDifferentialEquation main = new Linear(3, 0);
        ExpandableODE equation = new ExpandableODE(main);
        SecondaryODE secondary1 = new Linear(3, main.getDimension());
        int i1 = equation.addSecondaryEquations(secondary1);
        SecondaryODE secondary2 = new Linear(5, main.getDimension() + secondary1.getDimension());
        int i2 = equation.addSecondaryEquations(secondary2);
        Assertions.assertEquals(main.getDimension() + secondary1.getDimension() + secondary2.getDimension(),
                            equation.getMapper().getTotalDimension());
        Assertions.assertEquals(3, equation.getMapper().getNumberOfEquations());
        Assertions.assertEquals(1, i1);
        Assertions.assertEquals(2, i2);

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
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException dme) {
            // expected
        }
        try {
            equation.getMapper().mapStateAndDerivative(t0, complete, new double[completeDot.length + 1]);
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException dme) {
            // expected
        }
        ODEStateAndDerivative state = equation.getMapper().mapStateAndDerivative(t0, complete, completeDot);
        Assertions.assertEquals(2, state.getNumberOfSecondaryStates());
        Assertions.assertEquals(main.getDimension(),       state.getSecondaryStateDimension(0));
        Assertions.assertEquals(secondary1.getDimension(), state.getSecondaryStateDimension(i1));
        Assertions.assertEquals(secondary2.getDimension(), state.getSecondaryStateDimension(i2));

        double[] mainState             = state.getPrimaryState();
        double[] mainStateDot          = state.getPrimaryDerivative();
        double[] mainStateAlternate    = state.getSecondaryState(0);
        double[] mainStateDotAlternate = state.getSecondaryDerivative(0);
        Assertions.assertEquals(main.getDimension(), mainState.length);
        for (int i = 0; i < main.getDimension(); ++i) {
            Assertions.assertEquals(i, mainState[i],             1.0e-15);
            Assertions.assertEquals(i, mainStateDot[i],          1.0e-15);
            Assertions.assertEquals(i, mainStateAlternate[i],    1.0e-15);
            Assertions.assertEquals(i, mainStateDotAlternate[i], 1.0e-15);
            Assertions.assertEquals(i, completeDot[i],           1.0e-15);
        }

        double[] secondaryState1    = state.getSecondaryState(i1);
        double[] secondaryState1Dot = state.getSecondaryDerivative(i1);
        Assertions.assertEquals(secondary1.getDimension(), secondaryState1.length);
        for (int i = 0; i < secondary1.getDimension(); ++i) {
            Assertions.assertEquals(i + main.getDimension(), secondaryState1[i],   1.0e-15);
            Assertions.assertEquals(-i, secondaryState1Dot[i], 1.0e-15);
            Assertions.assertEquals(-i, completeDot[i + main.getDimension()],  1.0e-15);
        }

        double[] secondaryState2    = state.getSecondaryState(i2);
        double[] secondaryState2Dot = state.getSecondaryDerivative(i2);
        Assertions.assertEquals(secondary2.getDimension(), secondaryState2.length);
        for (int i = 0; i < secondary2.getDimension(); ++i) {
            Assertions.assertEquals(i + main.getDimension() + secondary1.getDimension(), secondaryState2[i],   1.0e-15);
            Assertions.assertEquals(-i, secondaryState2Dot[i], 1.0e-15);
            Assertions.assertEquals(-i, completeDot[i + main.getDimension() + secondary1.getDimension()],  1.0e-15);
        }

        double[] remappedState = state.getCompleteState();
        double[] remappedDerivative = state.getCompleteDerivative();
        Assertions.assertEquals(equation.getMapper().getTotalDimension(), remappedState.length);
        Assertions.assertEquals(equation.getMapper().getTotalDimension(), remappedDerivative.length);
        for (int i = 0; i < remappedState.length; ++i) {
            Assertions.assertEquals(complete[i],    remappedState[i],      1.0e-15);
            Assertions.assertEquals(completeDot[i], remappedDerivative[i], 1.0e-15);
        }
    }

    @Test
    public void testExtractDimensionMismatch() {
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
    public void testInsertTooShortComplete() {
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
    public void testInsertWrongEquationData() {
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
    public void testNegativeIndex() {
        assertThrows(MathIllegalArgumentException.class, () -> {

            OrdinaryDifferentialEquation main = new Linear(3, 0);
            ExpandableODE equation = new ExpandableODE(main);
            double[] complete = new double[equation.getMapper().getTotalDimension()];
            equation.getMapper().extractEquationData(-1, complete);
        });
    }

    @Test
    public void testTooLargeIndex() {
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
            Assertions.assertEquals(dimension, y0.length);
            Assertions.assertEquals(10.0,  t0, 1.0e-15);
            Assertions.assertEquals(100.0, finalTime, 1.0e-15);
            for (int i = 0; i < y0.length; ++i) {
                Assertions.assertEquals(i, y0[i], 1.0e-15);
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
            Assertions.assertEquals(dimension, secondary0.length);
            Assertions.assertEquals(10.0,  t0, 1.0e-15);
            Assertions.assertEquals(100.0, finalTime, 1.0e-15);
            for (int i = 0; i < primary0.length; ++i) {
                Assertions.assertEquals(i, primary0[i], 1.0e-15);
            }
            for (int i = 0; i < secondary0.length; ++i) {
                Assertions.assertEquals(start + i, secondary0[i], 1.0e-15);
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
