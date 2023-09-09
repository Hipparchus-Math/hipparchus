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
import org.junit.Assert;
import org.junit.Test;

public class ExpandableODETest {

    @Test
    public void testOnlyMainEquation() {
        OrdinaryDifferentialEquation main = new Linear(3, 0);
        ExpandableODE equation = new ExpandableODE(main);
        Assert.assertEquals(main.getDimension(), equation.getMapper().getTotalDimension());
        Assert.assertEquals(1, equation.getMapper().getNumberOfEquations());
        double t0 = 10;
        double t  = 100;
        double[] complete    = new double[equation.getMapper().getTotalDimension()];
        for (int i = 0; i < complete.length; ++i) {
            complete[i] = i;
        }
        double[] completeDot = equation.computeDerivatives(t0, complete);
        equation.init(equation.getMapper().mapStateAndDerivative(t0, complete, completeDot), t);
        ODEStateAndDerivative state = equation.getMapper().mapStateAndDerivative(t0, complete, completeDot);
        Assert.assertEquals(0, state.getNumberOfSecondaryStates());
        double[] mainState    = state.getPrimaryState();
        double[] mainStateDot = state.getPrimaryDerivative();
        Assert.assertEquals(main.getDimension(), mainState.length);
        for (int i = 0; i < main.getDimension(); ++i) {
            Assert.assertEquals(i, mainState[i],   1.0e-15);
            Assert.assertEquals(i, mainStateDot[i], 1.0e-15);
            Assert.assertEquals(i, completeDot[i],  1.0e-15);
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
        Assert.assertEquals(main.getDimension() + secondary1.getDimension() + secondary2.getDimension(),
                            equation.getMapper().getTotalDimension());
        Assert.assertEquals(3, equation.getMapper().getNumberOfEquations());
        Assert.assertEquals(1, i1);
        Assert.assertEquals(2, i2);

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
        Assert.assertEquals(main.getDimension(), mainState.length);
        for (int i = 0; i < main.getDimension(); ++i) {
            Assert.assertEquals(i, mainState[i],   1.0e-15);
            Assert.assertEquals(i, mainStateDot[i], 1.0e-15);
            Assert.assertEquals(i, completeDot[i],  1.0e-15);
        }

        double[] secondaryState1    = equation.getMapper().extractEquationData(i1,  complete);
        double[] secondaryState1Dot = equation.getMapper().extractEquationData(i1,  completeDot);
        Assert.assertEquals(secondary1.getDimension(), secondaryState1.length);
        for (int i = 0; i < secondary1.getDimension(); ++i) {
            Assert.assertEquals(i + main.getDimension(), secondaryState1[i],   1.0e-15);
            Assert.assertEquals(-i, secondaryState1Dot[i], 1.0e-15);
            Assert.assertEquals(-i, completeDot[i + main.getDimension()],  1.0e-15);
        }

        double[] secondaryState2    = equation.getMapper().extractEquationData(i2,  complete);
        double[] secondaryState2Dot = equation.getMapper().extractEquationData(i2,  completeDot);
        Assert.assertEquals(secondary2.getDimension(), secondaryState2.length);
        for (int i = 0; i < secondary2.getDimension(); ++i) {
            Assert.assertEquals(i + main.getDimension() + secondary1.getDimension(), secondaryState2[i],   1.0e-15);
            Assert.assertEquals(-i, secondaryState2Dot[i], 1.0e-15);
            Assert.assertEquals(-i, completeDot[i + main.getDimension() + secondary1.getDimension()],  1.0e-15);
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
        Assert.assertEquals(main.getDimension() + secondary1.getDimension() + secondary2.getDimension(),
                            equation.getMapper().getTotalDimension());
        Assert.assertEquals(3, equation.getMapper().getNumberOfEquations());
        Assert.assertEquals(1, i1);
        Assert.assertEquals(2, i2);

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
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException dme) {
            // expected
        }
        try {
            equation.getMapper().mapStateAndDerivative(t0, complete, new double[completeDot.length + 1]);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException dme) {
            // expected
        }
        ODEStateAndDerivative state = equation.getMapper().mapStateAndDerivative(t0, complete, completeDot);
        Assert.assertEquals(2, state.getNumberOfSecondaryStates());
        Assert.assertEquals(main.getDimension(),       state.getSecondaryStateDimension(0));
        Assert.assertEquals(secondary1.getDimension(), state.getSecondaryStateDimension(i1));
        Assert.assertEquals(secondary2.getDimension(), state.getSecondaryStateDimension(i2));

        double[] mainState             = state.getPrimaryState();
        double[] mainStateDot          = state.getPrimaryDerivative();
        double[] mainStateAlternate    = state.getSecondaryState(0);
        double[] mainStateDotAlternate = state.getSecondaryDerivative(0);
        Assert.assertEquals(main.getDimension(), mainState.length);
        for (int i = 0; i < main.getDimension(); ++i) {
            Assert.assertEquals(i, mainState[i],             1.0e-15);
            Assert.assertEquals(i, mainStateDot[i],          1.0e-15);
            Assert.assertEquals(i, mainStateAlternate[i],    1.0e-15);
            Assert.assertEquals(i, mainStateDotAlternate[i], 1.0e-15);
            Assert.assertEquals(i, completeDot[i],           1.0e-15);
        }

        double[] secondaryState1    = state.getSecondaryState(i1);
        double[] secondaryState1Dot = state.getSecondaryDerivative(i1);
        Assert.assertEquals(secondary1.getDimension(), secondaryState1.length);
        for (int i = 0; i < secondary1.getDimension(); ++i) {
            Assert.assertEquals(i + main.getDimension(), secondaryState1[i],   1.0e-15);
            Assert.assertEquals(-i, secondaryState1Dot[i], 1.0e-15);
            Assert.assertEquals(-i, completeDot[i + main.getDimension()],  1.0e-15);
        }

        double[] secondaryState2    = state.getSecondaryState(i2);
        double[] secondaryState2Dot = state.getSecondaryDerivative(i2);
        Assert.assertEquals(secondary2.getDimension(), secondaryState2.length);
        for (int i = 0; i < secondary2.getDimension(); ++i) {
            Assert.assertEquals(i + main.getDimension() + secondary1.getDimension(), secondaryState2[i],   1.0e-15);
            Assert.assertEquals(-i, secondaryState2Dot[i], 1.0e-15);
            Assert.assertEquals(-i, completeDot[i + main.getDimension() + secondary1.getDimension()],  1.0e-15);
        }

        double[] remappedState = state.getCompleteState();
        double[] remappedDerivative = state.getCompleteDerivative();
        Assert.assertEquals(equation.getMapper().getTotalDimension(), remappedState.length);
        Assert.assertEquals(equation.getMapper().getTotalDimension(), remappedDerivative.length);
        for (int i = 0; i < remappedState.length; ++i) {
            Assert.assertEquals(complete[i],    remappedState[i],      1.0e-15);
            Assert.assertEquals(completeDot[i], remappedDerivative[i], 1.0e-15);
        }
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testExtractDimensionMismatch() {
        OrdinaryDifferentialEquation main = new Linear(3, 0);
        ExpandableODE equation = new ExpandableODE(main);
        SecondaryODE secondary1 = new Linear(3, main.getDimension());
        int i1 = equation.addSecondaryEquations(secondary1);
        double[] tooShort    = new double[main.getDimension()];
        equation.getMapper().extractEquationData(i1, tooShort);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testInsertTooShortComplete() {
        OrdinaryDifferentialEquation main = new Linear(3, 0);
        ExpandableODE equation = new ExpandableODE(main);
        SecondaryODE secondary1 = new Linear(3, main.getDimension());
        int i1 = equation.addSecondaryEquations(secondary1);
        double[] equationData = new double[secondary1.getDimension()];
        double[] tooShort     = new double[main.getDimension()];
        equation.getMapper().insertEquationData(i1, equationData, tooShort);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testInsertWrongEquationData() {
        OrdinaryDifferentialEquation main = new Linear(3, 0);
        ExpandableODE equation = new ExpandableODE(main);
        SecondaryODE secondary1 = new Linear(3, main.getDimension());
        int i1 = equation.addSecondaryEquations(secondary1);
        double[] wrongEquationData = new double[secondary1.getDimension() + 1];
        double[] complete          = new double[equation.getMapper().getTotalDimension()];
        equation.getMapper().insertEquationData(i1, wrongEquationData, complete);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testNegativeIndex() {

        OrdinaryDifferentialEquation main = new Linear(3, 0);
        ExpandableODE equation = new ExpandableODE(main);
        double[] complete = new double[equation.getMapper().getTotalDimension()];
        equation.getMapper().extractEquationData(-1, complete);
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testTooLargeIndex() {
        OrdinaryDifferentialEquation main = new Linear(3, 0);
        ExpandableODE equation = new ExpandableODE(main);
        double[] complete = new double[equation.getMapper().getTotalDimension()];
        equation.getMapper().extractEquationData(+1, complete);
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
            Assert.assertEquals(dimension, y0.length);
            Assert.assertEquals(10.0,  t0, 1.0e-15);
            Assert.assertEquals(100.0, finalTime, 1.0e-15);
            for (int i = 0; i < y0.length; ++i) {
                Assert.assertEquals(i, y0[i], 1.0e-15);
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
            Assert.assertEquals(dimension, secondary0.length);
            Assert.assertEquals(10.0,  t0, 1.0e-15);
            Assert.assertEquals(100.0, finalTime, 1.0e-15);
            for (int i = 0; i < primary0.length; ++i) {
                Assert.assertEquals(i, primary0[i], 1.0e-15);
            }
            for (int i = 0; i < secondary0.length; ++i) {
                Assert.assertEquals(start + i, secondary0[i], 1.0e-15);
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
