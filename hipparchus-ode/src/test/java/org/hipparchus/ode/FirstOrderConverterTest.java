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
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.nonstiff.ClassicalRungeKuttaIntegrator;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;


class FirstOrderConverterTest {

    @Test
    void testDoubleDimension() {
        for (int i = 1; i < 10; ++i) {
            SecondOrderODE eqn2 = new Equations(i, 0.2);
            FirstOrderConverter eqn1 = new FirstOrderConverter(eqn2);
            assertEquals(eqn1.getDimension(), (2 * eqn2.getDimension()));
        }
    }

    @Test
    void testDecreasingSteps()
        throws MathIllegalArgumentException, MathIllegalStateException {

        double previousError = Double.NaN;
        for (int i = 0; i < 10; ++i) {

            double step  = FastMath.pow(2.0, -(i + 1));
            double error = integrateWithSpecifiedStep(4.0, 0.0, 1.0, step)
                            - FastMath.sin(4.0);
            if (i > 0) {
                assertTrue(FastMath.abs(error) < FastMath.abs(previousError));
            }
            previousError = error;

        }
    }

    @Test
    void testSmallStep()
        throws MathIllegalArgumentException, MathIllegalStateException {
        double error = integrateWithSpecifiedStep(4.0, 0.0, 1.0, 1.0e-4)
                        - FastMath.sin(4.0);
        assertTrue(FastMath.abs(error) < 1.0e-10);
    }

    @Test
    void testBigStep()
        throws MathIllegalArgumentException, MathIllegalStateException {
        double error = integrateWithSpecifiedStep(4.0, 0.0, 1.0, 0.5)
                        - FastMath.sin(4.0);
        assertTrue(FastMath.abs(error) > 0.1);
    }

    private static class Equations
    implements SecondOrderODE {

        private int n;

        private double omega2;

        public Equations(int n, double omega) {
            this.n = n;
            omega2 = omega * omega;
        }

        public int getDimension() {
            return n;
        }

        public double[] computeSecondDerivatives(double t, double[] y, double[] yDot) {
            final double[] yDDot = new double[n];
            for (int i = 0; i < n; ++i) {
                yDDot[i] = -omega2 * y[i];
            }
            return yDDot;
        }

    }

    private double integrateWithSpecifiedStep(double omega,
                                              double t0, double t,
                                              double step) throws MathIllegalArgumentException, MathIllegalStateException {
        double[] y0 = new double[2];
        y0[0] = FastMath.sin(omega * t0);
        y0[1] = omega * FastMath.cos(omega * t0);
        ClassicalRungeKuttaIntegrator i = new ClassicalRungeKuttaIntegrator(step);
        final ODEStateAndDerivative finalstate =
                        i.integrate(new FirstOrderConverter(new Equations(1, omega)), new ODEState(t0, y0), t);
        return finalstate.getPrimaryState()[0];
    }

}
