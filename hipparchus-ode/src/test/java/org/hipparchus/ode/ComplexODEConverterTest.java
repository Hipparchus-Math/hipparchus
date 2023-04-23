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

package org.hipparchus.ode;

import org.hipparchus.complex.Complex;
import org.hipparchus.ode.nonstiff.LutherIntegrator;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;


public class ComplexODEConverterTest {

    @Test
    public void testDoubleDimension() {
        ComplexODEConverter converter = new ComplexODEConverter();
        for (int i = 1; i < 10; ++i) {
            Assert.assertEquals(2 * i, converter.convertEquations(new Circle(i, 0.2)).getDimension());
        }
    }

    @Test
    public void testPrimaryEquation() {
        final double omega = 0.2;
        ComplexODEConverter converter = new ComplexODEConverter();
        ComplexOrdinaryDifferentialEquation circle = new Circle(1, omega);
        ComplexODEState initial = new ComplexODEState(0.0, new Complex[] { Complex.ONE });
        ComplexODEStateAndDerivative der = new ComplexODEStateAndDerivative(initial.getTime(),
                                                                            initial.getPrimaryState(),
                                                                            circle.computeDerivatives(initial.getTime(),
                                                                                                      initial.getPrimaryState()));
        Assert.assertEquals(initial.getTime(), der.getTime(), 1.0e-15);
        Assert.assertEquals(initial.getPrimaryState()[0], der.getPrimaryState()[0]);
        Assert.assertEquals(initial.getPrimaryState()[0], der.getCompleteState()[0]);
        Assert.assertEquals(initial.getPrimaryState()[0].multiply(((Circle) circle).iOmega),
                            der.getSecondaryDerivative(0)[0]);
        Assert.assertEquals(initial.getPrimaryState()[0].multiply(((Circle) circle).iOmega),
                            der.getCompleteDerivative()[0]);
        LutherIntegrator integrator = new LutherIntegrator(1.0e-3);
        final ComplexODEStateAndDerivative finalstate =
                        converter.convertState(integrator.integrate(converter.convertEquations(circle),
                                                                    converter.convertState(initial),
                                                                    FastMath.PI / omega));
        Assert.assertEquals(0, finalstate.getNumberOfSecondaryStates());
        Assert.assertEquals(1, finalstate.getPrimaryStateDimension());
        Assert.assertEquals(FastMath.PI / omega, finalstate.getTime(),                                1.0e-15);
        Assert.assertEquals(-1.0,                finalstate.getPrimaryState()[0].getReal(),           1.0e-12);
        Assert.assertEquals( 0.0,                finalstate.getPrimaryState()[0].getImaginary(),      1.0e-12);
        Assert.assertEquals( 0.0,                finalstate.getPrimaryDerivative()[0].getReal(),      1.0e-12);
        Assert.assertEquals(-omega,              finalstate.getPrimaryDerivative()[0].getImaginary(), 1.0e-12);
    }

    @Test
    public void testSecondaryEquation() {
        final double omegaP  = 0.2;
        final double omegaS0 = omegaP * 0.5;
        final double omegaS1 = omegaP * 2;
        ComplexODEConverter converter = new ComplexODEConverter();
        ComplexOrdinaryDifferentialEquation primary    = new Circle(1, omegaP);
        ComplexSecondaryODE                 secondary0 = new Circle(1, omegaS0);
        ComplexSecondaryODE                 secondary1 = new Circle(1, omegaS1);
        ExpandableODE expandable = new ExpandableODE(converter.convertEquations(primary));
        expandable.addSecondaryEquations(converter.convertSecondaryEquations(secondary0));
        expandable.addSecondaryEquations(converter.convertSecondaryEquations(secondary1));
        ComplexODEState initial = new ComplexODEState(0.0, new Complex[] { Complex.ONE },
                                                      new Complex[][] {
                                                          { Complex.ONE },
                                                          { Complex.ONE }
                                                      });
        ComplexODEStateAndDerivative der = new ComplexODEStateAndDerivative(initial.getTime(),
                                                                            initial.getPrimaryState(),
                                                                            primary.computeDerivatives(initial.getTime(),
                                                                                                      initial.getPrimaryState()),
                                                                            new Complex[][] {
                                                                                initial.getSecondaryState(0),
                                                                                initial.getSecondaryState(1)
                                                                            },
                                                                            new Complex[][] {
                                                                                secondary0.computeDerivatives(initial.getTime(),
                                                                                                              initial.getPrimaryState(),
                                                                                                              primary.computeDerivatives(initial.getTime(),
                                                                                                                                         initial.getPrimaryState()),
                                                                                                              initial.getSecondaryState(0)),
                                                                                secondary1.computeDerivatives(initial.getTime(),
                                                                                                              initial.getPrimaryState(),
                                                                                                              primary.computeDerivatives(initial.getTime(),
                                                                                                                                         initial.getPrimaryState()),
                                                                                                              initial.getSecondaryState(1))
                                                                            });
        Assert.assertEquals(initial.getTime(), der.getTime(), 1.0e-15);
        Assert.assertEquals(initial.getPrimaryState()[0], der.getPrimaryState()[0]);
        Assert.assertEquals(initial.getPrimaryState()[0], der.getCompleteState()[0]);
        Assert.assertEquals(initial.getSecondaryState(0)[0], der.getCompleteState()[1]);
        Assert.assertEquals(initial.getSecondaryState(1)[0], der.getCompleteState()[2]);
        Assert.assertEquals(initial.getPrimaryState()[0].multiply(((Circle) primary).iOmega),
                            der.getSecondaryDerivative(0)[0]);
        Assert.assertEquals(initial.getPrimaryState()[0].multiply(((Circle) primary).iOmega),
                            der.getCompleteDerivative()[0]);
        Assert.assertEquals(initial.getSecondaryState(0)[0].multiply(((Circle) secondary0).iOmega),
                            der.getCompleteDerivative()[1]);
        Assert.assertEquals(initial.getSecondaryState(1)[0].multiply(((Circle) secondary1).iOmega),
                            der.getCompleteDerivative()[2]);
        LutherIntegrator integrator = new LutherIntegrator(1.0e-3);
        final ComplexODEStateAndDerivative finalstate =
                        converter.convertState(integrator.integrate(expandable,
                                                                    converter.convertState(initial),
                                                                    FastMath.PI / omegaP));
        Assert.assertEquals(2, finalstate.getNumberOfSecondaryStates());
        Assert.assertEquals(1, finalstate.getPrimaryStateDimension());
        Assert.assertEquals(FastMath.PI / omegaP, finalstate.getTime(),                                   1.0e-15);
        Assert.assertEquals(-1.0,                 finalstate.getPrimaryState()[0].getReal(),              1.0e-12);
        Assert.assertEquals( 0.0,                 finalstate.getPrimaryState()[0].getImaginary(),         1.0e-12);
        Assert.assertEquals( 0.0,                 finalstate.getPrimaryDerivative()[0].getReal(),         1.0e-12);
        Assert.assertEquals(-omegaP,              finalstate.getPrimaryDerivative()[0].getImaginary(),    1.0e-12);
        Assert.assertEquals( 0.0,                 finalstate.getSecondaryState(1)[0].getReal(),           1.0e-12);
        Assert.assertEquals( 1.0,                 finalstate.getSecondaryState(1)[0].getImaginary(),      1.0e-12);
        Assert.assertEquals(-omegaS0,             finalstate.getSecondaryDerivative(1)[0].getReal(),      1.0e-12);
        Assert.assertEquals( 0.0,                 finalstate.getSecondaryDerivative(1)[0].getImaginary(), 1.0e-12);
        Assert.assertEquals( 1.0,                 finalstate.getSecondaryState(2)[0].getReal(),           1.0e-12);
        Assert.assertEquals( 0.0,                 finalstate.getSecondaryState(2)[0].getImaginary(),      2.0e-12);
        Assert.assertEquals( 0.0,                 finalstate.getSecondaryDerivative(2)[0].getReal(),      1.0e-12);
        Assert.assertEquals( omegaS1,             finalstate.getSecondaryDerivative(2)[0].getImaginary(), 1.0e-12);
    }

    private static class Circle
        implements ComplexOrdinaryDifferentialEquation, ComplexSecondaryODE {

        private int n;
        private Complex iOmega;

        public Circle(int n, double omega) {
            this.n     = n;
            this.iOmega = new Complex(0.0, omega);
        }

        public int getDimension() {
            return n;
        }

        public Complex[] computeDerivatives(double t, Complex[] y) {
            final Complex[] yDot = new Complex[y.length];
            for (int i = 0; i < yDot.length; ++i) {
                yDot[i] = iOmega.multiply(y[i]);
            }
            return yDot;
        }

        public Complex[] computeDerivatives(double t, Complex[] primary, Complex[] primaryDot, Complex[] secondary) {
            return computeDerivatives(t, secondary);
        }

    }

}
