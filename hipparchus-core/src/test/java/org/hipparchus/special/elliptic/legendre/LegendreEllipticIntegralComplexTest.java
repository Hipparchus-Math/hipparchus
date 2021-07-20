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
package org.hipparchus.special.elliptic.legendre;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.analysis.CalculusFieldUnivariateFunction;
import org.hipparchus.analysis.integration.IterativeLegendreGaussIntegrator;
import org.hipparchus.complex.Complex;
import org.hipparchus.complex.ComplexUnivariateIntegrator;
import org.junit.Test;

public class LegendreEllipticIntegralComplexTest extends LegendreEllipticIntegralAbstractComplexTest<Complex> {

    protected Complex buildComplex(double realPart) {
        return new Complex(realPart);
    }

    protected Complex buildComplex(double realPart, double imaginaryPart) {
        return new Complex(realPart, imaginaryPart);
    }

    protected Complex K(Complex k) {
        return LegendreEllipticIntegral.bigK(k);
    }

    protected Complex Kprime(Complex k) {
        return LegendreEllipticIntegral.bigKPrime(k);
    }

    protected Complex F(Complex phi, Complex k) {
        return LegendreEllipticIntegral.bigF(phi, k);
    }

    protected Complex E(Complex k) {
        return LegendreEllipticIntegral.bigE(k);
    }

    protected Complex E(Complex phi, Complex k) {
        return LegendreEllipticIntegral.bigE(phi, k);
    }

    protected Complex D(Complex k) {
        return LegendreEllipticIntegral.bigD(k);
    }

    protected Complex D(Complex phi, Complex k) {
        return LegendreEllipticIntegral.bigD(phi, k);
    }

    protected Complex Pi(Complex alpha2, Complex k) {
        return LegendreEllipticIntegral.bigPi(alpha2, k);
    }

    protected Complex Pi(Complex phi, Complex alpha2, Complex k) {
        return LegendreEllipticIntegral.bigPi(phi, alpha2, k);
    }

    @Test
    public void testIssueIncompleteSecondKindA() {
        final Complex phi = new Complex(1.2, 0.75);
        final Complex m   = new Complex(0.2, 0.6);
        final Complex ref = new Complex(1.4103674846223375296500, 0.644849758860533700396);
        final ComplexUnivariateIntegrator integrator =
                        new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                             1.0e-4,
                                                                                             1.0e-4));
        final CalculusFieldUnivariateFunction<Complex> f =
                        theta -> theta.sin().multiply(theta.sin()).multiply(m).negate().add(1).sqrt();
        final Complex integrated = integrator.integrate(100000, f, new Complex(1.0e-10, 1.0e-10), phi);
        UnitTestUtils.assertEquals(ref, integrated, 2.0e-10);
        UnitTestUtils.assertEquals(ref, LegendreEllipticIntegral.bigE(phi, m), 1.0e-10);
    }

    @Test
    public void testIssueIncompleteSecondKindB() {
        final Complex phi = new Complex(1.2, 0.0);
        final Complex m   = new Complex(2.3, -1.5);
        final Complex ref = new Complex(0.8591316843513079270009549421, 0.55423174445992167002660);
        final ComplexUnivariateIntegrator integrator =
                        new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                             1.0e-4,
                                                                                             1.0e-4));
        final CalculusFieldUnivariateFunction<Complex> f =
                        theta -> theta.sin().multiply(theta.sin()).multiply(m).negate().add(1).sqrt();
        final Complex integrated = integrator.integrate(100000, f, new Complex(1.0e-10, 1.0e-10), phi);
        UnitTestUtils.assertEquals(ref, integrated, 2.0e-10);
        UnitTestUtils.assertEquals(ref, LegendreEllipticIntegral.bigE(phi, m), 1.0e-10);
    }

    @Test
    public void testIssueIncompleteSecondKindC() {
        final Complex phi = new Complex(3, 2.5);
        final Complex m   = new Complex(2.3, -1.5);
        final Complex ref = new Complex(3.05969360032192938798, 11.16503469114870865999);
        final ComplexUnivariateIntegrator integrator =
                        new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                             1.0e-4,
                                                                                             1.0e-4));
        final CalculusFieldUnivariateFunction<Complex> f =
                        theta -> theta.sin().multiply(theta.sin()).multiply(m).negate().add(1).sqrt();
        // we have to use a specific path to get the correct result
        // integrating over a single straight line gives a completely wrong result
        final Complex integrated = integrator.integrate(100000, f, new Complex(1.0e-12, 1.0e-12),
                                                        new Complex(0, -1.5), phi);
        UnitTestUtils.assertEquals(ref, integrated, 2.0e-10);
        UnitTestUtils.assertEquals(ref, LegendreEllipticIntegral.bigE(phi, m), 1.0e-10);
    }

    @Test
    public void testIssueIncompleteSecondKindD() {
        final Complex phi = new Complex(-0.4, 2.5);
        final Complex m   = new Complex(2.3, -1.5);
        final Complex ref = new Complex(-1.68645030068870706703580773597, 9.176675281683098106653799);
        final ComplexUnivariateIntegrator integrator =
                        new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                             1.0e-4,
                                                                                             1.0e-4));
        final CalculusFieldUnivariateFunction<Complex> f =
                        theta -> theta.sin().multiply(theta.sin()).multiply(m).negate().add(1).sqrt();
        final Complex integrated = integrator.integrate(100000, f, new Complex(1.0e-10, 1.0e-10), phi);
        UnitTestUtils.assertEquals(ref, integrated, 2.0e-10);
        UnitTestUtils.assertEquals(ref, LegendreEllipticIntegral.bigE(phi, m), 1.0e-10);
    }

}
