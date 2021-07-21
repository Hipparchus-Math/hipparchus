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
import org.junit.Ignore;
import org.junit.Test;

public class LegendreEllipticIntegralComplexTest extends LegendreEllipticIntegralAbstractComplexTest<Complex> {

    protected Complex buildComplex(double realPart) {
        return new Complex(realPart);
    }

    protected Complex buildComplex(double realPart, double imaginaryPart) {
        return new Complex(realPart, imaginaryPart);
    }

    protected Complex K(Complex m) {
        return LegendreEllipticIntegral.bigK(m);
    }

    protected Complex Kprime(Complex m) {
        return LegendreEllipticIntegral.bigKPrime(m);
    }

    protected Complex F(Complex phi, Complex m) {
        return LegendreEllipticIntegral.bigF(phi, m);
    }

    protected Complex E(Complex m) {
        return LegendreEllipticIntegral.bigE(m);
    }

    protected Complex E(Complex phi, Complex m) {
        return LegendreEllipticIntegral.bigE(phi, m);
    }

    protected Complex D(Complex m) {
        return LegendreEllipticIntegral.bigD(m);
    }

    protected Complex D(Complex phi, Complex m) {
        return LegendreEllipticIntegral.bigD(phi, m);
    }

    protected Complex Pi(Complex n, Complex m) {
        return LegendreEllipticIntegral.bigPi(n, m);
    }

    protected Complex Pi(Complex n, Complex phi, Complex m) {
        return LegendreEllipticIntegral.bigPi(n, phi, m);
    }

    @Test
    public void testIssueIncompleteDifferenceA() {
        final Complex phi = new Complex(1.2, 0.75);
        final Complex m   = new Complex(0.2, 0.6);
        final Complex ref = F(phi, m).
                            subtract(E(phi, m)).
                            divide(m);
        final ComplexUnivariateIntegrator integrator =
                        new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                             1.0e-4,
                                                                                             1.0e-4));
        final Complex integrated = integrator.integrate(100000, new Difference(m),
                                                        new Complex(1.0e-10, 1.0e-10), phi);
        UnitTestUtils.assertEquals(ref, integrated, 2.0e-10);
        UnitTestUtils.assertEquals(ref, D(phi, m), 1.0e-10);
    }

    @Test
    public void testIssueIncompleteDifferenceB() {
        final Complex phi = new Complex(1.2, 0.0);
        final Complex m   = new Complex(2.3, -1.5);
        final Complex ref = F(phi, m).
                            subtract(E(phi, m)).
                            divide(m);
        final ComplexUnivariateIntegrator integrator =
                        new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                             1.0e-4,
                                                                                             1.0e-4));
        final Complex integrated = integrator.integrate(100000, new Difference(m),
                                                        new Complex(1.0e-10, 1.0e-10), phi);
        UnitTestUtils.assertEquals(ref, integrated, 2.0e-10);
        UnitTestUtils.assertEquals(ref, D(phi, m), 1.0e-10);
    }

    @Test
    public void testIssueIncompleteDifferenceC() {
        final Complex phi = new Complex(3, 2.5);
        final Complex m   = new Complex(2.3, -1.5);
        final Complex ref = F(phi, m).
                            subtract(E(phi, m)).
                            divide(m);
        final ComplexUnivariateIntegrator integrator =
                        new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                             1.0e-4,
                                                                                             1.0e-4));
        // we have to use a specific path to get the correct result
        // integrating over a single straight line gives a completely wrong result
        final Complex integrated = integrator.integrate(100000, new Difference(m),
                                                        new Complex(1.0e-12, 1.0e-12), new Complex(0, -1.5), phi);
        UnitTestUtils.assertEquals(ref, integrated, 2.0e-10);
        UnitTestUtils.assertEquals(ref, D(phi, m), 1.0e-10);
    }

    @Test
    public void testIssueIncompleteDifferenceD() {
        final Complex phi = new Complex(-0.4, 2.5);
        final Complex m   = new Complex(2.3, -1.5);
        final Complex ref = F(phi, m).
                            subtract(E(phi, m)).
                            divide(m);
        final ComplexUnivariateIntegrator integrator =
                        new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                             1.0e-4,
                                                                                             1.0e-4));
        final Complex integrated = integrator.integrate(100000, new Difference(m),
                                                        new Complex(1.0e-10, 1.0e-10), phi);
        UnitTestUtils.assertEquals(ref, integrated, 2.0e-10);
        UnitTestUtils.assertEquals(ref, D(phi, m), 1.0e-10);
    }

    @Test
    public void testIssueIncompleteFirstKindA() {
        final Complex phi = new Complex(1.2, 0.75);
        final Complex m   = new Complex(0.2, 0.6);
        final Complex ref = new Complex(1.00265860821563927579252866, 0.80128721521822408811217);
        final ComplexUnivariateIntegrator integrator =
                        new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                             1.0e-4,
                                                                                             1.0e-4));
        final Complex integrated = integrator.integrate(100000, new First(m),
                                                        new Complex(1.0e-10, 1.0e-10), phi);
        UnitTestUtils.assertEquals(ref, integrated, 2.0e-10);
        UnitTestUtils.assertEquals(ref, F(phi, m), 1.0e-10);
    }

    @Test
    public void testIssueIncompleteFirstKindB() {
        final Complex phi = new Complex(1.2, 0.0);
        final Complex m   = new Complex(2.3, -1.5);
        final Complex ref = new Complex(1.04335840461807753156026488, -0.5872679121672512828049797);
        final ComplexUnivariateIntegrator integrator =
                        new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                             1.0e-4,
                                                                                             1.0e-4));
        final Complex integrated = integrator.integrate(100000, new First(m),
                                                        new Complex(1.0e-10, 1.0e-10), phi);
        UnitTestUtils.assertEquals(ref, integrated, 2.0e-10);
        UnitTestUtils.assertEquals(ref, F(phi, m), 1.0e-10);
    }

    @Test
    public void testIssueIncompleteFirstKindC() {
        final Complex phi = new Complex(3, 2.5);
        final Complex m   = new Complex(2.3, -1.5);
        final Complex ref = new Complex(2.13626296176181376169951646, -0.573329373615824705851275203);
        final ComplexUnivariateIntegrator integrator =
                        new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                             1.0e-4,
                                                                                             1.0e-4));
        // we have to use a specific path to get the correct result
        // integrating over a single straight line gives a completely wrong result
        final Complex integrated = integrator.integrate(100000, new First(m),
                                                        new Complex(1.0e-12, 1.0e-12), new Complex(0, -1.5), phi);
        UnitTestUtils.assertEquals(ref, integrated, 2.0e-10);
        UnitTestUtils.assertEquals(ref, F(phi, m), 1.0e-10);
    }

    @Test
    public void testIssueIncompleteFirstKindD() {
        final Complex phi = new Complex(-0.4, 2.5);
        final Complex m   = new Complex(2.3, -1.5);
        final Complex ref = new Complex(-0.20646268947416273887690961, 1.0927692344374984107332330625089);
        final ComplexUnivariateIntegrator integrator =
                        new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                             1.0e-4,
                                                                                             1.0e-4));
        final Complex integrated = integrator.integrate(100000, new First(m),
                                                        new Complex(1.0e-10, 1.0e-10), phi);
        UnitTestUtils.assertEquals(ref, integrated, 2.0e-10);
        UnitTestUtils.assertEquals(ref, F(phi, m), 1.0e-10);
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
        final Complex integrated = integrator.integrate(100000, new Second(m),
                                                        new Complex(1.0e-10, 1.0e-10), phi);
        UnitTestUtils.assertEquals(ref, integrated, 2.0e-10);
        UnitTestUtils.assertEquals(ref, E(phi, m), 1.0e-10);
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
        final Complex integrated = integrator.integrate(100000, new Second(m),
                                                        new Complex(1.0e-10, 1.0e-10), phi);
        UnitTestUtils.assertEquals(ref, integrated, 2.0e-10);
        UnitTestUtils.assertEquals(ref, E(phi, m), 1.0e-10);
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
        // we have to use a specific path to get the correct result
        // integrating over a single straight line gives a completely wrong result
        final Complex integrated = integrator.integrate(100000, new Second(m),
                                                        new Complex(1.0e-12, 1.0e-12), new Complex(0, -1.5), phi);
        UnitTestUtils.assertEquals(ref, integrated, 2.0e-10);
        UnitTestUtils.assertEquals(ref, E(phi, m), 1.0e-10);
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
        final Complex integrated = integrator.integrate(100000, new Second(m),
                                                        new Complex(1.0e-10, 1.0e-10), phi);
        UnitTestUtils.assertEquals(ref, integrated, 2.0e-10);
        UnitTestUtils.assertEquals(ref, E(phi, m), 1.0e-10);
    }

    // TODO: this test fails and the Wolfram reference is consistent with numerical integral
    //       there is no argument reduction here
    @Ignore
    @Test
    public void testIssueIncompleteThirdKindA() {
        final Complex n      = new Complex(3.4, -1.3);
        final Complex phi    = new Complex(1.2, 0.75);
        final Complex m      = new Complex(0.2, 0.6);
        final Complex ref    = new Complex(-0.0377623234559622331654504518, 0.653477244699729422566224);
        final ComplexUnivariateIntegrator integrator =
                        new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                             1.0e-4,
                                                                                             1.0e-4));
        final Complex integrated = integrator.integrate(100000, new Third(n, m),
                                                        new Complex(1.0e-10, 1.0e-10), phi);
        System.out.println();
        System.out.println("n:                            " + n);
        System.out.println("φ:                            " + phi);
        System.out.println("m = k²:                       " + m);
        System.out.println("Wolfram ref:                  " + ref);
        System.out.println("numerical integration:        " + integrated);
        System.out.println("Carlson-based implementation: " + Pi(n, phi, m));
        UnitTestUtils.assertEquals(ref, integrated, 2.0e-10);
        UnitTestUtils.assertEquals(ref, Pi(n, phi, m), 1.0e-10);
    }

    @Test
    public void testIssueIncompleteThirdKindB() {
        final Complex n      = new Complex(3.4, -1.3);
        final Complex phi    = new Complex(1.2, 0.0);
        final Complex m      = new Complex(2.3, -1.5);
        final Complex ref    = new Complex(0.0549975664205737508390233052870987, -0.678563304934528449620935);
        final ComplexUnivariateIntegrator integrator =
                        new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                             1.0e-4,
                                                                                             1.0e-4));
        final Complex integrated = integrator.integrate(100000, new Third(n, m),
                                                        new Complex(1.0e-10, 1.0e-10), phi);
        UnitTestUtils.assertEquals(ref, integrated, 2.0e-10);
        UnitTestUtils.assertEquals(ref, Pi(n, phi, m), 1.0e-10);
    }

    // TODO: this test fails and the numerical integral is consistent with Carlson-based implementation
    //       only Wolfram reference is different!
    @Ignore
    @Test
    public void testIssueIncompleteThirdKindC() {
        final Complex n      = new Complex(3.4, -1.3);
        final Complex phi    = new Complex(3, 2.5);
        final Complex m      = new Complex(2.3, -1.5);
        final Complex ref    = new Complex(-0.08860226061236101143265025848085778, 0.47853763883046652697121849);
        final ComplexUnivariateIntegrator integrator =
                        new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                             1.0e-4,
                                                                                             1.0e-4));
        // we have to use a specific path to get the correct result
        // integrating over a single straight line gives a completely wrong result
        final Complex integrated = integrator.integrate(100000, new Third(n, m),
                                                        new Complex(0, 0), new Complex(2.3, 0.0), phi);
        System.out.println();
        System.out.println("n:                            " + n);
        System.out.println("φ:                            " + phi);
        System.out.println("m = k²:                       " + m);
        System.out.println("Wolfram ref:                  " + ref);
        System.out.println("numerical integration:        " + integrated);
        System.out.println("Carlson-based implementation: " + Pi(n, phi, m));
        UnitTestUtils.assertEquals(ref, integrated, 5.0e-4);
        UnitTestUtils.assertEquals(ref, Pi(n, phi, m), 1.0e-10);
    }

    @Test
    public void testIssueIncompleteThirdKindD() {
        final Complex n      = new Complex(3.4, -1.3);
        final Complex phi    = new Complex(-0.4, 2.5);
        final Complex m      = new Complex(2.3, -1.5);
        final Complex ref    = new Complex(-0.088785417225639387479764237202463094, 0.47856853147720156106978019898);
        final ComplexUnivariateIntegrator integrator =
                        new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                             1.0e-4,
                                                                                             1.0e-4));
        final Complex integrated = integrator.integrate(100000, new Third(n, m),
                                                        new Complex(1.0e-10, 1.0e-10), phi);
        UnitTestUtils.assertEquals(ref, integrated, 2.0e-10);
        UnitTestUtils.assertEquals(ref, Pi(n, phi, m), 1.0e-10);
    }

    private static class Difference implements CalculusFieldUnivariateFunction<Complex> {

        final Complex m;

        Difference(final Complex m) {
            this.m = m;
        }

        public Complex value(final Complex theta) {
            final Complex sin  = theta.sin();
            final Complex sin2 = sin.multiply(sin);
            return sin2.divide(sin2.multiply(m).negate().add(1).sqrt());
        }

    }

    private static class First implements CalculusFieldUnivariateFunction<Complex> {

        final Complex m;

        First(final Complex m) {
            this.m = m;
        }

        public Complex value(final Complex theta) {
            final Complex sin  = theta.sin();
            final Complex sin2 = sin.multiply(sin);
            return sin2.multiply(m).negate().add(1).sqrt().reciprocal();
        }

    }

    private static class Second implements CalculusFieldUnivariateFunction<Complex> {

        final Complex m;

        Second(final Complex m) {
            this.m = m;
        }

        public Complex value(final Complex theta) {
            final Complex sin = theta.sin();
            final Complex sin2 = sin.multiply(sin);
            return sin2.multiply(m).negate().add(1).sqrt();
        }

    }

    private static class Third implements CalculusFieldUnivariateFunction<Complex> {

        final Complex n;
        final Complex m;

        Third(final Complex n, final Complex m) {
            this.n = n;
            this.m = m;
        }

        public Complex value(final Complex theta) {
            final Complex sin  = theta.sin();
            final Complex sin2 = sin.multiply(sin);
            final Complex d1   = sin2.multiply(m).negate().add(1).sqrt();
            final Complex da   = sin2.multiply(n).negate().add(1);
            return d1.multiply(da).reciprocal();
        }

    }

}
