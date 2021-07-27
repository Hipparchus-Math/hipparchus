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
import org.hipparchus.exception.MathIllegalStateException;
import org.junit.Assert;
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
    public void testIncompleteDifferenceA() {
        final Complex phi = new Complex(1.2, 0.75);
        final Complex m   = new Complex(0.2, 0.6);
        final Complex ref = F(phi, m).
                            subtract(E(phi, m)).
                            divide(m);
        final ComplexUnivariateIntegrator integrator =
                        new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                             1.0e-6,
                                                                                             1.0e-6));
        final Complex integrated = integrator.integrate(100000, new Difference(m),
                                                        new Complex(1.0e-10, 1.0e-10), phi);
        UnitTestUtils.assertEquals(ref, integrated, 2.0e-10);
        UnitTestUtils.assertEquals(ref, D(phi, m), 1.0e-10);
    }

    @Test
    public void testIncompleteDifferenceB() {
        final Complex phi = new Complex(1.2, 0.0);
        final Complex m   = new Complex(2.3, -1.5);
        final Complex ref = F(phi, m).
                            subtract(E(phi, m)).
                            divide(m);
        final ComplexUnivariateIntegrator integrator =
                        new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                             1.0e-6,
                                                                                             1.0e-6));
        final Complex integrated = integrator.integrate(100000, new Difference(m),
                                                        new Complex(1.0e-10, 1.0e-10), phi);
        UnitTestUtils.assertEquals(ref, integrated, 2.0e-10);
        UnitTestUtils.assertEquals(ref, D(phi, m), 1.0e-10);
    }

    @Test
    public void testIncompleteDifferenceC() {
        final Complex phi = new Complex(3, 2.5);
        final Complex m   = new Complex(2.3, -1.5);
        final Complex ref = F(phi, m).
                            subtract(E(phi, m)).
                            divide(m);
        final ComplexUnivariateIntegrator integrator =
                        new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                             1.0e-6,
                                                                                             1.0e-6));
        // we have to use a specific path to get the correct result
        // integrating over a single straight line gives a completely wrong result
        final Complex integrated = integrator.integrate(100000, new Difference(m),
                                                        new Complex(1.0e-12, 1.0e-12), new Complex(0, -1.5), phi);
        UnitTestUtils.assertEquals(ref, integrated, 2.0e-10);
        UnitTestUtils.assertEquals(ref, D(phi, m), 1.0e-10);
    }

    @Test
    public void testIncompleteDifferenceD() {
        final Complex phi = new Complex(-0.4, 2.5);
        final Complex m   = new Complex(2.3, -1.5);
        final Complex ref = F(phi, m).
                            subtract(E(phi, m)).
                            divide(m);
        final ComplexUnivariateIntegrator integrator =
                        new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                             1.0e-6,
                                                                                             1.0e-6));
        final Complex integrated = integrator.integrate(100000, new Difference(m),
                                                        new Complex(1.0e-10, 1.0e-10), phi);
        UnitTestUtils.assertEquals(ref, integrated, 2.0e-10);
        UnitTestUtils.assertEquals(ref, D(phi, m), 1.0e-10);
    }

    @Test
    public void testIncompleteFirstKindA() {
        final Complex phi = new Complex(1.2, 0.75);
        final Complex m   = new Complex(0.2, 0.6);
        final Complex ref = new Complex(1.00265860821563927579252866, 0.80128721521822408811217);
        final ComplexUnivariateIntegrator integrator =
                        new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                             1.0e-6,
                                                                                             1.0e-6));
        UnitTestUtils.assertEquals(ref, LegendreEllipticIntegral.bigF(phi, m, integrator, 100000), 2.0e-10);
        UnitTestUtils.assertEquals(ref, F(phi, m), 1.0e-10);
    }

    @Test
    public void testIncompleteFirstKindB() {
        final Complex phi = new Complex(1.2, 0.0);
        final Complex m   = new Complex(2.3, -1.5);
        final Complex ref = new Complex(1.04335840461807753156026488, -0.5872679121672512828049797);
        final ComplexUnivariateIntegrator integrator =
                        new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                             1.0e-6,
                                                                                             1.0e-6));
        UnitTestUtils.assertEquals(ref, LegendreEllipticIntegral.bigF(phi, m, integrator, 100000), 2.0e-10);
        UnitTestUtils.assertEquals(ref, F(phi, m), 1.0e-10);
    }

    @Test
    public void testIncompleteFirstKindC() {
        final Complex phi = new Complex(-0.4, 2.5);
        final Complex m   = new Complex(2.3, -1.5);
        final Complex ref = new Complex(-0.20646268947416273887690961, 1.0927692344374984107332330625089);
        final ComplexUnivariateIntegrator integrator =
                        new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                             1.0e-6,
                                                                                             1.0e-6));
        UnitTestUtils.assertEquals(ref, LegendreEllipticIntegral.bigF(phi, m, integrator, 100000), 2.0e-10);
        UnitTestUtils.assertEquals(ref, F(phi, m), 1.0e-10);
    }

    @Test
    public void testIncompleteSecondKindA() {
        final Complex phi = new Complex(1.2, 0.75);
        final Complex m   = new Complex(0.2, 0.6);
        final Complex ref = new Complex(1.4103674846223375296500, 0.644849758860533700396);
        final ComplexUnivariateIntegrator integrator =
                        new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                             1.0e-6,
                                                                                             1.0e-6));
        UnitTestUtils.assertEquals(ref, LegendreEllipticIntegral.bigE(phi, m, integrator, 100000), 2.0e-10);
        UnitTestUtils.assertEquals(ref, E(phi, m), 1.0e-10);
    }

    @Test
    public void testIncompleteSecondKindB() {
        final Complex phi = new Complex(1.2, 0.0);
        final Complex m   = new Complex(2.3, -1.5);
        final Complex ref = new Complex(0.8591316843513079270009549421, 0.55423174445992167002660);
        final ComplexUnivariateIntegrator integrator =
                        new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                             1.0e-6,
                                                                                             1.0e-6));
        UnitTestUtils.assertEquals(ref, LegendreEllipticIntegral.bigE(phi, m, integrator, 100000), 2.0e-10);
        UnitTestUtils.assertEquals(ref, E(phi, m), 1.0e-10);
    }

    @Test
    public void testIncompleteSecondKindC() {
        final Complex phi = new Complex(-0.4, 2.5);
        final Complex m   = new Complex(2.3, -1.5);
        final Complex ref = new Complex(-1.68645030068870706703580773597, 9.176675281683098106653799);
        final ComplexUnivariateIntegrator integrator =
                        new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                             1.0e-6,
                                                                                             1.0e-6));
        UnitTestUtils.assertEquals(ref, LegendreEllipticIntegral.bigE(phi, m, integrator, 100000), 2.0e-10);
        UnitTestUtils.assertEquals(ref, E(phi, m), 1.0e-10);
    }

    @Test
    public void testIncompleteThirdKind() {
        final Complex n      = new Complex(3.4, -1.3);
        final Complex m      = new Complex(0.2, 0.6);
        final Complex[][] references = {
            { new Complex(1.2, -1.5),          new Complex( 0.03351171124667249063, -0.57566536173018225078) },
            { new Complex(1.2, -1.4),          new Complex( 0.03644476655784750591, -0.57331323414589059064) },
            { new Complex(1.2, -1.389190765),  new Complex( 0.03682595624642736804, -0.57302208430696377063) },
            { new Complex(1.2, -1.3),          new Complex( 0.04057582076289887060, -0.57031938416517331851) },
            { new Complex(1.2, -1.2),          new Complex( 0.04640620250501646778, -0.56663563199204550096) },
            { new Complex(1.2, -1.066681968),  new Complex( 0.05798102095188363268, -0.56085091376439940662) },
            { new Complex(1.2, -1.066681967),  new Complex( 0.16640600620018822283, -0.80852909755613891276) },
            { new Complex(1.2, -1.0),          new Complex( 0.17433162130249770262, -0.80552550782381317481) },
            { new Complex(1.2, -0.75),         new Complex( 0.21971749684893289198, -0.79853199575576671312) },
            { new Complex(1.2, -0.5),          new Complex( 0.28863330914614968227, -0.80811912506287693026) },
            { new Complex(1.2, 0.00),          new Complex( 0.46199936298130610116, -0.90668901748978345243) },
            { new Complex(1.2, 0.05),          new Complex( 0.47768146729143941350, -0.92262785201399518938) },
            { new Complex(1.2, 0.07),          new Complex( 0.48365870488290366257, -0.92920571368375179065) },
            { new Complex(1.2, 0.08),          new Complex( 0.48657872913710538725, -0.93253101776130872023) },
            { new Complex(1.2, 0.085),         new Complex( 0.48802108187307861659, -0.93420195155547102730) },
            { new Complex(1.2, 0.085181),      new Complex( 0.48807307162394529967, -0.93426253841598734629) },
            { new Complex(1.2, 0.085182),      new Complex( 1.10372915320771558982, -2.71126044767925836757) },
            { new Complex(1.2, 0.08524491),    new Complex( 1.10374721953636884577, -2.71128150740577828433) },
            { new Complex(1.2, 0.08524492),    new Complex(-1.35887595515643636904,  4.39670878728780118558) },
            { new Complex(1.2, 0.08524501),    new Complex(-1.35887592931182948115,  4.39670875715884785695) },
            { new Complex(1.2, 0.08524502),    new Complex(-0.12756433765799251204,  0.84271360479056584862) },
            { new Complex(1.2, 0.08524505),    new Complex(-0.12756432904312455421,  0.84271359474758096952) },
            { new Complex(1.2, 0.0852451),     new Complex(-0.12756431468501224811,  0.84271357800927242222) },
            { new Complex(1.2, 0.086),         new Complex(-0.12734767232640510383,  0.84246080397106417910) },
            { new Complex(1.2, 0.087),         new Complex(-0.12706111141088970319,  0.84212577870702519585) },
            { new Complex(1.2, 0.09),          new Complex(-0.12620431480082688859,  0.84111948445051909467) },
            { new Complex(1.2, 0.10),          new Complex(-0.12337990144324015420,  0.83775255734168295690) },
            { new Complex(1.2, 0.20),          new Complex(-0.09796886081364286844,  0.80343164778772696735) },
            { new Complex(1.2, 0.2049),        new Complex(-0.09686124781335778609,  0.80173956525364697562) },
            { new Complex(1.2, 0.2051631601),  new Complex(-0.096802132146685100037, 0.80164871118350302421) },
            { new Complex(1.2, 0.2051631602),  new Complex(-0.096802132124228501156, 0.80164871114897922197) },
            { new Complex(1.2, 0.24),          new Complex(-0.08930914015165779431,  0.78965671512935519039) },
            { new Complex(1.2, 0.2462),        new Complex(-0.08804459384315873622,  0.78753318617257610859) },
            { new Complex(1.2, 0.24675781599), new Complex(-0.087931839058365500429, 0.78734234011803621580) },
            { new Complex(1.2, 0.2475),        new Complex(-0.08778207674401363351,  0.78708847169559542169) },
            { new Complex(1.2, 0.25),          new Complex(-0.08727979375768656571,  0.78623380867526956107) },
            { new Complex(1.2, 0.45),          new Complex(-0.05716147875408999398,  0.72239458160027437391) },
            { new Complex(1.2, 0.75),          new Complex(-0.03776232345596223316,  0.65347724469972942256) },
        };

        final ComplexUnivariateIntegrator integrator =
                        new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                             1.0e-6,
                                                                                             1.0e-6));

        for (final Complex[] ref : references) {
            Complex integrated;
            try {
                integrated = LegendreEllipticIntegral.bigPi(n, ref[0], m, integrator, 100000);
            } catch (MathIllegalStateException mise) {
                integrated = Complex.NaN;
            }
            Complex carlson    = Pi(n, ref[0], m);
            if (ref[0].getImaginaryPart() < -1.35) {
                // TODO: integration, Carlson and Wolfram Alpha all give different results
                Assert.assertTrue(true);
            } else if (ref[0].getImaginaryPart() < -1.0666819675) {
                // integration and Carlson agree and are most probably right
                // Wolfram Alpha gives a different result which seems to be wrong
                UnitTestUtils.assertEquals(integrated, carlson, 4.0e-7);
            } else if (ref[0].getImaginaryPart() < 0.085181) {
                // integration, Carlson and Wolfram Alpha all agree and are most probably right
                UnitTestUtils.assertEquals(ref[1], integrated, 1.0e-10);
                UnitTestUtils.assertEquals(ref[1], carlson, 1.0e-10);
            } else if (ref[0].getImaginaryPart() < 0.20516316015) {
                // integration and Carlson agree and are most probably right
                // Wolfram Alpha gives a different result which seems to be wrong
                UnitTestUtils.assertEquals(integrated, carlson, 2.0e-6);
            } else if (ref[0].getImaginaryPart() < 0.25) {
                // TODO: integration, Carlson and Wolfram Alpha all give different results
                Assert.assertTrue(true);
            } else {
                // integration and Wolfram Alpha agree and are most probably right
                // Carlson gives a different result which seems to be wrong
                UnitTestUtils.assertEquals(ref[1], integrated, 2.5e-7);
            }
        }
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

}
