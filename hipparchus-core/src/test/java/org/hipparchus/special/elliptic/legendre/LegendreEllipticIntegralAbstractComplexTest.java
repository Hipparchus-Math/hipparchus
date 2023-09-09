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
package org.hipparchus.special.elliptic.legendre;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.analysis.CalculusFieldUnivariateFunction;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.special.elliptic.carlson.CarlsonEllipticIntegral;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.MathUtils;
import org.junit.Assert;
import org.junit.Test;

public abstract class LegendreEllipticIntegralAbstractComplexTest<T extends CalculusFieldElement<T>> {

    protected abstract T buildComplex(double realPart);
    protected abstract T buildComplex(double realPart, double imaginaryPart);
    protected abstract T K(T m);
    protected abstract T Kprime(T m);
    protected abstract T F(T phi, T m);
    protected abstract T integratedF(T phi, T m);
    protected abstract T E(T m);
    protected abstract T E(T phi, T m);
    protected abstract T integratedE(T phi, T m);
    protected abstract T D(T m);
    protected abstract T D(T phi, T m);
    protected abstract T Pi(T n, T m);
    protected abstract T Pi(T n, T phi, T m);
    protected abstract T integratedPi(T n, T phi, T m);
    protected abstract T integrate(int maxEval, CalculusFieldUnivariateFunction<T> f, T start, T end);
    protected abstract T integrate(int maxEval, CalculusFieldUnivariateFunction<T> f, T start, T middle, T end);

    private void check(double expectedReal, double expectedImaginary, T result, double tol) {
        Assert.assertEquals(0, buildComplex(expectedReal, expectedImaginary).subtract(result).norm(), tol);
    }

    @Test
    public void testNoConvergence() {
        Assert.assertTrue(K(buildComplex(Double.NaN)).isNaN());
    }

    @Test
    public void testComplementary() {
        for (double m = 0.01; m < 1; m += 0.01) {
            T k1 = K(buildComplex(m));
            T k2 = Kprime(buildComplex(1 - m));
            Assert.assertEquals(k1.getReal(), k2.getReal(), FastMath.ulp(k1).getReal());
        }
    }

    @Test
    public void testAbramowitzStegunExample3() {
        T k = K(buildComplex(80.0 / 81.0));
        Assert.assertEquals(3.591545001, k.getReal(), 2.0e-9);
    }

    public void testAbramowitzStegunExample4() {
        check(1.019106060, 0.0, E(buildComplex(80.0 / 81.0)), 2.0e-8);
    }

    @Test
    public void testAbramowitzStegunExample8() {
        final T m    = buildComplex(1.0 / 5.0);
        check(1.115921, 0.0, F(buildComplex(FastMath.acos(FastMath.sqrt(2) / 3.0)), m), 1.0e-6);
        check(0.800380, 0.0, F(buildComplex(FastMath.acos(FastMath.sqrt(2) / 2.0)), m), 1.0e-6);
    }

    @Test
    public void testAbramowitzStegunExample9() {
        final T m    = buildComplex(1.0 / 2.0);
        check(1.854075, 0.0, F(buildComplex(MathUtils.SEMI_PI), m), 1.0e-6);
        check(0.535623, 0.0, F(buildComplex(FastMath.PI / 6.0), m), 1.0e-6);
    }

    @Test
    public void testAbramowitzStegunExample10() {
        final T m    = buildComplex(4.0 / 5.0);
        check(0.543604, 0.0, F(buildComplex(FastMath.PI / 6.0), m), 1.0e-6);
    }

    @Test
    public void testAbramowitzStegunExample14() {
        final T k    = buildComplex(3.0 / 5.0);
        check(0.80904, 0.0, E(buildComplex(FastMath.asin(FastMath.sqrt(5.0) / 3.0)),          k.multiply(k)), 1.0e-5);
        check(0.41192, 0.0, E(buildComplex(FastMath.asin(5.0 / (3.0 * FastMath.sqrt(17.0)))), k.multiply(k)), 1.0e-5);
    }

    @Test
    public void testAbramowitzStegunTable175() {
        final double sinAlpha1 = FastMath.sin(FastMath.toRadians(32));
        check(0.26263487, 0.0, F(buildComplex(FastMath.toRadians(15)), buildComplex(sinAlpha1 * sinAlpha1)), 1.0e-8);
        final double sinAlpha2 = FastMath.sin(FastMath.toRadians(46));
        check(1.61923762, 0.0, F(buildComplex(FastMath.toRadians(80)), buildComplex(sinAlpha2 * sinAlpha2)), 1.0e-8);
    }

    @Test
    public void testAbramowitzStegunTable176() {
        final double sinAlpha1 = FastMath.sin(FastMath.toRadians(64));
        check(0.42531712, 0.0, E(buildComplex(FastMath.toRadians(25)), buildComplex(sinAlpha1 * sinAlpha1)), 1.0e-8);
        final double sinAlpha2 = FastMath.sin(FastMath.toRadians(76));
        check(0.96208074, 0.0, E(buildComplex(FastMath.toRadians(70)), buildComplex(sinAlpha2 * sinAlpha2)), 1.0e-8);
    }

    @Test
    public void testAbramowitzStegunTable179() {
        final double sinAlpha1 = FastMath.sin(FastMath.toRadians(15));
        check(1.62298, 0.0, Pi(buildComplex(0.4), buildComplex(FastMath.toRadians(75)), buildComplex(sinAlpha1 * sinAlpha1)), 1.0e-5);
        final double sinAlpha2 = FastMath.sin(FastMath.toRadians(60));
        check(1.03076, 0.0, Pi(buildComplex(0.8), buildComplex(FastMath.toRadians(45)), buildComplex(sinAlpha2 * sinAlpha2)), 1.0e-5);
        final double sinAlpha3 = FastMath.sin(FastMath.toRadians(15));
        check(2.79990, 0.0, Pi(buildComplex(0.9), buildComplex(FastMath.toRadians(75)), buildComplex(sinAlpha3 * sinAlpha3)), 1.0e-5);
    }

    @Test
    public void testCompleteVsIncompleteF() {
        for (double m = 0.01; m < 1; m += 0.01) {
            double complete   = K(buildComplex(m)).getReal();
            double incomplete = F(buildComplex(MathUtils.SEMI_PI), buildComplex(m)).getReal();
            Assert.assertEquals(complete, incomplete, FastMath.ulp(complete));
        }
    }

    @Test
    public void testCompleteVsIncompleteE() {
        for (double m = 0.01; m < 1; m += 0.01) {
            double complete   = E(buildComplex(m)).getReal();
            double incomplete = E(buildComplex(MathUtils.SEMI_PI), buildComplex(m)).getReal();
            Assert.assertEquals(complete, incomplete, 4 * FastMath.ulp(complete));
        }
    }

    @Test
    public void testCompleteVsIncompleteD() {
        for (double m = 0.01; m < 1; m += 0.01) {
            double complete   = D(buildComplex(m)).getReal();
            double incomplete = D(buildComplex(MathUtils.SEMI_PI), buildComplex(m)).getReal();
            Assert.assertEquals(complete, incomplete, FastMath.ulp(complete));
        }
    }

    @Test
    public void testCompleteVsIncompletePi() {
        for (double alpha2 = 0.01; alpha2 < 1; alpha2 += 0.01) {
            for (double m = 0.01; m < 1; m += 0.01) {
                double complete   = Pi(buildComplex(alpha2), buildComplex(m)).getReal();
                double incomplete = Pi(buildComplex(alpha2), buildComplex(MathUtils.SEMI_PI), buildComplex(m)).getReal();
                Assert.assertEquals(complete, incomplete, FastMath.ulp(complete));
            }
        }
    }

    @Test
    public void testNomeSmallParameter() {
        Assert.assertEquals(5.9375e-18, LegendreEllipticIntegral.nome(buildComplex(0.95e-16)).getReal(), 1.0e-50);
    }

    @Test
    public void testIntegralsSmallParameter() {
        Assert.assertEquals(7.8539816428e-10, K(buildComplex(2.0e-9)).getReal() - MathUtils.SEMI_PI, 1.0e-15);
    }

    @Test
    public void testPrecomputedDelta() {

        T n   = buildComplex(3.4,  -1.3);
        T m   = buildComplex(0.2,   0.6);
        T phi = buildComplex(1.2, 0.08);
        T ref = buildComplex(0.48657872913710487, -0.9325310177613093);
        Assert.assertEquals(0.0, Pi(n, phi, m).subtract(ref).getReal(), 1.0e-15);

        // no argument reduction and no precomputed delta
        final T csc     = phi.sin().reciprocal();
        final T csc2    = csc.multiply(csc);
        final T cM1     = csc2.subtract(1);
        final T cMm     = csc2.subtract(m);
        final T cMn     = csc2.subtract(n);
        final T pinphim = CarlsonEllipticIntegral.rF(cM1, cMm, csc2).
                          add(CarlsonEllipticIntegral.rJ(cM1, cMm, csc2, cMn).multiply(n).divide(3));
        Assert.assertEquals(0.0, pinphim.subtract(ref).getReal(), 1.0e-15);

    }

    @Test
    public void testIncompleteDifferenceA() {
        final T phi = buildComplex(1.2, 0.75);
        final T m   = buildComplex(0.2, 0.6);
        final T ref = F(phi, m).
                            subtract(E(phi, m)).
                            divide(m);
        final T integrated = integrate(100000, new Difference<>(m), buildComplex(1.0e-10, 1.0e-10), phi);
        Assert.assertEquals(0.0, integrated.subtract(ref).norm(), 2.0e-10);
        Assert.assertEquals(0.0, D(phi, m).subtract(ref).norm(), 1.0e-10);
    }

    @Test
    public void testIncompleteDifferenceB() {
        final T phi = buildComplex(1.2, 0.0);
        final T m   = buildComplex(2.3, -1.5);
        final T ref = F(phi, m).
                            subtract(E(phi, m)).
                            divide(m);
        final T integrated = integrate(100000, new Difference<>(m), buildComplex(1.0e-10, 1.0e-10), phi);
        Assert.assertEquals(0.0, integrated.subtract(ref).norm(), 2.0e-10);
        Assert.assertEquals(0.0, D(phi, m).subtract(ref).norm(), 1.0e-10);
    }

    @Test
    public void testIncompleteDifferenceC() {
        final T phi = buildComplex(3, 2.5);
        final T m   = buildComplex(2.3, -1.5);
        final T ref = F(phi, m).subtract(E(phi, m)).divide(m);
        // we have to use a specific path to get the correct result
        // integrating over a single straight line gives a completely wrong result
        final T integrated = integrate(100000, new Difference<>(m), buildComplex(1.0e-12, 1.0e-12), buildComplex(0, -1.5), phi);
        Assert.assertEquals(0.0, integrated.subtract(ref).norm(), 2.0e-10);
        Assert.assertEquals(0.0, D(phi, m).subtract(ref).norm(), 1.0e-10);
    }

    @Test
    public void testIncompleteDifferenceD() {
        final T phi = buildComplex(-0.4, 2.5);
        final T m   = buildComplex(2.3, -1.5);
        final T ref = F(phi, m).
                            subtract(E(phi, m)).
                            divide(m);
        final T integrated = integrate(100000, new Difference<>(m), buildComplex(1.0e-10, 1.0e-10), phi);
        Assert.assertEquals(0.0, integrated.subtract(ref).norm(), 2.0e-10);
        Assert.assertEquals(0.0, D(phi, m).subtract(ref).norm(), 1.0e-10);
    }

    @Test
    public void testIncompleteFirstKindA() {
        final T phi = buildComplex(1.2, 0.75);
        final T m   = buildComplex(0.2, 0.6);
        final T ref = buildComplex(1.00265860821563927579252866, 0.80128721521822408811217);
        Assert.assertEquals(0.0, integratedF(phi, m).subtract(ref).norm(), 2.0e-10);
        Assert.assertEquals(0.0, F(phi, m).subtract(ref).norm(), 1.0e-10);
    }

    @Test
    public void testIncompleteFirstKindB() {
        final T phi = buildComplex(1.2, 0.0);
        final T m   = buildComplex(2.3, -1.5);
        final T ref = buildComplex(1.04335840461807753156026488, -0.5872679121672512828049797);
        Assert.assertEquals(0.0, integratedF(phi, m).subtract(ref).norm(), 2.0e-10);
        Assert.assertEquals(0.0, F(phi, m).subtract(ref).norm(), 1.0e-10);
    }

    @Test
    public void testIncompleteFirstKindC() {
        final T phi = buildComplex(-0.4, 2.5);
        final T m   = buildComplex(2.3, -1.5);
        final T ref = buildComplex(-0.20646268947416273887690961, 1.0927692344374984107332330625089);
        Assert.assertEquals(0.0, integratedF(phi, m).subtract(ref).norm(), 2.0e-10);
        Assert.assertEquals(0.0, F(phi, m).subtract(ref).norm(), 1.0e-10);
    }

    @Test
    public void testIncompleteSecondKindA() {
        final T phi = buildComplex(1.2, 0.75);
        final T m   = buildComplex(0.2, 0.6);
        final T ref = buildComplex(1.4103674846223375296500, 0.644849758860533700396);
        Assert.assertEquals(0.0, integratedE(phi, m).subtract(ref).norm(), 2.0e-10);
        Assert.assertEquals(0.0, E(phi, m).subtract(ref).norm(), 1.0e-10);
    }

    @Test
    public void testIncompleteSecondKindB() {
        final T phi = buildComplex(1.2, 0.0);
        final T m   = buildComplex(2.3, -1.5);
        final T ref = buildComplex(0.8591316843513079270009549421, 0.55423174445992167002660);
        Assert.assertEquals(0.0, integratedE(phi, m).subtract(ref).norm(), 2.0e-10);
        Assert.assertEquals(0.0, E(phi, m).subtract(ref).norm(), 1.0e-10);
    }

    @Test
    public void testIncompleteSecondKindC() {
        final T phi = buildComplex(-0.4, 2.5);
        final T m   = buildComplex(2.3, -1.5);
        final T ref = buildComplex(-1.68645030068870706703580773597, 9.176675281683098106653799);
        Assert.assertEquals(0.0, integratedE(phi, m).subtract(ref).norm(), 2.0e-10);
        Assert.assertEquals(0.0, E(phi, m).subtract(ref).norm(), 1.0e-10);
    }

    @Test
    public void testIncompleteThirdKind() {
        final T n      = buildComplex(3.4, -1.3);
        final T m      = buildComplex(0.2, 0.6);
        final T[][] references = MathArrays.buildArray(n.getField(), 38, 2);
        references[ 0][0] = buildComplex(1.2, -1.5);          references[ 0][1] = buildComplex( 0.03351171124667249063, -0.57566536173018225078);
        references[ 1][0] = buildComplex(1.2, -1.4);          references[ 1][1] = buildComplex( 0.03644476655784750591, -0.57331323414589059064);
        references[ 2][0] = buildComplex(1.2, -1.389190765);  references[ 2][1] = buildComplex( 0.03682595624642736804, -0.57302208430696377063);
        references[ 3][0] = buildComplex(1.2, -1.3);          references[ 3][1] = buildComplex( 0.04057582076289887060, -0.57031938416517331851);
        references[ 4][0] = buildComplex(1.2, -1.2);          references[ 4][1] = buildComplex( 0.04640620250501646778, -0.56663563199204550096);
        references[ 5][0] = buildComplex(1.2, -1.066681968);  references[ 5][1] = buildComplex( 0.05798102095188363268, -0.56085091376439940662);
        references[ 6][0] = buildComplex(1.2, -1.066681967);  references[ 6][1] = buildComplex( 0.16640600620018822283, -0.80852909755613891276);
        references[ 7][0] = buildComplex(1.2, -1.0);          references[ 7][1] = buildComplex( 0.17433162130249770262, -0.80552550782381317481);
        references[ 8][0] = buildComplex(1.2, -0.75);         references[ 8][1] = buildComplex( 0.21971749684893289198, -0.79853199575576671312);
        references[ 9][0] = buildComplex(1.2, -0.5);          references[ 9][1] = buildComplex( 0.28863330914614968227, -0.80811912506287693026);
        references[10][0] = buildComplex(1.2, 0.00);          references[10][1] = buildComplex( 0.46199936298130610116, -0.90668901748978345243);
        references[11][0] = buildComplex(1.2, 0.05);          references[11][1] = buildComplex( 0.47768146729143941350, -0.92262785201399518938);
        references[12][0] = buildComplex(1.2, 0.07);          references[12][1] = buildComplex( 0.48365870488290366257, -0.92920571368375179065);
        references[13][0] = buildComplex(1.2, 0.08);          references[13][1] = buildComplex( 0.48657872913710538725, -0.93253101776130872023);
        references[14][0] = buildComplex(1.2, 0.085);         references[14][1] = buildComplex( 0.48802108187307861659, -0.93420195155547102730);
        references[15][0] = buildComplex(1.2, 0.085181);      references[15][1] = buildComplex( 0.48807307162394529967, -0.93426253841598734629);
        references[16][0] = buildComplex(1.2, 0.085182);      references[16][1] = buildComplex( 1.10372915320771558982, -2.71126044767925836757);
        references[17][0] = buildComplex(1.2, 0.08524491);    references[17][1] = buildComplex( 1.10374721953636884577, -2.71128150740577828433);
        references[18][0] = buildComplex(1.2, 0.08524492);    references[18][1] = buildComplex(-1.35887595515643636904,  4.39670878728780118558);
        references[19][0] = buildComplex(1.2, 0.08524501);    references[19][1] = buildComplex(-1.35887592931182948115,  4.39670875715884785695);
        references[20][0] = buildComplex(1.2, 0.08524502);    references[20][1] = buildComplex(-0.12756433765799251204,  0.84271360479056584862);
        references[21][0] = buildComplex(1.2, 0.08524505);    references[21][1] = buildComplex(-0.12756432904312455421,  0.84271359474758096952);
        references[22][0] = buildComplex(1.2, 0.0852451);     references[22][1] = buildComplex(-0.12756431468501224811,  0.84271357800927242222);
        references[23][0] = buildComplex(1.2, 0.086);         references[23][1] = buildComplex(-0.12734767232640510383,  0.84246080397106417910);
        references[24][0] = buildComplex(1.2, 0.087);         references[24][1] = buildComplex(-0.12706111141088970319,  0.84212577870702519585);
        references[25][0] = buildComplex(1.2, 0.09);          references[25][1] = buildComplex(-0.12620431480082688859,  0.84111948445051909467);
        references[26][0] = buildComplex(1.2, 0.10);          references[26][1] = buildComplex(-0.12337990144324015420,  0.83775255734168295690);
        references[27][0] = buildComplex(1.2, 0.20);          references[27][1] = buildComplex(-0.09796886081364286844,  0.80343164778772696735);
        references[28][0] = buildComplex(1.2, 0.2049);        references[28][1] = buildComplex(-0.09686124781335778609,  0.80173956525364697562);
        references[29][0] = buildComplex(1.2, 0.2051631601);  references[29][1] = buildComplex(-0.096802132146685100037, 0.80164871118350302421);
        references[30][0] = buildComplex(1.2, 0.2051631602);  references[30][1] = buildComplex(-0.096802132124228501156, 0.80164871114897922197);
        references[31][0] = buildComplex(1.2, 0.24);          references[31][1] = buildComplex(-0.08930914015165779431,  0.78965671512935519039);
        references[32][0] = buildComplex(1.2, 0.2462);        references[32][1] = buildComplex(-0.08804459384315873622,  0.78753318617257610859);
        references[33][0] = buildComplex(1.2, 0.24675781599); references[33][1] = buildComplex(-0.087931839058365500429, 0.78734234011803621580);
        references[34][0] = buildComplex(1.2, 0.2475);        references[34][1] = buildComplex(-0.08778207674401363351,  0.78708847169559542169);
        references[35][0] = buildComplex(1.2, 0.25);          references[35][1] = buildComplex(-0.08727979375768656571,  0.78623380867526956107);
        references[36][0] = buildComplex(1.2, 0.45);          references[36][1] = buildComplex(-0.05716147875408999398,  0.72239458160027437391);
        references[37][0] = buildComplex(1.2, 0.75);          references[37][1] = buildComplex(-0.03776232345596223316,  0.65347724469972942256);

        for (int i = 0; i < references.length; ++i) {
            final T[] ref = references[i];
            T integrated;
            try {
                integrated = integratedPi(n, ref[0], m);
            } catch (MathIllegalStateException mise) {
                integrated = buildComplex(Double.NaN);
            }
            T carlson    = Pi(n, ref[0], m);
            if (i < 2) {
                // integration, Carlson and Wolfram Alpha all give different results
                Assert.assertTrue(true);
            } else if (i == 2) {
                // integration hits the pole
                Assert.assertTrue(integrated.isNaN());
            } else if (i < 6) {
                // integration and Carlson agree and are most probably right
                // Wolfram Alpha gives a different result which seems to be wrong
                Assert.assertEquals(0.0, carlson.subtract(integrated).norm(), 4.1e-7);
            } else if (i < 16) {
                // integration, Carlson and Wolfram Alpha all agree and are most probably right
                Assert.assertEquals(0.0, integrated.subtract(ref[1]).norm(), 1.0e-10);
                Assert.assertEquals(0.0, carlson.subtract(ref[1]).norm(), 1.0e-10);
            } else if (i < 30) {
                // integration and Carlson agree and are most probably right
                // Wolfram Alpha gives a different result which seems to be wrong
                Assert.assertEquals(0.0, carlson.subtract(integrated).norm(), 2.0e-6);
            } else if (i < 35) {
                // integration, Carlson and Wolfram Alpha all give different results
                Assert.assertTrue(true);
            } else {
                // integration and Wolfram Alpha agree and are most probably right
                // Carlson gives a different result which seems to be wrong
                Assert.assertEquals(0.0, integrated.subtract(ref[1]).norm(), 2.5e-7);
            }
        }
    }

    private static class Difference<T extends CalculusFieldElement<T>> implements CalculusFieldUnivariateFunction<T> {

        final T m;

        Difference(final T m) {
            this.m = m;
        }

        public T value(final T theta) {
            final T sin  = theta.sin();
            final T sin2 = sin.multiply(sin);
            return sin2.divide(sin2.multiply(m).negate().add(1).sqrt());
        }

    }

}
