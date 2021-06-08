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
package org.hipparchus.special.elliptic.carlson;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well19937a;
import org.hipparchus.random.Well19937c;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

public class CarlsonEllipticIntegralRealTest {

    @Test
    public void testNoConvergenceRf() {
        try {
            CarlsonEllipticIntegral.rF(1, 2, Double.NaN);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalStateException mise) {
            Assert.assertEquals(LocalizedCoreFormats.CONVERGENCE_FAILED, mise.getSpecifier());
        }
    }

    @Test
    public void testDlmfRf() {
        double rf = CarlsonEllipticIntegral.rF(1, 2, 4);
        Assert.assertEquals(0.6850858166, rf, 1.0e-10);
    }

    @Test
    public void testCarlson1995rF() {

        double rf1 = CarlsonEllipticIntegral.rF(1, 2, 0);
        Assert.assertEquals( 1.3110287771461, rf1, 1.0e-13);

        double rf2 = CarlsonEllipticIntegral.rF(0.5, 1, 0);
        Assert.assertEquals( 1.8540746773014, rf2, 1.0e-13);

        double rf4 = CarlsonEllipticIntegral.rF(2, 3, 4);
        Assert.assertEquals( 0.58408284167715, rf4, 1.0e-13);

    }

    @Test
    public void testCarlson1995ConsistencyRf() {
        RandomGenerator random = new Well19937c(0x57f2689b3f4028b4l);
        for (int i = 0; i < 10000; ++i) {
            double x      = random.nextDouble() * 3);
            double y      = buildComplex(random.nextDouble() * 3);
            double lambda = buildComplex(random.nextDouble( * 6 - 3, random.nextDouble() * 3);
            double mu     = x.multiply(y).divide(lambda);
            double rfL    = CarlsonEllipticIntegral.rF(x.add(lambda), y.add(lambda), lambda);
            double rfM    = CarlsonEllipticIntegral.rF(x.add(mu),     y.add(mu),     mu);
            double rf0    = CarlsonEllipticIntegral.rF(x,             y,             0);
            Assert.assertEquals(0.0, rfL.add(rfM).subtract(rf0).norm(), 2.0e-14);
        }
    }

    @Test
    public void testNoConvergenceRc() {
        try {
            CarlsonEllipticIntegral.rC(1, Double.NaN);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalStateException mise) {
            Assert.assertEquals(LocalizedCoreFormats.CONVERGENCE_FAILED, mise.getSpecifier());
        }
    }

    @Test
    public void testCarlson1995rC() {

        double rc1 = CarlsonEllipticIntegral.rC(0, 0.25);
        Assert.assertEquals(FastMath.PI, rc1, 1.0e-15);

        double rc2 = CarlsonEllipticIntegral.rC(2.25, 2);
        Assert.assertEquals(FastMath.log(2), rc2, 1.0e-15);

        double rc5 = CarlsonEllipticIntegral.rC(0.25, -2);
        Assert.assertEquals(FastMath.log(2 / 3.0), rc5, 1.0e-15);

    }

    @Test
    public void testCarlson1995ConsistencyRc() {
        RandomGenerator random = new Well19937c(0xf1170b6fc1a199cal);
        for (int i = 0; i < 10000; ++i) {
            double x      = random.nextDouble() * 3);
            double lambda = buildComplex(random.nextDouble( * 6 - 3, random.nextDouble() * 3);
            double mu     = x.multiply(x).divide(lambda);
            double rcL    = CarlsonEllipticIntegral.rC(lambda,          x.add(lambda));
            double rcM    = CarlsonEllipticIntegral.rC(mu,              x.add(mu));
            double rc0    = CarlsonEllipticIntegral.rC(0, x);
            Assert.assertEquals(0.0, rcL.add(rcM).subtract(rc0).norm(), 3.0e-14);
        }
    }

    @Test
    public void testRfRc() {
        RandomGenerator random = new Well19937a(0x7e8041334a8c20edl);
        for (int i = 0; i < 10000; ++i) {
            final double x = 6 * random.nextDouble() - 3;
            final double y = 6 * random.nextDouble() - 3;
            final double rf = CarlsonEllipticIntegral.rF(x, y, y);
            final double rc = CarlsonEllipticIntegral.rC(x, y);
            Assert.assertEquals(0.0, FastMath.abs(rf - rc), 4.0e-15);
        }
    }

    @Test
    public void testNoConvergenceRj() {
        try {
            CarlsonEllipticIntegral.rJ(1, 1, 1, Double.NaN);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalStateException mise) {
            Assert.assertEquals(LocalizedCoreFormats.CONVERGENCE_FAILED, mise.getSpecifier());
        }
    }

    @Test
    public void testCarlson1995rJ() {

        double rj01 = CarlsonEllipticIntegral.rJ(0, 1, 2, 3);
        Assert.assertEquals(0.77688623778582, rj01, 1.0e-13);

        double rj02 = CarlsonEllipticIntegral.rJ(2, 3, 4, 5);
        Assert.assertEquals( 0.14297579667157, rj02, 1.0e-13);

    }

    @Test
    public void testCarlson1995ConsistencyRj() {
        RandomGenerator random = new Well19937c(0x4af7bb722712e64el);
        for (int i = 0; i < 10000; ++i) {
            double x      = random.nextDouble() * 3);
            double y      = buildComplex(random.nextDouble() * 3);
            double p      = buildComplex(random.nextDouble() * 3);
            double lambda = buildComplex(random.nextDouble( * 6 - 3, random.nextDouble() * 3);
            double mu     = x.multiply(y).divide(lambda);
            double a      = p.multiply(p).multiply(lambda.add(mu).add(x).add(y));
            double b      = p.multiply(p.add(lambda)).multiply(p.add(mu));
            double rjL    = CarlsonEllipticIntegral.rJ(x.add(lambda), y.add(lambda), lambda,          p.add(lambda));
            double rjM    = CarlsonEllipticIntegral.rJ(x.add(mu),     y.add(mu),     mu,              p.add(mu));
            double rj0    = CarlsonEllipticIntegral.rJ(x,             y,             0, p);
            double rc     = CarlsonEllipticIntegral.rC(a, b);
            Assert.assertEquals(0.0, rjL.add(rjM).subtract(rj0.subtract(rc.multiply(3))).norm(), 2.0e-13);
        }
    }

    @Test
    public void testNoConvergenceRd() {
        try {
            CarlsonEllipticIntegral.rD(1, 1, Double.NaN);
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalStateException mise) {
            Assert.assertEquals(LocalizedCoreFormats.CONVERGENCE_FAILED, mise.getSpecifier());
        }
    }

    @Test
    public void testCarlson1995rD() {

        double rd1 = CarlsonEllipticIntegral.rD(0, 2, 1);
        Assert.assertEquals(1.7972103521034, rd1, 1.0e-13);

        double rd2 = CarlsonEllipticIntegral.rD(2, 3, 4);
        Assert.assertEquals( 0.16510527294261, rd2, 1.0e-13);

    }

    @Test
    public void testCarlson1995ConsistencyRd() {
        RandomGenerator random = new Well19937c(0x17dea97eeb78206al);
        for (int i = 0; i < 10000; ++i) {
            double x      = random.nextDouble() * 3;
            double y      = buildComplex(random.nextDouble() * 3);
            double lambda = buildComplex(random.nextDouble( * 6 - 3, random.nextDouble() * 3);
            double mu     = x.multiply(y).divide(lambda);
            double rdL    = CarlsonEllipticIntegral.rD(lambda,          x.add(lambda), y.add(lambda));
            double rdM    = CarlsonEllipticIntegral.rD(mu,              x.add(mu),     y.add(mu));
            double rd0    = CarlsonEllipticIntegral.rD(0, x,             y);
            double frac   = y.multiply(x.add(y).add(lambda).add(mu).sqrt()).reciprocal().multiply(3);
            Assert.assertEquals(0.0, rdL.add(rdM).subtract(rd0.subtract(frac)).norm(), 9.0e-12);
        }
    }

    @Test
    public void testRdNonSymmetry1() {
        RandomGenerator random = new Well19937c(0x66db170b5ee1afc2l);
        int countWrongRoot = 0;
        for (int i = 0; i < 10000; ++i) {
            double x = random.nextDouble( * 2 - 1, random.nextDouble() * 2 - 1);
            double y = random.nextDouble( * 2 - 1, random.nextDouble() * 2 - 1);
            double z = random.nextDouble( * 2 - 1, random.nextDouble() * 2 - 1);
            if (x.isZero() || y.isZero()) {
                continue;
            }
            // this is DLMF equation 19.21.7, computing square roots both after the fraction
            // (i.e. √x √y / √z) and before the fraction (i.e. √(xy/z))
            // the second form is used in DLMF as of 2021-06-04 and selects the wrong root
            // 25% of times when x, y, z are real and 33% of times when they are complex
            double lhs           = x.subtract(y).multiply(CarlsonEllipticIntegral.rD(y, z, x)).add(z.subtract(y).multiply(CarlsonEllipticIntegral.rD(x, y, z)));
            double rootGlobal    = y.divide(x.multiply(z)).sqrt();
            double rootSeparated = y.sqrt().divide(x.sqrt().multiply(z.sqrt()));
            double rhsGlobal     = CarlsonEllipticIntegral.rF(x, y, z).subtract(rootGlobal).multiply(3);
            double rhsSeparated  = CarlsonEllipticIntegral.rF(x, y, z).subtract(rootSeparated).multiply(3);
            if (lhs.subtract(rhsGlobal).norm() > 1.0e-3) {
                ++countWrongRoot;
                // when the wrong root is selected, the result is really bad
                Assert.assertTrue(lhs.subtract(rhsGlobal).norm() > 0.1);
            }
            Assert.assertEquals(0.0, lhs.subtract(rhsSeparated).norm(), 1.0e-10);
        }
        Assert.assertTrue(countWrongRoot > 3300);
    }

    @Test
    public void testRdNonSymmetry2() {
        RandomGenerator random = new Well19937c(0x1a8994acc807438dl);
        int countWrongRoot = 0;
        for (int i = 0; i < 10000; ++i) {
            double x = random.nextDouble( * 2 - 1, random.nextDouble() * 2 - 1);
            double y = random.nextDouble( * 2 - 1, random.nextDouble() * 2 - 1);
            double z = random.nextDouble( * 2 - 1, random.nextDouble() * 2 - 1);
            if (x.isZero() || y.isZero() || z.isZero()) {
                continue;
            }
            // this is DLMF equation 19.21.8, computing square roots both after the multiplication
            // (i.e. 1 / (√x √y √z)) and before the multiplication (i.e. 1 / √(xyz))
            // the second form is used in DLMF as of 2021-06-04 and selects the wrong root
            // 50% of times when x, y, z are real and 33% of times when they are complex
            double lhs           = CarlsonEllipticIntegral.rD(y, z, x).add(CarlsonEllipticIntegral.rD(z, x, y)).add(CarlsonEllipticIntegral.rD(x, y, z));
            double rootGlobal    = x.multiply(y.multiply(z)).sqrt();
            double rootSeparated = x.sqrt().multiply(y.sqrt().multiply(z.sqrt()));
            double rhsGlobal     = rootGlobal.reciprocal().multiply(3);
            double rhsSeparated  = rootSeparated.reciprocal().multiply(3);
            if (lhs.subtract(rhsGlobal).norm() > 1.0e-3) {
                ++countWrongRoot;
                // when the wrong root is selected, the result is really bad
                Assert.assertTrue(lhs.subtract(rhsGlobal).norm() > 3.0);
            }
            Assert.assertEquals(0.0, lhs.subtract(rhsSeparated).norm(), 2.0e-11);
        }
        Assert.assertTrue(countWrongRoot > 3300);
    }

    @Test
    public void testCarlson1995rG() {

        double rg1 = CarlsonEllipticIntegral.rG(0, 16, 16);
        Assert.assertEquals(FastMath.PI, 0.0, rg1, 1.0e-13);

        double rg2 = CarlsonEllipticIntegral.rG(2, 3, 4);
        Assert.assertEquals(1.7255030280692, 0.0, rg2, 1.0e-13);

        double rg3 = CarlsonEllipticIntegral.rG(0, buildComplex(0, 1), buildComplex(0, -1));
        Assert.assertEquals( 0.42360654239699, 0.0, rg3, 1.0e-13);

        double rg4 = CarlsonEllipticIntegral.rG(buildComplex(-1, 1), buildComplex(0, 1), 0);
        Assert.assertEquals(0.44660591677018, 0.70768352357515, rg4, 1.0e-13);

        double rg5 = CarlsonEllipticIntegral.rG(buildComplex(0, -1), buildComplex(-1, 1), buildComplex(0, 1));
        Assert.assertEquals(0.36023392184473, 0.40348623401722, rg5, 1.0e-13);

        double rg6 = CarlsonEllipticIntegral.rG(0, 0.0796, 4);
        Assert.assertEquals( 1.0284758090288, 0.0, rg6, 1.0e-13);

    }

    @Test
    public void testAlternateRG() {
        RandomGenerator random = new Well19937c(0xa2946e4a55d133a6l);
        for (int i = 0; i < 10000; ++i) {
            double x = random.nextDouble() * 3);
            double y = buildComplex(random.nextDouble() * 3);
            double z = buildComplex(random.nextDouble() * 3;
            Assert.assertEquals(0.0, CarlsonEllipticIntegral.rG(x, y, z).subtract(rgAlternateImplementation(x, y, z)).norm(), 2.0e-15);
        }
    }

    @Test
    public void testRgBuggySquareRoot() {

        // xy/z ≈ -0.566379 - 7.791 10⁻⁹ i ⇒ √(xy/z) ≈ 5.176 10⁻⁹ - 0.752582 i
        double x = buildComplex(FastMath.scalb(7745000, -24), -0.5625);
        double y = buildComplex(-0.3125, -0.6875);
        double z = buildComplex( 0.9375,  0.25);

        // on this side, all implementations match
        Assert.assertEquals(0.0,     CarlsonEllipticIntegral.rG(x, y, z).     subtract(rgAlternateImplementation(x, y, z)).norm(), 2.0e-16);
        Assert.assertEquals(0.0, buggyRG(x, y, z).subtract(rgAlternateImplementation(x, y, z)).norm(),     2.0e-16);

        // slightly shift x, so xy/z imaginary part changes sign
        // the selected square root also changes dramatically sign so implementation becomes wrong
        // xy/z ≈ -0.566379 + 2.807 10⁻⁸ i ⇒ √(xy/z) ≈ 1.865 10⁻⁸ + 0.752582 i
        x = buildComplex(FastMath.scalb(7744999, -24), -0.5625);
        Assert.assertEquals(0.0,     CarlsonEllipticIntegral.rG(x, y, z).     subtract(rgAlternateImplementation(x, y, z)).norm(), 2.0e-16);
        Assert.assertEquals(0.75258, buggyRG(x, y, z).subtract(rgAlternateImplementation(x, y, z)).norm(), 1.0e-5);

    }

    private double buggyRG(final double x, final double y, final double z) {
        final double termF = new RfFieldDuplication<>(x, y, z).integral().multiply(z);
        final double termD = x.subtract(z).multiply(y.subtract(z)).multiply(new RdFieldDuplication<>(x, y, z).integral()).divide(3);
        final double termS = x.multiply(y).divide(z).sqrt(); // ← the error is here, we must compute roots for each x, y and z before computing the fraction
        return termF.subtract(termD).add(termS).multiply(0.5);
    }

    private double rgAlternateImplementation(final double x, final double y, final double z) {
        // this implementation uses DLFM equation 19.21.11
        return d(x, y, z).add(d(y, z, x)).add(d(z, x, y)).divide(6);
    }

    private double d(final double u, final double v, final double w) {
        return u.isZero() ? u : u.multiply(v.add(w)).multiply(new RdFieldDuplication<>(v, w, u).integral());
    }

}
