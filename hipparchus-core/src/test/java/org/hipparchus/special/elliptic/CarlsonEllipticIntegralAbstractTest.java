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
package org.hipparchus.special.elliptic;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well19937a;
import org.hipparchus.random.Well19937c;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

public abstract class CarlsonEllipticIntegralAbstractTest<T extends CalculusFieldElement<T>> {

    protected abstract T buildComplex(double realPart);
    protected abstract T buildComplex(double realPart, double imaginaryPart);
    protected abstract T rF(T x, T y, T z);
    protected abstract T rC(T x, T y);
    protected abstract T rJ(T x, T y, T z, T p);
    protected abstract T rD(T x, T y, T z);

    private void check(double expectedReal, double expectedImaginary, T result, double tol) {
        Assert.assertEquals(0, buildComplex(expectedReal, expectedImaginary).subtract(result).norm(), tol);
    }

    @Test
    public void testNoConvergenceRf() {
        try {
            rF(buildComplex(1), buildComplex(2), buildComplex(Double.NaN));
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalStateException mise) {
            Assert.assertEquals(LocalizedCoreFormats.CONVERGENCE_FAILED, mise.getSpecifier());
        }
    }

    @Test
    public void testDlmfRf() {
        T rf = rF(buildComplex(1), buildComplex(2), buildComplex(4));
        check(0.6850858166, 0.0, rf, 1.0e-10);
    }

    @Test
    public void testCarlson1995rF() {

        T rf1 = rF(buildComplex(1), buildComplex(2), buildComplex(0));
        check( 1.3110287771461, 0.0, rf1, 1.0e-13);

        T rf2 = rF(buildComplex(0.5), buildComplex(1), buildComplex(0));
        check( 1.8540746773014, 0.0, rf2, 1.0e-13);

        T rf3 = rF(buildComplex(-1, 1), buildComplex(0, 1), buildComplex(0));
        check( 0.79612586584234, -1.2138566698365, rf3, 1.0e-13);

        T rf4 = rF(buildComplex(2), buildComplex(3), buildComplex(4));
        check( 0.58408284167715, 0.0, rf4, 1.0e-13);

        T rf5 = rF(buildComplex(0, 1), buildComplex(0, -1), buildComplex(2));
        check( 1.0441445654064, 0.0, rf5, 1.0e-13);

        T rf6 = rF(buildComplex(-1, 1), buildComplex(0, 1), buildComplex(1, -1));
        check( 0.93912050218619, -0.53296252018635, rf6, 1.0e-13);

    }

    @Test
    public void testCarlson1995ConsistencyRf() {
        RandomGenerator random = new Well19937c(0x57f2689b3f4028b4l);
        for (int i = 0; i < 100000; ++i) {
            T x      = buildComplex(random.nextDouble() * 3);
            T y      = buildComplex(random.nextDouble() * 3);
            T lambda = buildComplex(random.nextDouble() * 6 - 3, random.nextDouble() * 3);
            T mu     = x.multiply(y).divide(lambda);
            T rfL    = rF(x.add(lambda), y.add(lambda), lambda);
            T rfM    = rF(x.add(mu),     y.add(mu),     mu);
            T rf0    = rF(x,             y,             buildComplex(0));
            Assert.assertEquals(0.0, rfL.add(rfM).subtract(rf0).norm(), 2.0e-14);
        }
    }

    @Test
    public void testNoConvergenceRc() {
        try {
            rC(buildComplex(1), buildComplex(Double.NaN));
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalStateException mise) {
            Assert.assertEquals(LocalizedCoreFormats.CONVERGENCE_FAILED, mise.getSpecifier());
        }
    }

    @Test
    public void testCarlson1995rC() {

        T rc1 = rC(buildComplex(0), buildComplex(0.25));
        check(FastMath.PI, 0.0, rc1, 1.0e-15);

        T rc2 = rC(buildComplex(2.25), buildComplex(2));
        check(FastMath.log(2), 0.0, rc2, 1.0e-15);

        T rc3 = rC(buildComplex(0), buildComplex(0, 1));
        check( 1.1107207345396, -1.1107207345396, rc3, 1.0e-13);

        T rc4 = rC(buildComplex(0, -1), buildComplex(0, 1));
        check( 1.2260849569072, -0.34471136988768, rc4, 1.0e-13);

        // the 1995 paper does not show any imaginary part for this case,
        // but it seems to be non-zero, which is consistent with R_F(x, y, y)
        T rc5 = rC(buildComplex(0.25), buildComplex(-2));
        check(FastMath.log(2) / 3.0, -1.047197551196598, rc5, 1.0e-15);

        // there is probably a confusion in the 1995 paper as the following
        // results were expected to be R_C(i, -1), but we found it was rather R_C(2, 1 - i)
        // or R_C(4 + 3i, -i)
        T rc6A = rC(buildComplex(2), buildComplex(1, -1));
        check( 0.77778596920447, 0.19832484993429, rc6A, 1.0e-13);
        T rc6B = rC(buildComplex(4, 3), buildComplex(0, -1));
        check( 0.77778596920447, 0.19832484993429, rc6B, 1.0e-13);

    }

    @Test
    public void testCarlson1995ConsistencyRc() {
        RandomGenerator random = new Well19937c(0xf1170b6fc1a199cal);
        for (int i = 0; i < 100000; ++i) {
            T x      = buildComplex(random.nextDouble() * 3);
            T lambda = buildComplex(random.nextDouble() * 6 - 3, random.nextDouble() * 3);
            T mu     = x.multiply(x).divide(lambda);
            T rcL    = rC(lambda,          x.add(lambda));
            T rcM    = rC(mu,              x.add(mu));
            T rc0    = rC(buildComplex(0), x);
            Assert.assertEquals(0.0, rcL.add(rcM).subtract(rc0).norm(), 2.0e-13);
        }
    }

    @Test
    public void testRfRc() {
        RandomGenerator random = new Well19937a(0x7e8041334a8c20edl);
        for (int i = 0; i < 100000; ++i) {
            final T x = buildComplex(6 * random.nextDouble() - 3,
                                          6 * random.nextDouble() - 3);
            final T y = buildComplex(6 * random.nextDouble() - 3,
                                          6 * random.nextDouble() - 3);
            final T rf = rF(x, y, y);
            final T rc = rC(x, y);
            Assert.assertEquals(0.0, rf.subtract(rc).norm(), 3.0e-14 * rf.norm());
        }
    }

    @Test
    public void testNoConvergenceRj() {
        try {
            rJ(buildComplex(1), buildComplex(1), buildComplex(1), buildComplex(Double.NaN));
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalStateException mise) {
            Assert.assertEquals(LocalizedCoreFormats.CONVERGENCE_FAILED, mise.getSpecifier());
        }
    }

    @Test
    public void testCarlson1995rJ() {

        T rj01 = rJ(buildComplex(0), buildComplex(1), buildComplex(2), buildComplex(3));
        check(0.77688623778582, 0.0, rj01, 1.0e-13);

        T rj02 = rJ(buildComplex(2), buildComplex(3), buildComplex(4), buildComplex(5));
        check( 0.14297579667157, 0.0, rj02, 1.0e-13);

        T rj03 = rJ(buildComplex(2), buildComplex(3), buildComplex(4), buildComplex(-1, 1));
        check( 0.13613945827771, -0.38207561624427, rj03, 1.0e-13);

        T rj04 = rJ(buildComplex(0, 1), buildComplex(0, -1), buildComplex(0), buildComplex(2));
        check( 1.6490011662711, 0.0, rj04, 1.0e-13);

        T rj05 = rJ(buildComplex(-1, 1), buildComplex(-1, -1), buildComplex(1), buildComplex(2));
        check( 0.94148358841220, 0.0, rj05, 1.0e-13);

        T rj06 = rJ(buildComplex(0, 1), buildComplex(0, -1), buildComplex(0), buildComplex(1, -1));
        check( 1.8260115229009, 1.2290661908643, rj06, 1.0e-13);

        T rj07 = rJ(buildComplex(-1, 1), buildComplex(-1, -1), buildComplex(1), buildComplex(-3, 1));
        check(-0.61127970812028, -1.0684038390007, rj07, 1.0e-13);

        T rj08 = rJ(buildComplex(-1, 1), buildComplex(-2, -1), buildComplex(0, -1), buildComplex(-1, 1));
        check( 1.8249027393704, -1.2218475784827, rj08, 1.0e-13);

        T rj09 = rJ(buildComplex(2), buildComplex(3), buildComplex(4), buildComplex(-0.5));
        check( 0.24723819703052, -0.7509842836891, rj09, 1.0e-13);

        T rj10 = rJ(buildComplex(2), buildComplex(3), buildComplex(4), buildComplex(-5));
        check(-0.12711230042964, -0.2099064885453, rj10, 1.0e-13);

    }

    @Test
    public void testCarlson1995ConsistencyRj() {
        RandomGenerator random = new Well19937c(0x4af7bb722712e64el);
        for (int i = 0; i < 100000; ++i) {
            T x      = buildComplex(random.nextDouble() * 3);
            T y      = buildComplex(random.nextDouble() * 3);
            T p      = buildComplex(random.nextDouble() * 3);
            T lambda = buildComplex(random.nextDouble() * 6 - 3, random.nextDouble() * 3);
            T mu     = x.multiply(y).divide(lambda);
            T a      = p.multiply(p).multiply(lambda.add(mu).add(x).add(y));
            T b      = p.multiply(p.add(lambda)).multiply(p.add(mu));
            T rjL    = rJ(x.add(lambda), y.add(lambda), lambda,          p.add(lambda));
            T rjM    = rJ(x.add(mu),     y.add(mu),     mu,              p.add(mu));
            T rj0    = rJ(x,             y,             buildComplex(0), p);
            T rc     = rC(a, b);
            Assert.assertEquals(0.0, rjL.add(rjM).subtract(rj0.subtract(rc.multiply(3))).norm(), 2.0e-12);
        }
    }

    @Test
    public void testNoConvergenceRd() {
        try {
            rD(buildComplex(1), buildComplex(1), buildComplex(Double.NaN));
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalStateException mise) {
            Assert.assertEquals(LocalizedCoreFormats.CONVERGENCE_FAILED, mise.getSpecifier());
        }
    }

    @Test
    public void testCarlson1995rD() {

        T rd1 = rD(buildComplex(0), buildComplex(2), buildComplex(1));
        check(1.7972103521034, 0.0, rd1, 1.0e-13);

        T rd2 = rD(buildComplex(2), buildComplex(3), buildComplex(4));
        check( 0.16510527294261, 0.0, rd2, 1.0e-13);

        T rd3 = rD(buildComplex(0, 1), buildComplex(0, -1), buildComplex(2));
        check( 0.65933854154220, 0.0, rd3, 1.0e-13);

        T rd4 = rD(buildComplex(0), buildComplex(0, 1), buildComplex(0, -1));
        check( 1.2708196271910, 2.7811120159521, rd4, 1.0e-13);

        T rd5 = rD(buildComplex(0), buildComplex(-1, 1), buildComplex(0, 1));
        check(-1.8577235439239, -0.96193450888830, rd5, 1.0e-13);

        T rd6 = rD(buildComplex(-2, -1), buildComplex(0, -1), buildComplex(-1, 1));
        check( 1.8249027393704, -1.2218475784827, rd6, 1.0e-13);

    }

    @Test
    public void testCarlson1995ConsistencyRd() {
        RandomGenerator random = new Well19937c(0x17dea97eeb78206al);
        for (int i = 0; i < 100000; ++i) {
            T x      = buildComplex(random.nextDouble() * 3);
            T y      = buildComplex(random.nextDouble() * 3);
            T lambda = buildComplex(random.nextDouble() * 6 - 3, random.nextDouble() * 3);
            T mu     = x.multiply(y).divide(lambda);
            T rdL    = rD(lambda,          x.add(lambda), y.add(lambda));
            T rdM    = rD(mu,              x.add(mu),     y.add(mu));
            T rd0    = rD(buildComplex(0), x,             y);
            T frac   = y.multiply(x.add(y).add(lambda).add(mu).sqrt()).reciprocal().multiply(3);
            Assert.assertEquals(0.0, rdL.add(rdM).subtract(rd0.subtract(frac)).norm(), 3.0e-9);
        }
    }

}
