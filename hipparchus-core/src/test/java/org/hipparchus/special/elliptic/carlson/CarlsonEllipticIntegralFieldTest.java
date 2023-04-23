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
package org.hipparchus.special.elliptic.carlson;

import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well19937a;
import org.hipparchus.random.Well19937c;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

public class CarlsonEllipticIntegralFieldTest {

    private Binary64 build(double d) {
        return new Binary64(d);
    }

    @Test
    public void testNoConvergenceRf() {
        Assert.assertTrue(CarlsonEllipticIntegral.rF(build(1), build(2), build(Double.NaN)).isNaN());
    }

    @Test
    public void testDlmfRf() {
        Binary64 rf = CarlsonEllipticIntegral.rF(build(1), build(2), build(4));
        Assert.assertEquals(0.6850858166, rf.getReal(), 1.0e-10);
    }

    @Test
    public void testCarlson1995rF() {

        Binary64 rf1 = CarlsonEllipticIntegral.rF(build(1), build(2), build(0));
        Assert.assertEquals( 1.3110287771461, rf1.getReal(), 1.0e-13);

        Binary64 rf2 = CarlsonEllipticIntegral.rF(build(0.5), build(1), build(0));
        Assert.assertEquals( 1.8540746773014, rf2.getReal(), 1.0e-13);

        Binary64 rf4 = CarlsonEllipticIntegral.rF(build(2), build(3), build(4));
        Assert.assertEquals( 0.58408284167715, rf4.getReal(), 1.0e-13);

    }

    @Test
    public void testCarlson1995ConsistencyRf() {
        RandomGenerator random = new Well19937c(0x57f2689b3f4028b4l);
        for (int i = 0; i < 10000; ++i) {
            Binary64 x      = build(random.nextDouble() * 3);
            Binary64 y      = build(random.nextDouble() * 3);
            Binary64 lambda = build(random.nextDouble() * 3);
            Binary64 mu     = x.multiply(y).divide(lambda);
            Binary64 rfL    = CarlsonEllipticIntegral.rF(x.add(lambda), y.add(lambda), lambda);
            Binary64 rfM    = CarlsonEllipticIntegral.rF(x.add(mu),     y.add(mu),     mu);
            Binary64 rf0    = CarlsonEllipticIntegral.rF(x,             y,             Binary64.ZERO);
            Assert.assertEquals(0.0, rfL.add(rfM).subtract(rf0).norm(), 2.0e-14);
        }
    }

    @Test
    public void testNoConvergenceRc() {
        Assert.assertTrue(CarlsonEllipticIntegral.rC(build(1), build(Double.NaN)).isNaN());
    }

    @Test
    public void testCarlson1995rC() {

        Binary64 rc1 = CarlsonEllipticIntegral.rC(build(0), build(0.25));
        Assert.assertEquals(FastMath.PI, rc1.getReal(), 1.0e-15);

        Binary64 rc2 = CarlsonEllipticIntegral.rC(build(2.25), build(2));
        Assert.assertEquals(FastMath.log(2), rc2.getReal(), 1.0e-15);

        Binary64 rc5 = CarlsonEllipticIntegral.rC(build(0.25), build(-2));
        Assert.assertEquals(FastMath.log(2) / 3.0, rc5.getReal(), 1.0e-15);

    }

    @Test
    public void testCarlson1995ConsistencyRc() {
        RandomGenerator random = new Well19937c(0xf1170b6fc1a199cal);
        for (int i = 0; i < 10000; ++i) {
            Binary64 x      = build(random.nextDouble() * 3);
            Binary64 lambda = build(random.nextDouble() * 3);
            Binary64 mu     = x.multiply(x).divide(lambda);
            Binary64 rcL    = CarlsonEllipticIntegral.rC(lambda,         x.add(lambda));
            Binary64 rcM    = CarlsonEllipticIntegral.rC(mu,             x.add(mu));
            Binary64 rc0    = CarlsonEllipticIntegral.rC(Binary64.ZERO, x);
            Assert.assertEquals(0.0, rcL.add(rcM).subtract(rc0).norm(), 3.0e-14);
        }
    }

    @Test
    public void testRfRc() {
        RandomGenerator random = new Well19937a(0x7e8041334a8c20edl);
        for (int i = 0; i < 10000; ++i) {
            final Binary64 x = build(3 * random.nextDouble());
            final Binary64 y = build(3 * random.nextDouble());
            final Binary64 rf = CarlsonEllipticIntegral.rF(x, y, y);
            final Binary64 rc = CarlsonEllipticIntegral.rC(x, y);
            Assert.assertEquals(0.0, rf.subtract(rc).norm(), 4.0e-15);
        }
    }

    @Test
    public void testNoConvergenceRj() {
        Assert.assertTrue(CarlsonEllipticIntegral.rJ(build(1), build(1), build(1), build(Double.NaN)).isNaN());
    }

    @Test
    public void testCarlson1995rJ() {

        Binary64 rj01 = CarlsonEllipticIntegral.rJ(build(0), build(1), build(2), build(3));
        Assert.assertEquals(0.77688623778582, rj01.getReal(), 1.0e-13);

        Binary64 rj02 = CarlsonEllipticIntegral.rJ(build(2), build(3), build(4), build(5));
        Assert.assertEquals( 0.14297579667157, rj02.getReal(), 1.0e-13);

    }

    @Test
    public void testCarlson1995ConsistencyRj() {
        RandomGenerator random = new Well19937c(0x4af7bb722712e64el);
        for (int i = 0; i < 10000; ++i) {
            Binary64 x      = build(random.nextDouble() * 3);
            Binary64 y      = build(random.nextDouble() * 3);
            Binary64 p      = build(random.nextDouble() * 3);
            Binary64 lambda = build(random.nextDouble() * 3);
            Binary64 mu     = x.multiply(y).divide(lambda);
            Binary64 a      = p.multiply(p).multiply(lambda.add(mu).add(x).add(y));
            Binary64 b      = p.multiply(p.add(lambda)).multiply(p.add(mu));
            Binary64 rjL    = CarlsonEllipticIntegral.rJ(x.add(lambda), y.add(lambda), lambda,         p.add(lambda));
            Binary64 rjM    = CarlsonEllipticIntegral.rJ(x.add(mu),     y.add(mu),     mu,             p.add(mu));
            Binary64 rj0    = CarlsonEllipticIntegral.rJ(x,             y,             Binary64.ZERO, p);
            Binary64 rc     = CarlsonEllipticIntegral.rC(a, b);
            Assert.assertEquals(0.0, rjL.add(rjM).subtract(rj0.subtract(rc.multiply(3))).norm(), 3.0e-13);
        }
    }

    @Test
    public void testNoConvergenceRd() {
        Assert.assertTrue(CarlsonEllipticIntegral.rD(build(1), build(1), build(Double.NaN)).isNaN());
    }

    @Test
    public void testCarlson1995rD() {

        Binary64 rd1 = CarlsonEllipticIntegral.rD(build(0), build(2), build(1));
        Assert.assertEquals(1.7972103521034, rd1.getReal(), 1.0e-13);

        Binary64 rd2 = CarlsonEllipticIntegral.rD(build(2), build(3), build(4));
        Assert.assertEquals( 0.16510527294261, rd2.getReal(), 1.0e-13);

    }

    @Test
    public void testCarlson1995ConsistencyRd() {
        RandomGenerator random = new Well19937c(0x17dea97eeb78206al);
        for (int i = 0; i < 10000; ++i) {
            Binary64 x      = build(random.nextDouble() * 3);
            Binary64 y      = build(random.nextDouble() * 3);
            Binary64 lambda = build(random.nextDouble() * 3);
            Binary64 mu     = x.multiply(y).divide(lambda);
            Binary64 rdL    = CarlsonEllipticIntegral.rD(lambda,         x.add(lambda), y.add(lambda));
            Binary64 rdM    = CarlsonEllipticIntegral.rD(mu,             x.add(mu),     y.add(mu));
            Binary64 rd0    = CarlsonEllipticIntegral.rD(Binary64.ZERO, x,             y);
            Binary64 frac   = y.multiply(x.add(y).add(lambda).add(mu).sqrt()).reciprocal().multiply(3);
            Assert.assertEquals(0.0, rdL.add(rdM).subtract(rd0.subtract(frac)).norm(), 9.0e-12);
        }
    }

    @Test
    public void testRdNonSymmetry1() {
        RandomGenerator random = new Well19937c(0x66db170b5ee1afc2l);
        for (int i = 0; i < 10000; ++i) {
            Binary64 x = build(random.nextDouble());
            Binary64 y = build(random.nextDouble());
            Binary64 z = build(random.nextDouble());
            if (x.isZero() || y.isZero()) {
                continue;
            }
            // this is DLMF equation 19.21.7
            Binary64 lhs = x.subtract(y).multiply(CarlsonEllipticIntegral.rD(y, z, x)).
                            add(z.subtract(y).multiply(CarlsonEllipticIntegral.rD(x, y, z)));
            Binary64 rhs = CarlsonEllipticIntegral.rF(x, y, z).subtract(y.divide(x.multiply(z)).sqrt()).multiply(3);
            Assert.assertEquals(0.0, lhs.subtract(rhs).norm(), 1.0e-10);
        }
    }

    @Test
    public void testRdNonSymmetry2() {
        RandomGenerator random = new Well19937c(0x1a8994acc807438dl);
        for (int i = 0; i < 10000; ++i) {
            Binary64 x = build(random.nextDouble());
            Binary64 y = build(random.nextDouble());
            Binary64 z = build(random.nextDouble());
            if (x.isZero() || y.isZero() || z.isZero()) {
                continue;
            }
            // this is DLMF equation 19.21.8
            Binary64 lhs = CarlsonEllipticIntegral.rD(y, z, x).
                            add(CarlsonEllipticIntegral.rD(z, x, y)).
                            add(CarlsonEllipticIntegral.rD(x, y, z));
            Binary64 rhs = x.multiply(y.multiply(z)).sqrt().reciprocal().multiply(3);
            Assert.assertEquals(0.0, lhs.subtract(rhs).norm(), 2.0e-11);
        }
    }

    @Test
    public void testCarlson1995rG() {

        Binary64 rg1 = CarlsonEllipticIntegral.rG(build(0), build(16), build(16));
        Assert.assertEquals(FastMath.PI, rg1.getReal(), 1.0e-13);

        Binary64 rg2 = CarlsonEllipticIntegral.rG(build(2), build(3), build(4));
        Assert.assertEquals(1.7255030280692, rg2.getReal(), 1.0e-13);

        Binary64 rg6 = CarlsonEllipticIntegral.rG(build(0), build(0.0796), build(4));
        Assert.assertEquals( 1.0284758090288, rg6.getReal(), 1.0e-13);

    }

    @Test
    public void testAlternateRG() {
        RandomGenerator random = new Well19937c(0xa2946e4a55d133a6l);
        for (int i = 0; i < 10000; ++i) {
            Binary64 x = build(random.nextDouble() * 3);
            Binary64 y = build(random.nextDouble() * 3);
            Binary64 z = build(random.nextDouble() * 3);
            Assert.assertEquals(0.0, CarlsonEllipticIntegral.rG(x, y, z).subtract(rgAlternateImplementation(x, y, z)).norm(), 2.0e-15);
        }
    }

    private Binary64 rgAlternateImplementation(final Binary64 x, final Binary64 y, final Binary64 z) {
        // this implementation uses DLFM equation 19.21.11
        return d(x, y, z).add(d(y, z, x)).add(d(z, x, y)).divide(6);
    }

    private Binary64 d(final Binary64 u, final Binary64 v, final Binary64 w) {
        return u.isZero() ? u : u.multiply(v.add(w)).multiply(new RdFieldDuplication<>(v, w, u).integral());
    }

}
