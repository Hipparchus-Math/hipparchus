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
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

public class CarlsonEllipticIntegralRealTest {

    @Test
    public void testNoConvergenceRf() {
        Assert.assertTrue(Double.isNaN(CarlsonEllipticIntegral.rF(1, 2, Double.NaN)));
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
            double x      = random.nextDouble() * 3;
            double y      = random.nextDouble() * 3;
            double lambda = random.nextDouble() * 3;
            double mu     = x * y / lambda;
            double rfL    = CarlsonEllipticIntegral.rF(x + lambda, y + lambda, lambda);
            double rfM    = CarlsonEllipticIntegral.rF(x + mu,     y + mu,     mu);
            double rf0    = CarlsonEllipticIntegral.rF(x,             y,             0);
            Assert.assertEquals(0.0, FastMath.abs(rfL + rfM - rf0), 2.0e-14);
        }
    }

    @Test
    public void testNoConvergenceRc() {
        Assert.assertTrue(Double.isNaN(CarlsonEllipticIntegral.rC(1, Double.NaN)));
    }

    @Test
    public void testCarlson1995rC() {

        double rc1 = CarlsonEllipticIntegral.rC(0, 0.25);
        Assert.assertEquals(FastMath.PI, rc1, 1.0e-15);

        double rc2 = CarlsonEllipticIntegral.rC(2.25, 2);
        Assert.assertEquals(FastMath.log(2), rc2, 1.0e-15);

        double rc5 = CarlsonEllipticIntegral.rC(0.25, -2);
        Assert.assertEquals(FastMath.log(2) / 3.0, rc5, 1.0e-15);

    }

    @Test
    public void testCarlson1995ConsistencyRc() {
        RandomGenerator random = new Well19937c(0xf1170b6fc1a199cal);
        for (int i = 0; i < 10000; ++i) {
            double x      = random.nextDouble() * 3;
            double lambda = random.nextDouble() * 3;
            double mu     = x * x / lambda;
            double rcL    = CarlsonEllipticIntegral.rC(lambda,          x + lambda);
            double rcM    = CarlsonEllipticIntegral.rC(mu,              x + mu);
            double rc0    = CarlsonEllipticIntegral.rC(0, x);
            Assert.assertEquals(0.0, FastMath.abs(rcL + rcM - rc0), 3.0e-14);
        }
    }

    @Test
    public void testRfRc() {
        RandomGenerator random = new Well19937a(0x7e8041334a8c20edl);
        for (int i = 0; i < 10000; ++i) {
            final double x = 3 * random.nextDouble();
            final double y = 3 * random.nextDouble();
            final double rf = CarlsonEllipticIntegral.rF(x, y, y);
            final double rc = CarlsonEllipticIntegral.rC(x, y);
            Assert.assertEquals(0.0, FastMath.abs(rf - rc), 4.0e-15);
        }
    }

    @Test
    public void testNoConvergenceRj() {
        Assert.assertTrue(Double.isNaN(CarlsonEllipticIntegral.rJ(1, 1, 1, Double.NaN)));
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
            double x      = random.nextDouble() * 3;
            double y      = random.nextDouble() * 3;
            double p      = random.nextDouble() * 3;
            double lambda = random.nextDouble() * 3;
            double mu     = x * y / lambda;
            double a      = p * p * (lambda + mu + x + y);
            double b      = p * (p + lambda) * (p + mu);
            double rjL    = CarlsonEllipticIntegral.rJ(x + lambda, y + lambda, lambda,  p + lambda);
            double rjM    = CarlsonEllipticIntegral.rJ(x + mu,     y + mu,     mu,      p + mu);
            double rj0    = CarlsonEllipticIntegral.rJ(x,          y,          0,       p);
            double rc     = CarlsonEllipticIntegral.rC(a, b);
            Assert.assertEquals(0.0, FastMath.abs(rjL + rjM - (rj0 - rc * 3)), 3.0e-13);
        }
    }

    @Test
    public void testNoConvergenceRd() {
        Assert.assertTrue(Double.isNaN(CarlsonEllipticIntegral.rD(1, 1, Double.NaN)));
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
            double y      = random.nextDouble() * 3;
            double lambda = random.nextDouble() * 3;
            double mu     = x * y / lambda;
            double rdL    = CarlsonEllipticIntegral.rD(lambda,          x + lambda, y + lambda);
            double rdM    = CarlsonEllipticIntegral.rD(mu,              x + mu,     y + mu);
            double rd0    = CarlsonEllipticIntegral.rD(0,               x,          y);
            double frac   = 3 / (y * FastMath.sqrt(x + y + lambda + mu));
            Assert.assertEquals(0.0, FastMath.abs(rdL + rdM - rd0 + frac), 9.0e-12);
        }
    }

    @Test
    public void testRdNonSymmetry1() {
        RandomGenerator random = new Well19937c(0x66db170b5ee1afc2l);
        for (int i = 0; i < 10000; ++i) {
            double x = random.nextDouble();
            double y = random.nextDouble();
            double z = random.nextDouble();
            if (x == 0 || y == 0) {
                continue;
            }
            // this is DLMF equation 19.21.7
            double lhs = (x - y) * CarlsonEllipticIntegral.rD(y, z, x) + (z - y) * CarlsonEllipticIntegral.rD(x, y, z);
            double rhs = (CarlsonEllipticIntegral.rF(x, y, z) - FastMath.sqrt(y / (x * z))) * 3;
            Assert.assertEquals(0.0, FastMath.abs(lhs - rhs), 1.0e-10);
        }
    }

    @Test
    public void testRdNonSymmetry2() {
        RandomGenerator random = new Well19937c(0x1a8994acc807438dl);
        for (int i = 0; i < 10000; ++i) {
            double x = random.nextDouble();
            double y = random.nextDouble();
            double z = random.nextDouble();
            if (x == 0 || y == 0 || z == 0) {
                continue;
            }
            // this is DLMF equation 19.21.8
            double lhs = CarlsonEllipticIntegral.rD(y, z, x) + CarlsonEllipticIntegral.rD(z, x, y) + CarlsonEllipticIntegral.rD(x, y, z);
            double rhs = 3 / FastMath.sqrt(x * y * z);
            Assert.assertEquals(0.0, FastMath.abs(lhs - rhs), 2.0e-11);
        }
    }

    @Test
    public void testCarlson1995rG() {

        double rg1 = CarlsonEllipticIntegral.rG(0, 16, 16);
        Assert.assertEquals(FastMath.PI, rg1, 1.0e-13);

        double rg2 = CarlsonEllipticIntegral.rG(2, 3, 4);
        Assert.assertEquals(1.7255030280692, rg2, 1.0e-13);

        double rg6 = CarlsonEllipticIntegral.rG(0, 0.0796, 4);
        Assert.assertEquals( 1.0284758090288, rg6, 1.0e-13);

    }

    @Test
    public void testAlternateRG() {
        RandomGenerator random = new Well19937c(0xa2946e4a55d133a6l);
        for (int i = 0; i < 10000; ++i) {
            double x = random.nextDouble() * 3;
            double y = random.nextDouble() * 3;
            double z = random.nextDouble() * 3;
            Assert.assertEquals(0.0, FastMath.abs(CarlsonEllipticIntegral.rG(x, y, z) - rgAlternateImplementation(x, y, z)), 2.0e-15);
        }
    }

    private double rgAlternateImplementation(final double x, final double y, final double z) {
        // this implementation uses DLFM equation 19.21.11
        return (d(x, y, z) + d(y, z, x) + d(z, x, y)) / 6;
    }

    private double d(final double u, final double v, final double w) {
        return u == 0 ? u : u * (v + w) * (new RdRealDuplication(v, w, u).integral());
    }

}
