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
package org.hipparchus.special.elliptic.jacobi;

import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

public class JacobiEllipticTest {

    @Test
    public void testCircular() {
        for (double m : new double[] { -1.0e-10, 0.0, 1.0e-10 }) {
         final double eps = 3 * FastMath.max(1.0e-14, FastMath.abs(m));
            final JacobiElliptic je = JacobiEllipticBuilder.build(m);
            for (double t = -10; t < 10; t += 0.01) {
                final CopolarN n = je.valuesN(t);
                Assert.assertEquals(FastMath.sin(t), n.sn(), eps);
                Assert.assertEquals(FastMath.cos(t), n.cn(), eps);
                Assert.assertEquals(1.0,             n.dn(), eps);
            }
        }
    }

    @Test
    public void testHyperbolic() {
        for (double m1 : new double[] { -1.0e-12, 0.0, 1.0e-12 }) {
            final double eps = 3 * FastMath.max(1.0e-14, FastMath.abs(m1));
            final JacobiElliptic je = JacobiEllipticBuilder.build(1.0 - m1);
            for (double t = -3; t < 3; t += 0.01) {
                final CopolarN n = je.valuesN(t);
                Assert.assertEquals(FastMath.tanh(t),       n.sn(), eps);
                Assert.assertEquals(1.0 / FastMath.cosh(t), n.cn(), eps);
                Assert.assertEquals(1.0 / FastMath.cosh(t), n.dn(), eps);
            }
        }
    }

    @Test
    public void testNoConvergence() {
        Assert.assertTrue(Double.isNaN(JacobiEllipticBuilder.build(Double.NaN).valuesS(0.0).cs()));
    }

    @Test
    public void testNegativeParameter() {
        Assert.assertEquals(0.49781366219021166315, JacobiEllipticBuilder.build(-4.5).valuesN(8.3).sn(), 1.5e-10);
        Assert.assertEquals(0.86728401215332559984, JacobiEllipticBuilder.build(-4.5).valuesN(8.3).cn(), 1.5e-10);
        Assert.assertEquals(1.45436686918553524215, JacobiEllipticBuilder.build(-4.5).valuesN(8.3).dn(), 1.5e-10);
    }

    @Test
    public void testAbramowitzStegunExample1() {
        // Abramowitz and Stegun give a result of -1667, but Wolfram Alpha gives the following value
        Assert.assertEquals(-1392.11114434139393839735, JacobiEllipticBuilder.build(0.64).valuesC(1.99650).nc(), 6.0e-10);
    }

    @Test
    public void testAbramowitzStegunExample2() {
        Assert.assertEquals(0.996253, JacobiEllipticBuilder.build(0.19).valuesN(0.20).dn(), 1.0e-6);
    }

    @Test
    public void testAbramowitzStegunExample3() {
        Assert.assertEquals(0.984056, JacobiEllipticBuilder.build(0.81).valuesN(0.20).dn(), 1.0e-6);
    }

    @Test
    public void testAbramowitzStegunExample4() {
        Assert.assertEquals(0.980278, JacobiEllipticBuilder.build(0.81).valuesN(0.20).cn(), 1.0e-6);
    }

    @Test
    public void testAbramowitzStegunExample5() {
        Assert.assertEquals(0.60952, JacobiEllipticBuilder.build(0.36).valuesN(0.672).sn(), 1.0e-5);
        Assert.assertEquals(1.1740, JacobiEllipticBuilder.build(0.36).valuesC(0.672).dc(), 1.0e-4);
    }

    @Test
    public void testAbramowitzStegunExample7() {
        Assert.assertEquals(1.6918083, JacobiEllipticBuilder.build(0.09).valuesS(0.5360162).cs(), 1.0e-7);
    }

    @Test
    public void testAbramowitzStegunExample8() {
        Assert.assertEquals(0.56458, JacobiEllipticBuilder.build(0.5).valuesN(0.61802).sn(), 1.0e-5);
    }

    @Test
    public void testAbramowitzStegunExample9() {
        Assert.assertEquals(0.68402, JacobiEllipticBuilder.build(0.5).valuesC(0.61802).sc(), 1.0e-5);
    }

    @Test
    public void testAllFunctions() {
        // reference was computed from Wolfram Alpha, using the square relations
        // from Abramowitz and Stegun section 16.9 for the functions Wolfram Alpha
        // did not understood (i.e. for the sake of validation we did *not* use the
        // relations from section 16.3 which are used in the Hipparchus implementation)
        final double u = 1.4;
        final double m = 0.7;
        final double[] reference = {
              0.92516138673582827365, 0.37957398289798418747, 0.63312991237590996850,
              0.41027866958131945027, 0.68434537093007175683, 1.08089249544689572795,
              1.66800134071905681841, 2.63453251554589286796, 2.43736775548306830513,
              1.57945467502452678756, 1.46125047743207819361, 0.59951990180590090343
        };
        final JacobiElliptic je = JacobiEllipticBuilder.build(m);
        Assert.assertEquals(reference[ 0], je.valuesN(u).sn(), 4 * FastMath.ulp(reference[ 0]));
        Assert.assertEquals(reference[ 1], je.valuesN(u).cn(), 4 * FastMath.ulp(reference[ 1]));
        Assert.assertEquals(reference[ 2], je.valuesN(u).dn(), 4 * FastMath.ulp(reference[ 2]));
        Assert.assertEquals(reference[ 3], je.valuesS(u).cs(), 4 * FastMath.ulp(reference[ 3]));
        Assert.assertEquals(reference[ 4], je.valuesS(u).ds(), 4 * FastMath.ulp(reference[ 4]));
        Assert.assertEquals(reference[ 5], je.valuesS(u).ns(), 4 * FastMath.ulp(reference[ 5]));
        Assert.assertEquals(reference[ 6], je.valuesC(u).dc(), 4 * FastMath.ulp(reference[ 6]));
        Assert.assertEquals(reference[ 7], je.valuesC(u).nc(), 4 * FastMath.ulp(reference[ 7]));
        Assert.assertEquals(reference[ 8], je.valuesC(u).sc(), 4 * FastMath.ulp(reference[ 8]));
        Assert.assertEquals(reference[ 9], je.valuesD(u).nd(), 4 * FastMath.ulp(reference[ 9]));
        Assert.assertEquals(reference[10], je.valuesD(u).sd(), 4 * FastMath.ulp(reference[10]));
        Assert.assertEquals(reference[11], je.valuesD(u).cd(), 4 * FastMath.ulp(reference[11]));
    }

    @Test
    public void testInverseCopolarN() {
        final double m = 0.7;
        final JacobiElliptic je = JacobiEllipticBuilder.build(m);
        doTestInverse(-0.80,  0.80, 100, u -> je.valuesN(u).sn(), x -> je.arcsn(x), 1.0e-14);
        doTestInverse(-1.00,  1.00, 100, u -> je.valuesN(u).cn(), x -> je.arccn(x), 1.0e-14);
        doTestInverse( 0.55,  1.00, 100, u -> je.valuesN(u).dn(), x -> je.arcdn(x), 1.0e-14);
    }

    @Test
    public void testInverseCopolarS() {
        final double m = 0.7;
        final JacobiElliptic je = JacobiEllipticBuilder.build(m);
        doTestInverse(-2.00,  2.00, 100, u -> je.valuesS(u).cs(), x -> je.arccs(x), 1.0e-14);
        doTestInverse( 0.55,  2.00, 100, u -> je.valuesS(u).ds(), x -> je.arcds(x), 1.0e-14);
        doTestInverse(-2.00, -0.55, 100, u -> je.valuesS(u).ds(), x -> je.arcds(x), 1.0e-14);
        doTestInverse( 1.00,  2.00, 100, u -> je.valuesS(u).ns(), x -> je.arcns(x), 1.0e-11);
        doTestInverse(-2.00, -1.00, 100, u -> je.valuesS(u).ns(), x -> je.arcns(x), 1.0e-11);
    }

    @Test
    public void testInverseCopolarC() {
        final double m = 0.7;
        final JacobiElliptic je = JacobiEllipticBuilder.build(m);
        doTestInverse( 1.00,  2.00, 100, u -> je.valuesC(u).dc(), x -> je.arcdc(x), 1.0e-14);
        doTestInverse(-2.00, -1.00, 100, u -> je.valuesC(u).dc(), x -> je.arcdc(x), 1.0e-14);
        doTestInverse( 1.00,  2.00, 100, u -> je.valuesC(u).nc(), x -> je.arcnc(x), 1.0e-14);
        doTestInverse(-2.00, -1.00, 100, u -> je.valuesC(u).nc(), x -> je.arcnc(x), 1.0e-14);
        doTestInverse(-2.00,  2.00, 100, u -> je.valuesC(u).sc(), x -> je.arcsc(x), 1.0e-14);
    }

    @Test
    public void testInverseCopolarD() {
        final double m = 0.7;
        final JacobiElliptic je = JacobiEllipticBuilder.build(m);
        doTestInverse( 1.00,  1.80, 100, u -> je.valuesD(u).nd(), x -> je.arcnd(x), 1.0e-14);
        doTestInverse(-1.80,  1.80, 100, u -> je.valuesD(u).sd(), x -> je.arcsd(x), 1.0e-14);
        doTestInverse(-1.00,  1.00, 100, u -> je.valuesD(u).cd(), x -> je.arccd(x), 1.0e-14);
    }

    private void doTestInverse(final double xMin, final double xMax, final int n,
                               final UnivariateFunction direct, final UnivariateFunction inverse,
                               final double tolerance) {
        for (int i = 0; i < n; ++i) {
            final double x        = xMin + i * (xMax - xMin) / (n - 1);
            final double xRebuilt = direct.value(inverse.value(x));
            Assert.assertEquals(x, xRebuilt, tolerance);
        }
    }

}
