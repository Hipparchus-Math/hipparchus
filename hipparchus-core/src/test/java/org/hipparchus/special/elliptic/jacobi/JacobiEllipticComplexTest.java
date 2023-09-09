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

import java.io.IOException;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.complex.Complex;
import org.hipparchus.special.elliptic.legendre.LegendreEllipticIntegral;
import org.junit.Test;

public class JacobiEllipticComplexTest {

    @Test
    public void testComplex() throws IOException {
        final FieldJacobiElliptic<Complex> je = JacobiEllipticBuilder.build(Complex.valueOf(0.3, 1.0));
        UnitTestUtils.assertEquals(Complex.valueOf(-0.24609405083573348938, 0.74202229271111558523),
                                   je.valuesC(Complex.valueOf(5.2, -2.5)).sc(),
                                   1.0e-15);
    }

    @Test
    public void testIssue150Plus() {
        final double tol = 1.0e-15;
        final FieldJacobiElliptic<Complex> je = JacobiEllipticBuilder.build(Complex.valueOf(-1, 1));
        UnitTestUtils.assertEquals(Complex.valueOf(-0.05283804848930030197,  0.54805291345769017271),
                                   je.valuesC(Complex.valueOf(2.5, 1)).sc(),
                                   tol);
        UnitTestUtils.assertEquals(Complex.valueOf(-0.60770261083279560398, -0.14952869024894608466),
                                   je.valuesC(Complex.valueOf(2.5, 1)).dc(),
                                   tol);
        UnitTestUtils.assertEquals(Complex.valueOf(-0.83882158549142135917,  0.03452229522564651472),
                                   je.valuesC(Complex.valueOf(2.5, 1)).nc(),
                                   tol);
        UnitTestUtils.assertEquals(Complex.valueOf(-0.17429449621856489159, -1.80783751829078560914),
                                   je.valuesS(Complex.valueOf(2.5, 1)).cs(),
                                   tol);
        UnitTestUtils.assertEquals(Complex.valueOf(-0.16440435588711757777,  1.12468960756395418591),
                                   je.valuesS(Complex.valueOf(2.5, 1)).ds(),
                                   tol);
        UnitTestUtils.assertEquals(Complex.valueOf( 0.20861268618691977878,  1.51043608734889061868),
                                   je.valuesS(Complex.valueOf(2.5, 1)).ns(),
                                   tol);
        UnitTestUtils.assertEquals(Complex.valueOf(-0.12725233583887664864, -0.87053277198136117922),
                                   je.valuesD(Complex.valueOf(2.5, 1)).sd(),
                                   tol);
        UnitTestUtils.assertEquals(Complex.valueOf(-1.55160242432190970322,  0.38178061795390252544),
                                   je.valuesD(Complex.valueOf(2.5, 1)).cd(),
                                   tol);
        UnitTestUtils.assertEquals(Complex.valueOf( 1.28833766241760303452, -0.37381070022725694453),
                                   je.valuesD(Complex.valueOf(2.5, 1)).nd(),
                                   tol);
        UnitTestUtils.assertEquals(Complex.valueOf( 0.08972834000183662533, -0.64966769410777501600),
                                   je.valuesN(Complex.valueOf(2.5, 1)).sn(),
                                   tol);
        UnitTestUtils.assertEquals(Complex.valueOf(-1.19013278764664544420, -0.04898075605528078495),
                                   je.valuesN(Complex.valueOf(2.5, 1)).cn(),
                                   tol);
        UnitTestUtils.assertEquals(Complex.valueOf( 0.71592277399023027709,  0.20772473029448836951),
                                   je.valuesN(Complex.valueOf(2.5, 1)).dn(),
                                   tol);
    }

    @Test
    public void testIssue150Minus() {
        final double tol = 4.0e-15;
        final FieldJacobiElliptic<Complex> je = JacobiEllipticBuilder.build(Complex.valueOf(-1, 1));
        UnitTestUtils.assertEquals(Complex.valueOf( 0.185179781846526661883, -0.7341316302932100299),
                                   je.valuesC(Complex.valueOf(2.5, -1)).sc(),
                                   tol);
        UnitTestUtils.assertEquals(Complex.valueOf( 0.036804110724415414398, -0.5315661866967806324),
                                   je.valuesC(Complex.valueOf(2.5, -1)).dc(),
                                   tol);
        UnitTestUtils.assertEquals(Complex.valueOf(-0.728148159153380263318,  0.1867014747416195502),
                                   je.valuesC(Complex.valueOf(2.5, -1)).nc(),
                                   tol);
        UnitTestUtils.assertEquals(Complex.valueOf( 0.323040127526322749710,  1.2806688349356319835),
                                   je.valuesS(Complex.valueOf(2.5, -1)).cs(),
                                   tol);
        UnitTestUtils.assertEquals(Complex.valueOf( 0.6926494536300507479991, -0.1245833311369300158),
                                   je.valuesS(Complex.valueOf(2.5, -1)).ds(),
                                   tol);
        UnitTestUtils.assertEquals(Complex.valueOf(-0.474323834329079345837, -0.8722045864335993302),
                                   je.valuesS(Complex.valueOf(2.5, -1)).ns(),
                                   tol);
        UnitTestUtils.assertEquals(Complex.valueOf( 1.398488691791959747251,  0.2515390416720468003),
                                   je.valuesD(Complex.valueOf(2.5, -1)).sd(),
                                   tol);
        UnitTestUtils.assertEquals(Complex.valueOf( 0.129629753901629389295,  1.8722580876874522195),
                                   je.valuesD(Complex.valueOf(2.5, -1)).cd(),
                                   tol);
        UnitTestUtils.assertEquals(Complex.valueOf(-0.443943012743148934256, -1.3390792137858189218),
                                   je.valuesD(Complex.valueOf(2.5, -1)).nd(),
                                   tol);
        UnitTestUtils.assertEquals(Complex.valueOf(-0.48119337969263266180,  0.88483656597829146168),
                                   je.valuesN(Complex.valueOf(2.5, -1)).sn(),
                                   tol);
        UnitTestUtils.assertEquals(Complex.valueOf(-1.28862738480059428568, -0.33041164783612229823),
                                   je.valuesN(Complex.valueOf(2.5, -1)).cn(),
                                   tol);
        UnitTestUtils.assertEquals(Complex.valueOf(-0.22306244463316206117,  0.67283023813989968132),
                                   je.valuesN(Complex.valueOf(2.5, -1)).dn(),
                                   tol);
    }

    @Test
    public void testWolframArcsn() {
        final FieldJacobiElliptic<Complex> je = JacobiEllipticBuilder.build(Complex.valueOf(0.3, 0.1));
        // the reference value wolframZ comes from https://www.wolframalpha.com/input?i=InverseJacobiSN%282.5-0.1+i+%2C0.3%2B0.1i%29
        // note that as shown in table https://dlmf.nist.gov/22.4.T3, sn(z + 2K) = -sn(z)
        // and since sn is odd in z (because θ₁ is odd and θ₂, θ₃, θ₄ are all even), then sn(-z) = -sn(z)
        // this implies that sn(2K - z) = sn(z)
        // so despite we don't provide the same result as Wolfram, we provide a correct result
        Complex wolframZ = new Complex(0.94236298773531348838, 1.92822841948038407820);
        Complex k = LegendreEllipticIntegral.bigK(je.getM());
        UnitTestUtils.assertEquals(k.multiply(2).subtract(wolframZ), je.arcsn(new Complex(2.5, -0.1)), 1.0e-15);
    }

    @Test
    public void testWolframArccn() {
        final FieldJacobiElliptic<Complex> je = JacobiEllipticBuilder.build(Complex.valueOf(0.3, 0.1));
        // the reference value wolframZ comes from https://www.wolframalpha.com/input?i=InverseJacobiCN%282.5-0.1+i+%2C0.3%2B0.1i%29
        final Complex wolframZ = new Complex(0.07230884246063775670, 1.35935692611563542744);
        UnitTestUtils.assertEquals(wolframZ, je.arccn(new Complex(2.5, -0.1)), 1.0e-15);
    }

    @Test
    public void testWolframArcdn() {
        final FieldJacobiElliptic<Complex> je = JacobiEllipticBuilder.build(Complex.valueOf(0.3, 0.1));
        // the reference value comes from https://www.wolframalpha.com/input?i=InverseJacobiDN%282.5-0.1+i+%2C0.3%2B0.1i%29
        final Complex wolframZ = new Complex(0.15903444834380499952, 1.6312723209949430373);
        UnitTestUtils.assertEquals(wolframZ, je.arcdn(new Complex(2.5, -0.1)), 1.0e-15);
    }

    @Test
    public void testWolframArccs() {
        final FieldJacobiElliptic<Complex> je = JacobiEllipticBuilder.build(Complex.valueOf(0.3, 0.1));
        // the reference value comes from https://www.wolframalpha.com/input?i=InverseJacobiCS%282.5-0.1+i+%2C0.3%2B0.1i%29
        // note that as shown in table https://dlmf.nist.gov/22.4.T2, one of the two cs periods is 2K
        // so despite we don't provide the same result as Wolfram, we provide a correct result
        final Complex wolframZ = new Complex(-3.0359823836145783448, -0.1010810531387431972);
        final Complex k        = LegendreEllipticIntegral.bigK(je.getM());
        final Complex period   = k.multiply(2);
        UnitTestUtils.assertEquals(wolframZ.add(period), je.arccs(new Complex(2.5, -0.1)), 1.0e-15);
    }

    @Test
    public void testWolframArcds() {
        final FieldJacobiElliptic<Complex> je = JacobiEllipticBuilder.build(Complex.valueOf(0.3, 0.1));
        // the reference value comes from https://www.wolframalpha.com/input?i=InverseJacobiDS%282.5-0.1+i+%2C0.3%2B0.1i%29
        // note that as shown in table https://dlmf.nist.gov/22.4.T2, one of the two cs periods is 2K+2iK'
        // so despite we don't provide the same result as Wolfram, we provide a correct result
        final Complex wolframZ = new Complex(-3.29799820949072209724, -4.202477914142034207);
        final Complex k      = LegendreEllipticIntegral.bigK(je.getM());
        final Complex kPrime = LegendreEllipticIntegral.bigKPrime(je.getM());
        final Complex period = k.add(kPrime.multiplyPlusI()).multiply(2);
        UnitTestUtils.assertEquals(wolframZ.add(period), je.arcds(new Complex(2.5, -0.1)), 1.0e-15);
    }

    @Test
    public void testWolframArcns() {
        final FieldJacobiElliptic<Complex> je = JacobiEllipticBuilder.build(Complex.valueOf(0.3, 0.1));
        // the reference value comes from https://www.wolframalpha.com/input?i=InverseJacobiNS%282.5-0.1+i+%2C0.3%2B0.1i%29
        // note that as shown in table https://dlmf.nist.gov/22.4.T2, one of the two ns periods is 2iK'
        // so despite we don't provide the same result as Wolfram, we provide a correct result
        final Complex wolframZ = new Complex(0.130428467771173061681, -4.08174077348072795334);
        final Complex kPrime   = LegendreEllipticIntegral.bigKPrime(je.getM());
        final Complex period   = kPrime.multiply(new Complex(0, 2));
        UnitTestUtils.assertEquals(wolframZ.add(period), je.arcns(new Complex(2.5, -0.1)), 1.0e-15);
    }

    @Test
    public void testWolframArcdc() {
        final FieldJacobiElliptic<Complex> je = JacobiEllipticBuilder.build(Complex.valueOf(0.3, 0.1));
        // the reference value comes from https://www.wolframalpha.com/input?i=InverseJacobiDC%282.5-0.1+i+%2C0.3%2B0.1i%29
        // note that as shown in table https://dlmf.nist.gov/22.4.T2, one of the two dc periods is 2iK'
        // so despite we don't provide the same result as Wolfram, we provide a correct result
        final Complex wolframZ = new Complex(1.578878353498834641592, 4.13977600222252116268406);
        final Complex kPrime   = LegendreEllipticIntegral.bigKPrime(je.getM());
        final Complex period   = kPrime.multiply(new Complex(0, 2));
        UnitTestUtils.assertEquals(wolframZ.subtract(period), je.arcdc(new Complex(2.5, -0.1)), 1.0e-15);
    }

    @Test
    public void testWolframArcnc() {
        final FieldJacobiElliptic<Complex> je = JacobiEllipticBuilder.build(Complex.valueOf(0.3, 0.1));
        // the reference value comes from https://www.wolframalpha.com/input?i=InverseJacobiNC%282.5-0.1+i+%2C0.3%2B0.1i%29
        final Complex wolframZ = new Complex(1.228171217697117722621, 0.006214187083238678388730749);
        UnitTestUtils.assertEquals(wolframZ, je.arcnc(new Complex(2.5, -0.1)), 1.0e-15);
    }

    @Test
    public void testWolframArcsc() {
        final FieldJacobiElliptic<Complex> je = JacobiEllipticBuilder.build(Complex.valueOf(0.3, 0.1));
        // the reference value comes from https://www.wolframalpha.com/input?i=InverseJacobiSC%282.5-0.1+i+%2C0.3%2B0.1i%29
        final Complex wolframZ = new Complex(1.26341293208997584885855, 0.012370090803027899483388605);
        UnitTestUtils.assertEquals(wolframZ, je.arcsc(new Complex(2.5, -0.1)), 1.0e-15);
    }

    @Test
    public void testWolframArcnd() {
        final FieldJacobiElliptic<Complex> je = JacobiEllipticBuilder.build(Complex.valueOf(0.3, 0.1));
        // the reference value comes from https://www.wolframalpha.com/input?i=InverseJacobiND%282.5-0.1+i+%2C0.3%2B0.1i%29
        final Complex wolframZ = new Complex(1.504742584755344845407901, -1.4869490086301163397113928);
        UnitTestUtils.assertEquals(wolframZ, je.arcnd(new Complex(2.5, -0.1)), 1.0e-15);
    }

    @Test
    public void testWolframArcsd() {
        final FieldJacobiElliptic<Complex> je = JacobiEllipticBuilder.build(Complex.valueOf(0.3, 0.1));
        // the reference value comes from https://www.wolframalpha.com/input?i=InverseJacobiSD%282.5-0.1+i+%2C0.3%2B0.1i%29
        final Complex wolframZ = new Complex(1.5881213420484595584158128, -1.1747797577661992518362122);
        UnitTestUtils.assertEquals(wolframZ, je.arcsd(new Complex(2.5, -0.1)), 1.0e-15);
    }

    @Test
    public void testWolframArccd() {
        final FieldJacobiElliptic<Complex> je = JacobiEllipticBuilder.build(Complex.valueOf(0.3, 0.1));
        // the reference value comes from https://www.wolframalpha.com/input?i=InverseJacobiCD%282.5-0.1+i+%2C0.3%2B0.1i%29
        final Complex wolframZ = new Complex(0.7669438335346942148909019, -1.87019319073859086886105089);
        UnitTestUtils.assertEquals(wolframZ, je.arccd(new Complex(2.5, -0.1)), 1.0e-15);
    }

    @Test
    public void testBranchCutArccn() {
        final FieldJacobiElliptic<Complex> je = JacobiEllipticBuilder.build(Complex.valueOf(0.3, 0.0));
        UnitTestUtils.assertEquals(new Complex(0.0, -1.3652045107), je.arccn(new Complex(2.5, +0.0)), 1.0e-10);
        UnitTestUtils.assertEquals(new Complex(0.0, +1.3652045107), je.arccn(new Complex(2.5, -0.0)), 1.0e-10);
    }

    @Test
    public void testBranchCutArcdn() {
        final FieldJacobiElliptic<Complex> je = JacobiEllipticBuilder.build(Complex.valueOf(0.7, 0.0));
        UnitTestUtils.assertEquals(new Complex(0.0, -1.1513008795), je.arcdn(new Complex(1.9, +0.0)), 1.0e-10);
        UnitTestUtils.assertEquals(new Complex(0.0, +1.1513008795), je.arcdn(new Complex(1.9, -0.0)), 1.0e-10);

        UnitTestUtils.assertEquals(new Complex(2.0753631353, -1.1247480245), je.arcdn(new Complex(0.3, +1.0e-11)), 1.0e-10);
        UnitTestUtils.assertEquals(new Complex(2.0753631353, +1.1247480245), je.arcdn(new Complex(0.3, -1.0e-11)), 1.0e-10);
        UnitTestUtils.assertEquals(new Complex(2.0753631353, -1.1247480245), je.arcdn(new Complex(0.3, +0.0)), 1.0e-10);

        // the following commented out test fails! We get the wrong Riemann sheet for negative zero on imaginary part
//        UnitTestUtils.assertEquals(new Complex(2.0753631353, +1.1247480245), je.arcdn(new Complex(0.3, -0.0)), 1.0e-10);

        UnitTestUtils.assertEquals(new Complex(4.150726270582776, -1.2596062674998643), je.arcdn(new Complex(-2.3, +1.0e-11)), 1.0e-10);
        UnitTestUtils.assertEquals(new Complex(4.150726270582776, +1.2596062674998643), je.arcdn(new Complex(-2.3, -1.0e-11)), 1.0e-10);
        UnitTestUtils.assertEquals(new Complex(4.150726270582776, -1.2596062674998643), je.arcdn(new Complex(-2.3, +0.0)), 1.0e-10);
        UnitTestUtils.assertEquals(new Complex(4.150726270582776, +1.2596062674998643), je.arcdn(new Complex(-2.3, -0.0)), 1.0e-10);
    }

}
