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
package org.hipparchus.special.elliptic.jacobi;

import java.io.IOException;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.complex.Complex;
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

}
