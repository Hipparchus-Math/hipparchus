/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */

package org.hipparchus.complex;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

/**
 */
public class ComplexUtilsTest {

    private final double inf = Double.POSITIVE_INFINITY;
    private final double negInf = Double.NEGATIVE_INFINITY;
    private final double nan = Double.NaN;
    private final double pi = FastMath.PI;

    private final Complex negInfInf = new Complex(negInf, inf);
    private final Complex infNegInf = new Complex(inf, negInf);
    private final Complex infInf = new Complex(inf, inf);
    private final Complex negInfNegInf = new Complex(negInf, negInf);
    private final Complex infNaN = new Complex(inf, nan);

    @Test
    public void testPolar2Complex() {
        UnitTestUtils.assertEquals(Complex.ONE,
                ComplexUtils.polar2Complex(1, 0), 10e-12);
        UnitTestUtils.assertEquals(Complex.ZERO,
                ComplexUtils.polar2Complex(0, 1), 10e-12);
        UnitTestUtils.assertEquals(Complex.ZERO,
                ComplexUtils.polar2Complex(0, -1), 10e-12);
        UnitTestUtils.assertEquals(Complex.I,
                ComplexUtils.polar2Complex(1, pi/2), 10e-12);
        UnitTestUtils.assertEquals(Complex.I.negate(),
                ComplexUtils.polar2Complex(1, -pi/2), 10e-12);
        double r = 0;
        for (int i = 0; i < 5; i++) {
          r += i;
          double theta = 0;
          for (int j =0; j < 20; j++) {
              theta += pi / 6;
              UnitTestUtils.assertEquals(altPolar(r, theta),
                      ComplexUtils.polar2Complex(r, theta), 10e-12);
          }
          theta = -2 * pi;
          for (int j =0; j < 20; j++) {
              theta -= pi / 6;
              UnitTestUtils.assertEquals(altPolar(r, theta),
                      ComplexUtils.polar2Complex(r, theta), 10e-12);
          }
        }
    }

    protected Complex altPolar(double r, double theta) {
        return Complex.I.multiply(new Complex(theta, 0)).exp().multiply(new Complex(r, 0));
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testPolar2ComplexIllegalModulus() {
        ComplexUtils.polar2Complex(-1, 0);
    }

    @Test
    public void testPolar2ComplexNaN() {
        UnitTestUtils.assertSame(Complex.NaN, ComplexUtils.polar2Complex(nan, 1));
        UnitTestUtils.assertSame(Complex.NaN, ComplexUtils.polar2Complex(1, nan));
        UnitTestUtils.assertSame(Complex.NaN,
                ComplexUtils.polar2Complex(nan, nan));
    }

    @Test
    public void testPolar2ComplexInf() {
        UnitTestUtils.assertSame(Complex.NaN, ComplexUtils.polar2Complex(1, inf));
        UnitTestUtils.assertSame(Complex.NaN,
                ComplexUtils.polar2Complex(1, negInf));
        UnitTestUtils.assertSame(Complex.NaN, ComplexUtils.polar2Complex(inf, inf));
        UnitTestUtils.assertSame(Complex.NaN,
                ComplexUtils.polar2Complex(inf, negInf));
        UnitTestUtils.assertSame(infInf, ComplexUtils.polar2Complex(inf, pi/4));
        UnitTestUtils.assertSame(infNaN, ComplexUtils.polar2Complex(inf, 0));
        UnitTestUtils.assertSame(infNegInf, ComplexUtils.polar2Complex(inf, -pi/4));
        UnitTestUtils.assertSame(negInfInf, ComplexUtils.polar2Complex(inf, 3*pi/4));
        UnitTestUtils.assertSame(negInfNegInf, ComplexUtils.polar2Complex(inf, 5*pi/4));
    }

    @Test
    public void testConvertToComplex() {
        final double[] real = new double[] { negInf, -123.45, 0, 1, 234.56, pi, inf };
        final Complex[] complex = ComplexUtils.convertToComplex(real);

        for (int i = 0; i < real.length; i++) {
            Assert.assertEquals(real[i], complex[i].getReal(), 0d);
        }
    }
}
