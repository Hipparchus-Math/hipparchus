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

import org.hipparchus.complex.Complex;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well1024a;
import org.hipparchus.util.FastMath;
import org.junit.Test;

public class CarlsonEllipticIntegralTest extends CarlsonEllipticIntegralAbstractTest<Complex> {

    protected Complex buildComplex(double realPart) {
        return new Complex(realPart);
    }

    protected Complex buildComplex(double realPart, double imaginaryPart) {
        return new Complex(realPart, imaginaryPart);
    }

    protected Complex rF(Complex x, Complex y, Complex z) {
        return CarlsonEllipticIntegral.rF(x, y, z);
    }

    protected Complex rC(Complex x, Complex y) {
        return CarlsonEllipticIntegral.rC(x, y);
    }

    protected Complex rJ(Complex x, Complex y, Complex z, Complex p) {
        return CarlsonEllipticIntegral.rJ(x, y, z, p);
    }

    protected Complex rD(Complex x, Complex y, Complex z) {
        return CarlsonEllipticIntegral.rD(x, y, z);
    }

    protected Complex rG(Complex x, Complex y, Complex z) {
        return CarlsonEllipticIntegral.rG(x, y, z);
    }

    @Test
    public void testTmp() {
        Complex x = new Complex(0.4375, -0.58);
        Complex y = new Complex(-0.33, -0.70);
        Complex z = new Complex(0.96, 0.25);
        for (double dxR = 0; dxR < 1; dxR += FastMath.scalb(1.0, -7)) {
            CarlsonEllipticIntegral.checkRG(x.add(dxR), y, z);
        }
//        RandomGenerator random = new Well1024a(0x123456l);
//        for (int i = 0; i < 20; ++i) {
//            Complex x = new Complex(2 * random.nextDouble() - 1, 2 * random.nextDouble() - 1);
//            Complex y = new Complex(2 * random.nextDouble() - 1, 2 * random.nextDouble() - 1);
//            Complex z = new Complex(2 * random.nextDouble() - 1, 2 * random.nextDouble() - 1);
//            CarlsonEllipticIntegral.checkRG(x, y, z);
//        }
    }
}
