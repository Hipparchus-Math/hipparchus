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

import org.hipparchus.complex.Complex;

public class CarlsonEllipticIntegralComplexTest extends CarlsonEllipticIntegralAbstractComplexTest<Complex> {

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

}
