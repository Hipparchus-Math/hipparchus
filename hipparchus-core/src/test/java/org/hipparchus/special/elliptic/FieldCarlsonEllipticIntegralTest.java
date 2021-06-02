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

import org.hipparchus.complex.FieldComplex;
import org.hipparchus.util.Decimal64;

public class FieldCarlsonEllipticIntegralTest extends CarlsonEllipticIntegralAbstractTest<FieldComplex<Decimal64>> {

    protected FieldComplex<Decimal64> buildComplex(double realPart) {
        return new FieldComplex<>(new Decimal64(realPart));
    }

    protected FieldComplex<Decimal64> buildComplex(double realPart, double imaginaryPart) {
        return new FieldComplex<>(new Decimal64(realPart), new Decimal64(imaginaryPart));
    }

    protected FieldComplex<Decimal64> rF(FieldComplex<Decimal64> x, FieldComplex<Decimal64> y, FieldComplex<Decimal64> z) {
        return CarlsonEllipticIntegral.rF(x, y, z);
    }

    protected FieldComplex<Decimal64> rC(FieldComplex<Decimal64> x, FieldComplex<Decimal64> y) {
        return CarlsonEllipticIntegral.rC(x, y);
    }

    protected FieldComplex<Decimal64> rJ(FieldComplex<Decimal64> x, FieldComplex<Decimal64> y,
                                         FieldComplex<Decimal64> z, FieldComplex<Decimal64> p) {
        return CarlsonEllipticIntegral.rJ(x, y, z, p);
    }

    protected FieldComplex<Decimal64> rD(FieldComplex<Decimal64> x, FieldComplex<Decimal64> y, FieldComplex<Decimal64> z) {
        return CarlsonEllipticIntegral.rD(x, y, z);
    }

    protected FieldComplex<Decimal64> rG(FieldComplex<Decimal64> x, FieldComplex<Decimal64> y, FieldComplex<Decimal64> z) {
        return CarlsonEllipticIntegral.rG(x, y, z);
    }

}
