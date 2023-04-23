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

import org.hipparchus.complex.FieldComplex;
import org.hipparchus.util.Binary64;

public class CarlsonEllipticIntegralFieldComplexTest extends CarlsonEllipticIntegralAbstractComplexTest<FieldComplex<Binary64>> {

    protected FieldComplex<Binary64> buildComplex(double realPart) {
        return new FieldComplex<>(new Binary64(realPart));
    }

    protected FieldComplex<Binary64> buildComplex(double realPart, double imaginaryPart) {
        return new FieldComplex<>(new Binary64(realPart), new Binary64(imaginaryPart));
    }

    protected FieldComplex<Binary64> rF(FieldComplex<Binary64> x, FieldComplex<Binary64> y, FieldComplex<Binary64> z) {
        return CarlsonEllipticIntegral.rF(x, y, z);
    }

    protected FieldComplex<Binary64> rC(FieldComplex<Binary64> x, FieldComplex<Binary64> y) {
        return CarlsonEllipticIntegral.rC(x, y);
    }

    protected FieldComplex<Binary64> rJ(FieldComplex<Binary64> x, FieldComplex<Binary64> y,
                                         FieldComplex<Binary64> z, FieldComplex<Binary64> p) {
        return CarlsonEllipticIntegral.rJ(x, y, z, p);
    }

    protected FieldComplex<Binary64> rD(FieldComplex<Binary64> x, FieldComplex<Binary64> y, FieldComplex<Binary64> z) {
        return CarlsonEllipticIntegral.rD(x, y, z);
    }

    protected FieldComplex<Binary64> rG(FieldComplex<Binary64> x, FieldComplex<Binary64> y, FieldComplex<Binary64> z) {
        return CarlsonEllipticIntegral.rG(x, y, z);
    }

}
