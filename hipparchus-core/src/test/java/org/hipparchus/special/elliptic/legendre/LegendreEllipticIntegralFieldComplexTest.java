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
package org.hipparchus.special.elliptic.legendre;

import org.hipparchus.complex.FieldComplex;
import org.hipparchus.util.Decimal64;

public class LegendreEllipticIntegralFieldComplexTest extends LegendreEllipticIntegralAbstractComplexTest<FieldComplex<Decimal64>> {

    protected FieldComplex<Decimal64> buildComplex(double realPart) {
        return new FieldComplex<>(new Decimal64(realPart));
    }

    protected FieldComplex<Decimal64> buildComplex(double realPart, double imaginaryPart) {
        return new FieldComplex<>(new Decimal64(realPart), new Decimal64(imaginaryPart));
    }

    protected FieldComplex<Decimal64> K(FieldComplex<Decimal64> k) {
        return LegendreEllipticIntegral.bigK(k);
    }

    protected FieldComplex<Decimal64> Kprime(FieldComplex<Decimal64> k) {
        return LegendreEllipticIntegral.bigKPrime(k);
    }

    protected FieldComplex<Decimal64> F(FieldComplex<Decimal64> phi, FieldComplex<Decimal64> k) {
        return LegendreEllipticIntegral.bigF(phi, k);
    }

    protected FieldComplex<Decimal64> E(FieldComplex<Decimal64> k) {
        return LegendreEllipticIntegral.bigE(k);
    }

    protected FieldComplex<Decimal64> E(FieldComplex<Decimal64> phi, FieldComplex<Decimal64> k) {
        return LegendreEllipticIntegral.bigE(phi, k);
    }

    protected FieldComplex<Decimal64> D(FieldComplex<Decimal64> k) {
        return LegendreEllipticIntegral.bigD(k);
    }

    protected FieldComplex<Decimal64> D(FieldComplex<Decimal64> phi, FieldComplex<Decimal64> k) {
        return LegendreEllipticIntegral.bigD(phi, k);
    }

    protected FieldComplex<Decimal64> Pi(FieldComplex<Decimal64> alpha2, FieldComplex<Decimal64> k) {
        return LegendreEllipticIntegral.bigPi(alpha2, k);
    }

    protected FieldComplex<Decimal64> Pi(FieldComplex<Decimal64> phi, FieldComplex<Decimal64> alpha2, FieldComplex<Decimal64> k) {
        return LegendreEllipticIntegral.bigPi(phi, alpha2, k);
    }

}
