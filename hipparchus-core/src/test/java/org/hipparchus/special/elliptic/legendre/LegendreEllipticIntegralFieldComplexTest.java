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
package org.hipparchus.special.elliptic.legendre;

import org.hipparchus.analysis.CalculusFieldUnivariateFunction;
import org.hipparchus.analysis.integration.IterativeLegendreFieldGaussIntegrator;
import org.hipparchus.complex.FieldComplex;
import org.hipparchus.complex.FieldComplexUnivariateIntegrator;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.Binary64Field;

public class LegendreEllipticIntegralFieldComplexTest extends LegendreEllipticIntegralAbstractComplexTest<FieldComplex<Binary64>> {

    private FieldComplexUnivariateIntegrator<Binary64> integrator() {
        return new FieldComplexUnivariateIntegrator<>(new IterativeLegendreFieldGaussIntegrator<>(Binary64Field.getInstance(),
                                                                                                  24,
                                                                                                  1.0e-6,
                                                                                                  1.0e-6));
    }

    protected FieldComplex<Binary64> buildComplex(double realPart) {
        return new FieldComplex<>(new Binary64(realPart));
    }

    protected FieldComplex<Binary64> buildComplex(double realPart, double imaginaryPart) {
        return new FieldComplex<>(new Binary64(realPart), new Binary64(imaginaryPart));
    }

    protected FieldComplex<Binary64> K(FieldComplex<Binary64> m) {
        return LegendreEllipticIntegral.bigK(m);
    }

    protected FieldComplex<Binary64> Kprime(FieldComplex<Binary64> m) {
        return LegendreEllipticIntegral.bigKPrime(m);
    }

    protected FieldComplex<Binary64> F(FieldComplex<Binary64> phi, FieldComplex<Binary64> m) {
        return LegendreEllipticIntegral.bigF(phi, m);
    }

    protected FieldComplex<Binary64> integratedF(FieldComplex<Binary64> phi, FieldComplex<Binary64> m) {
        return LegendreEllipticIntegral.bigF(phi, m, integrator(), 100000);
    }

    protected FieldComplex<Binary64> E(FieldComplex<Binary64> m) {
        return LegendreEllipticIntegral.bigE(m);
    }

    protected FieldComplex<Binary64> E(FieldComplex<Binary64> phi, FieldComplex<Binary64> m) {
        return LegendreEllipticIntegral.bigE(phi, m);
    }

    protected FieldComplex<Binary64> integratedE(FieldComplex<Binary64> phi, FieldComplex<Binary64> m) {
        return LegendreEllipticIntegral.bigE(phi, m, integrator(), 100000);
    }

    protected FieldComplex<Binary64> D(FieldComplex<Binary64> m) {
        return LegendreEllipticIntegral.bigD(m);
    }

    protected FieldComplex<Binary64> D(FieldComplex<Binary64> phi, FieldComplex<Binary64> m) {
        return LegendreEllipticIntegral.bigD(phi, m);
    }

    protected FieldComplex<Binary64> Pi(FieldComplex<Binary64> n, FieldComplex<Binary64> m) {
        return LegendreEllipticIntegral.bigPi(n, m);
    }

    protected FieldComplex<Binary64> Pi(FieldComplex<Binary64> n, FieldComplex<Binary64> phi, FieldComplex<Binary64> m) {
        return LegendreEllipticIntegral.bigPi(n, phi, m);
    }

    protected FieldComplex<Binary64> integratedPi(FieldComplex<Binary64> n, FieldComplex<Binary64> phi, FieldComplex<Binary64> m) {
        return LegendreEllipticIntegral.bigPi(n, phi, m, integrator(), 100000);
    }

    protected FieldComplex<Binary64> integrate(int maxEval, CalculusFieldUnivariateFunction<FieldComplex<Binary64>> f,
                                                FieldComplex<Binary64> start, FieldComplex<Binary64> end) {
        return integrator().integrate(maxEval, f, start, end);
    }

    @SuppressWarnings("unchecked")
    protected FieldComplex<Binary64> integrate(int maxEval, CalculusFieldUnivariateFunction<FieldComplex<Binary64>> f,
                                                FieldComplex<Binary64> start, FieldComplex<Binary64> middle, FieldComplex<Binary64> end) {
        return integrator().integrate(maxEval, f, start, middle, end);
    }

}
