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
import org.hipparchus.analysis.integration.IterativeLegendreGaussIntegrator;
import org.hipparchus.complex.Complex;
import org.hipparchus.complex.ComplexUnivariateIntegrator;

public class LegendreEllipticIntegralComplexTest extends LegendreEllipticIntegralAbstractComplexTest<Complex> {

    private ComplexUnivariateIntegrator integrator() {
        return new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                    1.0e-6,
                                                                                    1.0e-6));
    }

    protected Complex buildComplex(double realPart) {
        return new Complex(realPart);
    }

    protected Complex buildComplex(double realPart, double imaginaryPart) {
        return new Complex(realPart, imaginaryPart);
    }

    protected Complex K(Complex m) {
        return LegendreEllipticIntegral.bigK(m);
    }

    protected Complex Kprime(Complex m) {
        return LegendreEllipticIntegral.bigKPrime(m);
    }

    protected Complex F(Complex phi, Complex m) {
        return LegendreEllipticIntegral.bigF(phi, m);
    }

    protected Complex integratedF(Complex phi, Complex m) {
        return LegendreEllipticIntegral.bigF(phi, m, integrator(), 100000);
    }

    protected Complex E(Complex m) {
        return LegendreEllipticIntegral.bigE(m);
    }

    protected Complex E(Complex phi, Complex m) {
        return LegendreEllipticIntegral.bigE(phi, m);
    }

    protected Complex integratedE(Complex phi, Complex m) {
        return LegendreEllipticIntegral.bigE(phi, m, integrator(), 100000);
    }

    protected Complex D(Complex m) {
        return LegendreEllipticIntegral.bigD(m);
    }

    protected Complex D(Complex phi, Complex m) {
        return LegendreEllipticIntegral.bigD(phi, m);
    }

    protected Complex Pi(Complex n, Complex m) {
        return LegendreEllipticIntegral.bigPi(n, m);
    }

    protected Complex Pi(Complex n, Complex phi, Complex m) {
        return LegendreEllipticIntegral.bigPi(n, phi, m);
    }

    protected Complex integratedPi(Complex n, Complex phi, Complex m) {
        return LegendreEllipticIntegral.bigPi(n, phi, m, integrator(), 100000);
    }

    protected Complex integrate(int maxEval, CalculusFieldUnivariateFunction<Complex> f, Complex start, Complex end) {
        return integrator().integrate(maxEval, f, start, end);
    }

    protected Complex integrate(int maxEval, CalculusFieldUnivariateFunction<Complex> f, Complex start, Complex middle, Complex end) {
        return integrator().integrate(maxEval, f, start, middle, end);
    }

}
