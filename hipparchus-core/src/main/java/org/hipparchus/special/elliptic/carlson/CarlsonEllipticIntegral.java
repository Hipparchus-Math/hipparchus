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

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.complex.Complex;
import org.hipparchus.complex.FieldComplex;
import org.hipparchus.util.FastMath;

/** Elliptic integrals in Carlson symmetric form.
 * <p>
 * This utility class computes the various symmetric elliptic
 * integrals defined as:
 * \[
 *   \left\{\begin{align}
 *   R_F(x,y,z)   &amp;= \frac{1}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{s(t)}\\
 *   R_J(x,y,z,p) &amp;= \frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{s(t)(t+p)}\\
 *   R_G(x,y,z)   &amp;= \frac{1}{4}\int_{0}^{\infty}\frac{1}{s(t)}
                     \left(\frac{x}{t+x}+\frac{y}{t+y}+\frac{z}{t+z}\right)t\mathrm{d}t\\
 *   R_D(x,y,z)   &amp;= R_J(x,y,z,z)\\
 *   R_C(x,y)     &amp;= R_F(x,y,y)
 *   \end{align}\right.
 * \]
 * </p>
 * <p>
 * where
 * \[
 *   s(t) = \sqrt{t+x}\sqrt{t+y}\sqrt{t+z}
 * \]
 * </p>
 * <p>
 * The algorithms used are based on the duplication method as described in
 * B. C. Carlson 1995 paper "Numerical computation of real or complex
 * elliptic integrals", with the improvements described in the appendix
 * of B. C. Carlson and James FitzSimons 2000 paper "Reduction theorems
 * for elliptic integrands with the square root of two quadratic factors".
 * They are also described in <a href="https://dlmf.nist.gov/19.36#i">section 19.36(i)</a>
 * of Digital Library of Mathematical Functions.
 * </p>
 * <p>
 * <em>
 * Beware that when computing elliptic integrals in the complex plane,
 * many issues arise due to branch cuts. See the
 * <a href="https://www.hipparchus.org/hipparchus-core/special.html#Elliptic_functions_and_integrals">user guide</a>
 * for a thorough explanation.
 * </em>
 * </p>
 * @since 2.0
 */
public class CarlsonEllipticIntegral {

    /** Private constructor for a utility class.
     */
    private CarlsonEllipticIntegral() {
    }

    /** Compute Carlson elliptic integral R<sub>C</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>C</sub>is defined as
     * \[
     *   R_C(x,y,z)=R_F(x,y,y)=\frac{1}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}(t+y)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @return Carlson elliptic integral R<sub>C</sub>
     */
    public static double rC(final double x, final double y) {
        if (y < 0) {
            // y is on the branch cut, we must use a transformation to get the Cauchy principal value
            // see equation 2.14 in Carlson[1995]
            final double xMy = x - y;
            return FastMath.sqrt(x / xMy) * new RcRealDuplication(xMy, -y).integral();
        } else {
            return new RcRealDuplication(x, y).integral();
        }
    }

    /** Compute Carlson elliptic integral R<sub>C</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>C</sub>is defined as
     * \[
     *   R_C(x,y,z)=R_F(x,y,y)=\frac{1}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}(t+y)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param <T> type of the field elements
     * @return Carlson elliptic integral R<sub>C</sub>
     */
    public static <T extends CalculusFieldElement<T>> T rC(final T x, final T y) {
        if (y.getReal() < 0) {
            // y is on the branch cut, we must use a transformation to get the Cauchy principal value
            // see equation 2.14 in Carlson[1995]
            final T xMy = x.subtract(y);
            return FastMath.sqrt(x.divide(xMy)).multiply(new RcFieldDuplication<>(xMy, y.negate()).integral());
        } else {
            return new RcFieldDuplication<>(x, y).integral();
        }
    }

    /** Compute Carlson elliptic integral R<sub>C</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>C</sub>is defined as
     * \[
     *   R_C(x,y,z)=R_F(x,y,y)=\frac{1}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}(t+y)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @return Carlson elliptic integral R<sub>C</sub>
     */
    public static Complex rC(final Complex x, final Complex y) {
        if (y.getImaginaryPart() == 0 && y.getRealPart() < 0) {
            // y is on the branch cut, we must use a transformation to get the Cauchy principal value
            // see equation 2.14 in Carlson[1995]
            final Complex xMy = x.subtract(y);
            return FastMath.sqrt(x.divide(xMy)).multiply(new RcFieldDuplication<>(xMy, y.negate()).integral());
        } else {
            return new RcFieldDuplication<>(x, y).integral();
        }
    }

    /** Compute Carlson elliptic integral R<sub>C</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>C</sub>is defined as
     * \[
     *   R_C(x,y,z)=R_F(x,y,y)=\frac{1}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}(t+y)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param <T> type of the field elements
     * @return Carlson elliptic integral R<sub>C</sub>
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> rC(final FieldComplex<T> x, final FieldComplex<T> y) {
        if (y.getImaginaryPart().isZero() && y.getRealPart().getReal() < 0) {
            // y is on the branch cut, we must use a transformation to get the Cauchy principal value
            // see equation 2.14 in Carlson[1995]
            final FieldComplex<T> xMy = x.subtract(y);
            return FastMath.sqrt(x.divide(xMy)).multiply(new RcFieldDuplication<>(xMy, y.negate()).integral());
        } else {
            return new RcFieldDuplication<>(x, y).integral();
        }
    }

    /** Compute Carlson elliptic integral R<sub>F</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>F</sub> is defined as
     * \[
     *   R_F(x,y,z)=\frac{1}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @return Carlson elliptic integral R<sub>F</sub>
     */
    public static double rF(final double x, final double y, final double z) {
        return new RfRealDuplication(x, y, z).integral();
    }

    /** Compute Carlson elliptic integral R<sub>F</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>F</sub> is defined as
     * \[
     *   R_F(x,y,z)=\frac{1}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param <T> type of the field elements
     * @return Carlson elliptic integral R<sub>F</sub>
     */
    public static <T extends CalculusFieldElement<T>> T rF(final T x, final T y, final T z) {
        return new RfFieldDuplication<>(x, y, z).integral();
    }

    /** Compute Carlson elliptic integral R<sub>F</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>F</sub> is defined as
     * \[
     *   R_F(x,y,z)=\frac{1}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @return Carlson elliptic integral R<sub>F</sub>
     */
    public static Complex rF(final Complex x, final Complex y, final Complex z) {
        return new RfFieldDuplication<>(x, y, z).integral();
    }

    /** Compute Carlson elliptic integral R<sub>F</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>F</sub> is defined as
     * \[
     *   R_F(x,y,z)=\frac{1}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param <T> type of the field elements
     * @return Carlson elliptic integral R<sub>F</sub>
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> rF(final FieldComplex<T> x, final FieldComplex<T> y, final FieldComplex<T> z) {
        return new RfFieldDuplication<>(x, y, z).integral();
    }

    /** Compute Carlson elliptic integral R<sub>J</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>J</sub> is defined as
     * \[
     *   R_J(x,y,z,p)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+p)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param p fourth <em>not</em> symmetric variable of the integral
     * @return Carlson elliptic integral R<sub>J</sub>
     */
    public static double rJ(final double x, final double y, final double z, final double p) {
        final double delta = (p - x) * (p - y) * (p - z);
        return rJ(x, y, z, p, delta);
    }

    /** Compute Carlson elliptic integral R<sub>J</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>J</sub> is defined as
     * \[
     *   R_J(x,y,z,p)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+p)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param p fourth <em>not</em> symmetric variable of the integral
     * @param delta precomputed value of (p-x)(p-y)(p-z)
     * @return Carlson elliptic integral R<sub>J</sub>
     */
    public static double rJ(final double x, final double y, final double z, final double p, final double delta) {
        return new RjRealDuplication(x, y, z, p, delta).integral();
    }

    /** Compute Carlson elliptic integral R<sub>J</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>J</sub> is defined as
     * \[
     *   R_J(x,y,z,p)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+p)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param p fourth <em>not</em> symmetric variable of the integral
     * @param <T> type of the field elements
     * @return Carlson elliptic integral R<sub>J</sub>
     */
    public static <T extends CalculusFieldElement<T>> T rJ(final T x, final T y, final T z, final T p) {
        final T delta = p.subtract(x).multiply(p.subtract(y)).multiply(p.subtract(z));
        return new RjFieldDuplication<>(x, y, z, p, delta). integral();
    }

    /** Compute Carlson elliptic integral R<sub>J</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>J</sub> is defined as
     * \[
     *   R_J(x,y,z,p)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+p)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param p fourth <em>not</em> symmetric variable of the integral
     * @param delta precomputed value of (p-x)(p-y)(p-z)
     * @param <T> type of the field elements
     * @return Carlson elliptic integral R<sub>J</sub>
     */
    public static <T extends CalculusFieldElement<T>> T rJ(final T x, final T y, final T z, final T p, final T delta) {
        return new RjFieldDuplication<>(x, y, z, p, delta).integral();
    }

    /** Compute Carlson elliptic integral R<sub>J</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>J</sub> is defined as
     * \[
     *   R_J(x,y,z,p)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+p)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param p fourth <em>not</em> symmetric variable of the integral
     * @return Carlson elliptic integral R<sub>J</sub>
     */
    public static Complex rJ(final Complex x, final Complex y, final Complex z, final Complex p) {
        final Complex delta = p.subtract(x).multiply(p.subtract(y)).multiply(p.subtract(z));
        return new RjFieldDuplication<>(x, y, z, p, delta).integral();
    }

    /** Compute Carlson elliptic integral R<sub>J</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>J</sub> is defined as
     * \[
     *   R_J(x,y,z,p)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+p)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param p fourth <em>not</em> symmetric variable of the integral
     * @param delta precomputed value of (p-x)(p-y)(p-z)
     * @return Carlson elliptic integral R<sub>J</sub>
     */
    public static Complex rJ(final Complex x, final Complex y, final Complex z, final Complex p, final Complex delta) {
        return new RjFieldDuplication<>(x, y, z, p, delta).integral();
    }

    /** Compute Carlson elliptic integral R<sub>J</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>J</sub> is defined as
     * \[
     *   R_J(x,y,z,p)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+p)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param p fourth <em>not</em> symmetric variable of the integral
     * @param <T> type of the field elements
     * @return Carlson elliptic integral R<sub>J</sub>
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> rJ(final FieldComplex<T> x, final FieldComplex<T> y,
                                                                         final FieldComplex<T> z, final FieldComplex<T> p) {
        final FieldComplex<T> delta = p.subtract(x).multiply(p.subtract(y)).multiply(p.subtract(z));
        return new RjFieldDuplication<>(x, y, z, p, delta).integral();
    }

    /** Compute Carlson elliptic integral R<sub>J</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>J</sub> is defined as
     * \[
     *   R_J(x,y,z,p)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+p)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param p fourth <em>not</em> symmetric variable of the integral
     * @param delta precomputed value of (p-x)(p-y)(p-z)
     * @param <T> type of the field elements
     * @return Carlson elliptic integral R<sub>J</sub>
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> rJ(final FieldComplex<T> x, final FieldComplex<T> y,
                                                                         final FieldComplex<T> z, final FieldComplex<T> p,
                                                                         final FieldComplex<T> delta) {
        return new RjFieldDuplication<>(x, y, z, p, delta).integral();
    }

    /** Compute Carlson elliptic integral R<sub>D</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>D</sub> is defined as
     * \[
     *   R_D(x,y,z)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+z)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @return Carlson elliptic integral R<sub>D</sub>
     */
    public static double rD(final double x, final double y, final double z) {
        return new RdRealDuplication(x, y, z).integral();
    }

    /** Compute Carlson elliptic integral R<sub>D</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>D</sub> is defined as
     * \[
     *   R_D(x,y,z)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+z)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param <T> type of the field elements
     * @return Carlson elliptic integral R<sub>D</sub>
     */
    public static <T extends CalculusFieldElement<T>> T rD(final T x, final T y, final T z) {
        return new RdFieldDuplication<>(x, y, z).integral();
    }

    /** Compute Carlson elliptic integral R<sub>D</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>D</sub> is defined as
     * \[
     *   R_D(x,y,z)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+z)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @return Carlson elliptic integral R<sub>D</sub>
     */
    public static Complex rD(final Complex x, final Complex y, final Complex z) {
        return new RdFieldDuplication<>(x, y, z).integral();
    }

    /** Compute Carlson elliptic integral R<sub>D</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>D</sub> is defined as
     * \[
     *   R_D(x,y,z)=\frac{3}{2}\int_{0}^{\infty}\frac{\mathrm{d}t}{\sqrt{t+x}\sqrt{t+y}\sqrt{t+z}(t+z)}
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param <T> type of the field elements
     * @return Carlson elliptic integral R<sub>D</sub>
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> rD(final FieldComplex<T> x, final FieldComplex<T> y,
                                                                         final FieldComplex<T> z) {
        return new RdFieldDuplication<>(x, y, z).integral();
    }

    /** Compute Carlson elliptic integral R<sub>G</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>G</sub>is defined as
     * \[
     *   R_{G}(x,y,z)=\frac{1}{4}\int_{0}^{\infty}\frac{1}{s(t)}
     *                \left(\frac{x}{t+x}+\frac{y}{t+y}+\frac{z}{t+z}\right)t\mathrm{d}t
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z second symmetric variable of the integral
     * @return Carlson elliptic integral R<sub>G</sub>
     */
    public static double rG(final double x, final double y, final double z) {
        return generalComputeRg(x, y, z);
    }

    /** Compute Carlson elliptic integral R<sub>G</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>G</sub>is defined as
     * \[
     *   R_{G}(x,y,z)=\frac{1}{4}\int_{0}^{\infty}\frac{1}{s(t)}
     *                \left(\frac{x}{t+x}+\frac{y}{t+y}+\frac{z}{t+z}\right)t\mathrm{d}t
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z second symmetric variable of the integral
     * @param <T> type of the field elements
     * @return Carlson elliptic integral R<sub>G</sub>
     */
    public static <T extends CalculusFieldElement<T>> T rG(final T x, final T y, final T z) {
        return generalComputeRg(x, y, z);
    }

    /** Compute Carlson elliptic integral R<sub>G</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>G</sub>is defined as
     * \[
     *   R_{G}(x,y,z)=\frac{1}{4}\int_{0}^{\infty}\frac{1}{s(t)}
     *                \left(\frac{x}{t+x}+\frac{y}{t+y}+\frac{z}{t+z}\right)t\mathrm{d}t
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z second symmetric variable of the integral
     * @return Carlson elliptic integral R<sub>G</sub>
     */
    public static Complex rG(final Complex x, final Complex y, final Complex z) {
        return generalComputeRg(x, y, z);
    }

    /** Compute Carlson elliptic integral R<sub>G</sub>.
     * <p>
     * The Carlson elliptic integral R<sub>G</sub>is defined as
     * \[
     *   R_{G}(x,y,z)=\frac{1}{4}\int_{0}^{\infty}\frac{1}{s(t)}
     *                \left(\frac{x}{t+x}+\frac{y}{t+y}+\frac{z}{t+z}\right)t\mathrm{d}t
     * \]
     * </p>
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z second symmetric variable of the integral
     * @param <T> type of the field elements
     * @return Carlson elliptic integral R<sub>G</sub>
     */
    public static <T extends CalculusFieldElement<T>> FieldComplex<T> rG(final FieldComplex<T> x,
                                                                         final FieldComplex<T> y,
                                                                         final FieldComplex<T> z) {
        return generalComputeRg(x, y, z);
    }

    /** Compute Carlson elliptic integral R<sub>G</sub> in the general case.
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @return Carlson elliptic integral R<sub>G</sub>
     */
    private static double generalComputeRg(final double x, final double y, final double z) {
        // permute parameters if needed to avoid cancellations
        if (x <= y) {
            if (y <= z) {
                // x ≤ y ≤ z
                return permutedComputeRg(x, z, y);
            } else if (x <= z) {
                // x ≤ z < y
                return permutedComputeRg(x, y, z);
            } else {
                // z < x ≤ y
                return permutedComputeRg(z, y, x);
            }
        } else if (x <= z) {
            // y < x ≤ z
            return permutedComputeRg(y, z, x);
        } else if (y <= z) {
            // y ≤ z < x
            return permutedComputeRg(y, x, z);
        } else {
            // z < y < x
            return permutedComputeRg(z, x, y);
        }
    }

    /** Compute Carlson elliptic integral R<sub>G</sub> in the general case.
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param <T> type of the field elements (really {@link Complex} or {@link FieldComplex})
     * @return Carlson elliptic integral R<sub>G</sub>
     */
    private static <T extends CalculusFieldElement<T>> T generalComputeRg(final T x, final T y, final T z) {
        // permute parameters if needed to avoid cancellations
        final double xR = x.getReal();
        final double yR = y.getReal();
        final double zR = z.getReal();
        if (xR <= yR) {
            if (yR <= zR) {
                // x ≤ y ≤ z
                return permutedComputeRg(x, z, y);
            } else if (xR <= zR) {
                // x ≤ z < y
                return permutedComputeRg(x, y, z);
            } else {
                // z < x ≤ y
                return permutedComputeRg(z, y, x);
            }
        } else if (xR <= zR) {
            // y < x ≤ z
            return permutedComputeRg(y, z, x);
        } else if (yR <= zR) {
            // y ≤ z < x
            return permutedComputeRg(y, x, z);
        } else {
            // z < y < x
            return permutedComputeRg(z, x, y);
        }
    }

    /** Compute Carlson elliptic integral R<sub>G</sub> with already permuted variables to avoid cancellations.
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @return Carlson elliptic integral R<sub>G</sub>
     */
    private static double permutedComputeRg(final double x, final double y, final double z) {
        // permute parameters if needed to avoid divisions by zero
        if (z == 0) {
            return x == 0 ? safeComputeRg(z, x, y) : safeComputeRg(y, z, x);
        } else {
            return safeComputeRg(x, y, z);
        }
    }

    /** Compute Carlson elliptic integral R<sub>G</sub> with already permuted variables to avoid cancellations.
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param <T> type of the field elements (really {@link Complex} or {@link FieldComplex})
     * @return Carlson elliptic integral R<sub>G</sub>
     */
    private static <T extends CalculusFieldElement<T>> T permutedComputeRg(final T x, final T y, final T z) {
        // permute parameters if needed to avoid divisions by zero
        if (z.isZero()) {
            return x.isZero() ? safeComputeRg(z, x, y) : safeComputeRg(y, z, x);
        } else {
            return safeComputeRg(x, y, z);
        }
    }

    /** Compute Carlson elliptic integral R<sub>G</sub> with non-zero third variable.
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @see <a href="https://dlmf.nist.gov/19.21#E10">Digital Library of Mathematical Functions, equation 19.21.10</a>
     * @return Carlson elliptic integral R<sub>G</sub>
     */
    private static double safeComputeRg(final double x, final double y, final double z) {

        // contribution of the R_F integral
        final double termF = new RfRealDuplication(x, y, z).integral() * z;

        // contribution of the R_D integral
        final double termD = (x - z) * (y - z) * new RdRealDuplication(x, y, z).integral() / 3;

        // contribution of the square roots
        final double termS = FastMath.sqrt(x * y / z);

        // equation 19.21.10
        return (termF - termD + termS) * 0.5;

    }

    /** Compute Carlson elliptic integral R<sub>G</sub> with non-zero third variable.
     * @param x first symmetric variable of the integral
     * @param y second symmetric variable of the integral
     * @param z third symmetric variable of the integral
     * @param <T> type of the field elements (really {@link Complex} or {@link FieldComplex})
     * @see <a href="https://dlmf.nist.gov/19.21#E10">Digital Library of Mathematical Functions, equation 19.21.10</a>
     * @return Carlson elliptic integral R<sub>G</sub>
     */
    private static <T extends CalculusFieldElement<T>> T safeComputeRg(final T x, final T y, final T z) {

        // contribution of the R_F integral
        final T termF = new RfFieldDuplication<>(x, y, z).integral().multiply(z);

        // contribution of the R_D integral
        final T termD = x.subtract(z).multiply(y.subtract(z)).multiply(new RdFieldDuplication<>(x, y, z).integral()).divide(3);

        // contribution of the square roots
        // BEWARE: this term MUST be computed as √x√y/√z with all square roots selected with positive real part
        // and NOT as √(xy/z), otherwise sign errors may occur
        final T termS = x.sqrt().multiply(y.sqrt()).divide(z.sqrt());

        // equation 19.21.10
        return termF.subtract(termD).add(termS).multiply(0.5);

    }

}
