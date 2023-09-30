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
package org.hipparchus.analysis.differentiation;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.FieldSinCos;
import org.hipparchus.util.FieldSinhCosh;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.MathUtils;
import org.hipparchus.util.SinCos;
import org.hipparchus.util.SinhCosh;

/** Class representing both the value and the differentials of a function.
 * <p>This class is a stripped-down version of {@link DerivativeStructure}
 * with only one {@link DerivativeStructure#getFreeParameters() free parameter}
 * and {@link DerivativeStructure#getOrder() derivation order} also limited to two.
 * It should have less overhead than {@link DerivativeStructure} in its domain.</p>
 * <p>This class is an implementation of Rall's numbers. Rall's numbers are an
 * extension to the real numbers used throughout mathematical expressions; they hold
 * the derivative together with the value of a function.</p>
 * <p>{@link UnivariateDerivative2} instances can be used directly thanks to
 * the arithmetic operators to the mathematical functions provided as
 * methods by this class (+, -, *, /, %, sin, cos ...).</p>
 * <p>Implementing complex expressions by hand using these classes is
 * a tedious and error-prone task but has the advantage of having no limitation
 * on the derivation order despite not requiring users to compute the derivatives by
 * themselves.</p>
 * <p>Instances of this class are guaranteed to be immutable.</p>
 * @see DerivativeStructure
 * @see UnivariateDerivative2
 * @see Gradient
 * @see FieldDerivativeStructure
 * @see FieldUnivariateDerivative2
 * @see FieldUnivariateDerivative2
 * @see FieldGradient
 * @since 1.7
 */
public class UnivariateDerivative2 extends UnivariateDerivative<UnivariateDerivative2> {

    /** The constant value of π as a {@code UnivariateDerivative2}.
     * @since 2.0
     */
    public static final UnivariateDerivative2 PI = new UnivariateDerivative2(FastMath.PI, 0.0, 0.0);

    /** Serializable UID. */
    private static final long serialVersionUID = 20200520L;

    /** Value of the function. */
    private final double f0;

    /** First derivative of the function. */
    private final double f1;

    /** Second derivative of the function. */
    private final double f2;

    /** Build an instance with values and derivative.
     * @param f0 value of the function
     * @param f1 first derivative of the function
     * @param f2 second derivative of the function
     */
    public UnivariateDerivative2(final double f0, final double f1, final double f2) {
        this.f0 = f0;
        this.f1 = f1;
        this.f2 = f2;
    }

    /** Build an instance from a {@link DerivativeStructure}.
     * @param ds derivative structure
     * @exception MathIllegalArgumentException if either {@code ds} parameters
     * is not 1 or {@code ds} order is not 2
     */
    public UnivariateDerivative2(final DerivativeStructure ds) throws MathIllegalArgumentException {
        MathUtils.checkDimension(ds.getFreeParameters(), 1);
        MathUtils.checkDimension(ds.getOrder(), 2);
        this.f0 = ds.getValue();
        this.f1 = ds.getPartialDerivative(1);
        this.f2 = ds.getPartialDerivative(2);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 newInstance(final double value) {
        return new UnivariateDerivative2(value, 0.0, 0.0);
    }

    /** {@inheritDoc} */
    @Override
    public double getReal() {
        return getValue();
    }

    /** {@inheritDoc} */
    @Override
    public double getValue() {
        return f0;
    }

    /** {@inheritDoc} */
    @Override
    public double getDerivative(final int n) {
        switch (n) {
            case 0 :
                return f0;
            case 1 :
                return f1;
            case 2 :
                return f2;
            default :
                throw new MathIllegalArgumentException(LocalizedCoreFormats.DERIVATION_ORDER_NOT_ALLOWED, n);
        }
    }

    /** {@inheritDoc} */
    @Override
    public int getOrder() {
        return 2;
    }

    /** Get the first derivative.
     * @return first derivative
     * @see #getValue()
     * @see #getSecondDerivative()
     */
    public double getFirstDerivative() {
        return f1;
    }

    /** Get the second derivative.
     * @return second derivative
     * @see #getValue()
     * @see #getFirstDerivative()
     */
    public double getSecondDerivative() {
        return f2;
    }

    /** {@inheritDoc} */
    @Override
    public DerivativeStructure toDerivativeStructure() {
        return getField().getConversionFactory().build(f0, f1, f2);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 add(final double a) {
        return new UnivariateDerivative2(f0 + a, f1, f2);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 add(final UnivariateDerivative2 a) {
        return new UnivariateDerivative2(f0 + a.f0, f1 + a.f1, f2 + a.f2);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 subtract(final double a) {
        return new UnivariateDerivative2(f0 - a, f1, f2);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 subtract(final UnivariateDerivative2 a) {
        return new UnivariateDerivative2(f0 - a.f0, f1 - a.f1, f2 - a.f2);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 multiply(final int n) {
        return new UnivariateDerivative2(f0 * n, f1 * n, f2 * n);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 multiply(final double a) {
        return new UnivariateDerivative2(f0 * a, f1 * a, f2 * a);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 multiply(final UnivariateDerivative2 a) {
        return new UnivariateDerivative2(f0 * a.f0,
                                         MathArrays.linearCombination(f1, a.f0, f0, a.f1),
                                         MathArrays.linearCombination(f2, a.f0, 2 * f1, a.f1, f0, a.f2));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 divide(final double a) {
        final double inv1 = 1.0 / a;
        return new UnivariateDerivative2(f0 * inv1, f1 * inv1, f2 * inv1);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 divide(final UnivariateDerivative2 a) {
        final double inv1 = 1.0 / a.f0;
        final double inv2 = inv1 * inv1;
        final double inv3 = inv1 * inv2;
        return new UnivariateDerivative2(f0 * inv1,
                                         MathArrays.linearCombination(f1, a.f0, -f0, a.f1) * inv2,
                                         MathArrays.linearCombination(f2, a.f0 * a.f0,
                                                                      -2 * f1, a.f0 * a.f1,
                                                                       2 * f0, a.f1 * a.f1,
                                                                       -f0, a.f0 * a.f2) * inv3);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 remainder(final double a) {
        return new UnivariateDerivative2(FastMath.IEEEremainder(f0, a), f1, f2);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 remainder(final UnivariateDerivative2 a) {

        // compute k such that lhs % rhs = lhs - k rhs
        final double rem = FastMath.IEEEremainder(f0, a.f0);
        final double k   = FastMath.rint((f0 - rem) / a.f0);

        return new UnivariateDerivative2(rem, f1 - k * a.f1, f2 - k * a.f2);

    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 negate() {
        return new UnivariateDerivative2(-f0, -f1, -f2);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 abs() {
        if (Double.doubleToLongBits(f0) < 0) {
            // we use the bits representation to also handle -0.0
            return negate();
        } else {
            return this;
        }
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 ceil() {
        return new UnivariateDerivative2(FastMath.ceil(f0), 0.0, 0.0);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 floor() {
        return new UnivariateDerivative2(FastMath.floor(f0), 0.0, 0.0);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 rint() {
        return new UnivariateDerivative2(FastMath.rint(f0), 0.0, 0.0);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 sign() {
        return new UnivariateDerivative2(FastMath.signum(f0), 0.0, 0.0);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 copySign(final UnivariateDerivative2 sign) {
        long m = Double.doubleToLongBits(f0);
        long s = Double.doubleToLongBits(sign.f0);
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
            return this;
        }
        return negate(); // flip sign
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 copySign(final double sign) {
        long m = Double.doubleToLongBits(f0);
        long s = Double.doubleToLongBits(sign);
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
            return this;
        }
        return negate(); // flip sign
    }

    /** {@inheritDoc} */
    @Override
    public int getExponent() {
        return FastMath.getExponent(f0);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 scalb(final int n) {
        return new UnivariateDerivative2(FastMath.scalb(f0, n), FastMath.scalb(f1, n), FastMath.scalb(f2, n));
    }

    /** {@inheritDoc}
     * <p>
     * The {@code ulp} function is a step function, hence all its derivatives are 0.
     * </p>
     * @since 2.0
     */
    @Override
    public UnivariateDerivative2 ulp() {
        return new UnivariateDerivative2(FastMath.ulp(f0), 0.0, 0.0);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 hypot(final UnivariateDerivative2 y) {

        if (Double.isInfinite(f0) || Double.isInfinite(y.f0)) {
            return new UnivariateDerivative2(Double.POSITIVE_INFINITY, 0.0, 0.0);
        } else if (Double.isNaN(f0) || Double.isNaN(y.f0)) {
            return new UnivariateDerivative2(Double.NaN, 0.0, 0.0);
        } else {

            final int expX = getExponent();
            final int expY = y.getExponent();
            if (expX > expY + 27) {
                // y is neglectible with respect to x
                return abs();
            } else if (expY > expX + 27) {
                // x is neglectible with respect to y
                return y.abs();
            } else {

                // find an intermediate scale to avoid both overflow and underflow
                final int middleExp = (expX + expY) / 2;

                // scale parameters without losing precision
                final UnivariateDerivative2 scaledX = scalb(-middleExp);
                final UnivariateDerivative2 scaledY = y.scalb(-middleExp);

                // compute scaled hypotenuse
                final UnivariateDerivative2 scaledH =
                        scaledX.multiply(scaledX).add(scaledY.multiply(scaledY)).sqrt();

                // remove scaling
                return scaledH.scalb(middleExp);

            }

        }
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 reciprocal() {
        final double inv1 = 1.0 / f0;
        final double inv2 = inv1 * inv1;
        final double inv3 = inv1 * inv2;
        return new UnivariateDerivative2(inv1, -f1 * inv2, MathArrays.linearCombination(2 * f1, f1, -f0, f2) * inv3);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 compose(final double... f) {
        MathUtils.checkDimension(f.length, getOrder() + 1);
        return new UnivariateDerivative2(f[0],
                                         f[1] * f1,
                                         MathArrays.linearCombination(f[1], f2, f[2], f1 * f1));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 sqrt() {
        final double s0 = FastMath.sqrt(f0);
        final double s0twice = 2. * s0;
        final double s1 = f1 / s0twice;
        final double s2 = (f2 - 2. * s1 * s1) / s0twice;
        return new UnivariateDerivative2(s0, s1, s2);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 cbrt() {
        final double c  = FastMath.cbrt(f0);
        final double c2 = c * c;
        return compose(c, 1 / (3 * c2), -1 / (4.5 * c2 * f0));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 rootN(final int n) {
        if (n == 2) {
            return sqrt();
        } else if (n == 3) {
            return cbrt();
        } else {
            final double r = FastMath.pow(f0, 1.0 / n);
            final double z = n * FastMath.pow(r, n - 1);
            return compose(r, 1 / z, (1 - n) / (z * z * r));
        }
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2Field getField() {
        return UnivariateDerivative2Field.getInstance();
    }

    /** Compute a<sup>x</sup> where a is a double and x a {@link UnivariateDerivative2}
     * @param a number to exponentiate
     * @param x power to apply
     * @return a<sup>x</sup>
     */
    public static UnivariateDerivative2 pow(final double a, final UnivariateDerivative2 x) {
        if (a == 0) {
            return x.getField().getZero();
        } else {
            final double aX    = FastMath.pow(a, x.f0);
            final double lnA   = FastMath.log(a);
            final double aXlnA = aX * lnA;
            return new UnivariateDerivative2(aX, aXlnA * x.f1, aXlnA * (x.f1 * x.f1 * lnA + x.f2));
        }
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 pow(final double p) {
        if (p == 0) {
            return getField().getOne();
        } else {
            final double f0Pm2 = FastMath.pow(f0, p - 2);
            final double f0Pm1 = f0Pm2 * f0;
            final double f0P   = f0Pm1 * f0;
            return compose(f0P, p * f0Pm1, p * (p - 1) * f0Pm2);
        }
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 pow(final int n) {
        if (n == 0) {
            return getField().getOne();
        } else {
            final double f0Nm2 = FastMath.pow(f0, n - 2);
            final double f0Nm1 = f0Nm2 * f0;
            final double f0N   = f0Nm1 * f0;
            return compose(f0N, n * f0Nm1, n * (n - 1) * f0Nm2);
        }
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 pow(final UnivariateDerivative2 e) {
        return log().multiply(e).exp();
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 exp() {
        final double exp = FastMath.exp(f0);
        return compose(exp, exp, exp);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 expm1() {
        final double exp   = FastMath.exp(f0);
        final double expM1 = FastMath.expm1(f0);
        return compose(expM1, exp, exp);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 log() {
        final double inv = 1 / f0;
        return compose(FastMath.log(f0), inv, -inv * inv);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 log1p() {
        final double inv = 1 / (1 + f0);
        return compose(FastMath.log1p(f0), inv, -inv * inv);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 log10() {
        final double invF0 = 1 / f0;
        final double inv = invF0 / FastMath.log(10.0);
        return compose(FastMath.log10(f0), inv, -inv * invF0);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 cos() {
        final SinCos sinCos = FastMath.sinCos(f0);
        return compose(sinCos.cos(), -sinCos.sin(), -sinCos.cos());
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 sin() {
        final SinCos sinCos = FastMath.sinCos(f0);
        return compose(sinCos.sin(), sinCos.cos(), -sinCos.sin());
    }

    /** {@inheritDoc} */
    @Override
    public FieldSinCos<UnivariateDerivative2> sinCos() {
        final SinCos sinCos = FastMath.sinCos(f0);
        return new FieldSinCos<>(compose(sinCos.sin(),  sinCos.cos(), -sinCos.sin()),
                                 compose(sinCos.cos(), -sinCos.sin(), -sinCos.cos()));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 tan() {
        final double tan  = FastMath.tan(f0);
        final double sec2 = 1 + tan * tan;
        return compose(tan, sec2, 2 * sec2 * tan);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 acos() {
        final double inv = 1.0 / (1 - f0 * f0);
        final double mS  = -FastMath.sqrt(inv);
        return compose(FastMath.acos(f0), mS, mS * f0 * inv);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 asin() {
        final double inv = 1.0 / (1 - f0 * f0);
        final double s   = FastMath.sqrt(inv);
        return compose(FastMath.asin(f0), s, s * f0 * inv);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 atan() {
        final double inv = 1 / (1 + f0 * f0);
        return compose(FastMath.atan(f0), inv, -2 * f0 * inv * inv);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 atan2(final UnivariateDerivative2 x) {
        final double x2    = x.f0 * x.f0;
        final double f02   = f0 + f0;
        final double inv   = 1.0 / (f0 * f0 + x2);
        final double atan0 = FastMath.atan2(f0, x.f0);
        final double atan1 = MathArrays.linearCombination(x.f0, f1, -x.f1, f0) * inv;
        final double c     = MathArrays.linearCombination(f2, x2,
                                                          -2 * f1, x.f0 * x.f1,
                                                          f02, x.f1 * x.f1,
                                                          -f0, x.f0 * x.f2) * inv;
        return new UnivariateDerivative2(atan0, atan1, (c - f02 * atan1 * atan1) / x.f0);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 cosh() {
        final double c = FastMath.cosh(f0);
        final double s = FastMath.sinh(f0);
        return compose(c, s, c);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 sinh() {
        final double c = FastMath.cosh(f0);
        final double s = FastMath.sinh(f0);
        return compose(s, c, s);
    }

    /** {@inheritDoc} */
    @Override
    public FieldSinhCosh<UnivariateDerivative2> sinhCosh() {
        final SinhCosh sinhCosh = FastMath.sinhCosh(f0);
        return new FieldSinhCosh<>(compose(sinhCosh.sinh(), sinhCosh.cosh(), sinhCosh.sinh()),
                                   compose(sinhCosh.cosh(), sinhCosh.sinh(), sinhCosh.cosh()));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 tanh() {
        final double tanh  = FastMath.tanh(f0);
        final double sech2 = 1 - tanh * tanh;
        return compose(tanh, sech2, -2 * sech2 * tanh);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 acosh() {
        final double inv = 1 / (f0 * f0 - 1);
        final double s   = FastMath.sqrt(inv);
        return compose(FastMath.acosh(f0), s, -f0 * s * inv);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 asinh() {
        final double inv = 1 / (f0 * f0 + 1);
        final double s   = FastMath.sqrt(inv);
        return compose(FastMath.asinh(f0), s, -f0 * s * inv);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 atanh() {
        final double inv = 1 / (1 - f0 * f0);
        return compose(FastMath.atanh(f0), inv, 2 * f0 * inv * inv);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 toDegrees() {
        return new UnivariateDerivative2(FastMath.toDegrees(f0), FastMath.toDegrees(f1), FastMath.toDegrees(f2));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 toRadians() {
        return new UnivariateDerivative2(FastMath.toRadians(f0), FastMath.toRadians(f1), FastMath.toRadians(f2));
    }

    /** Evaluate Taylor expansion a univariate derivative.
     * @param delta parameter offset Δx
     * @return value of the Taylor expansion at x + Δx
     */
    public double taylor(final double delta) {
        return f0 + delta * (f1 + 0.5 * delta * f2);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 linearCombination(final UnivariateDerivative2[] a, final UnivariateDerivative2[] b) {

        // extract values and derivatives
        final int      n  = a.length;
        final double[] a0 = new double[n];
        final double[] b0 = new double[n];
        final double[] a1 = new double[2 * n];
        final double[] b1 = new double[2 * n];
        final double[] a2 = new double[3 * n];
        final double[] b2 = new double[3 * n];
        for (int i = 0; i < n; ++i) {
            final UnivariateDerivative2 ai = a[i];
            final UnivariateDerivative2 bi = b[i];
            a0[i]         = ai.f0;
            b0[i]         = bi.f0;
            a1[2 * i]     = ai.f0;
            a1[2 * i + 1] = ai.f1;
            b1[2 * i]     = bi.f1;
            b1[2 * i + 1] = bi.f0;
            a2[3 * i]     = ai.f0;
            a2[3 * i + 1] = ai.f1 + ai.f1;
            a2[3 * i + 2] = ai.f2;
            b2[3 * i]     = bi.f2;
            b2[3 * i + 1] = bi.f1;
            b2[3 * i + 2] = bi.f0;
        }

        return new UnivariateDerivative2(MathArrays.linearCombination(a0, b0),
                                         MathArrays.linearCombination(a1, b1),
                                         MathArrays.linearCombination(a2, b2));

    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 linearCombination(final double[] a, final UnivariateDerivative2[] b) {

        // extract values and derivatives
        final int      n  = b.length;
        final double[] b0 = new double[n];
        final double[] b1 = new double[n];
        final double[] b2 = new double[n];
        for (int i = 0; i < n; ++i) {
            b0[i] = b[i].f0;
            b1[i] = b[i].f1;
            b2[i] = b[i].f2;
        }

        return new UnivariateDerivative2(MathArrays.linearCombination(a, b0),
                                         MathArrays.linearCombination(a, b1),
                                         MathArrays.linearCombination(a, b2));

    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 linearCombination(final UnivariateDerivative2 a1, final UnivariateDerivative2 b1,
                                                   final UnivariateDerivative2 a2, final UnivariateDerivative2 b2) {
        return new UnivariateDerivative2(MathArrays.linearCombination(a1.f0, b1.f0,
                                                                      a2.f0, b2.f0),
                                         MathArrays.linearCombination(a1.f0, b1.f1,
                                                                      a1.f1, b1.f0,
                                                                      a2.f0, b2.f1,
                                                                      a2.f1, b2.f0),
                                         MathArrays.linearCombination(new double[] {
                                                                          a1.f0, 2 * a1.f1, a1.f2,
                                                                          a2.f0, 2 * a2.f1, a2.f2
                                                                      }, new double[] {
                                                                          b1.f2, b1.f1, b1.f0,
                                                                          b2.f2, b2.f1, b2.f0
                                                                      }));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 linearCombination(final double a1, final UnivariateDerivative2 b1,
                                                   final double a2, final UnivariateDerivative2 b2) {
        return new UnivariateDerivative2(MathArrays.linearCombination(a1, b1.f0,
                                                                      a2, b2.f0),
                                         MathArrays.linearCombination(a1, b1.f1,
                                                                      a2, b2.f1),
                                         MathArrays.linearCombination(a1, b1.f2,
                                                                      a2, b2.f2));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 linearCombination(final UnivariateDerivative2 a1, final UnivariateDerivative2 b1,
                                                   final UnivariateDerivative2 a2, final UnivariateDerivative2 b2,
                                                   final UnivariateDerivative2 a3, final UnivariateDerivative2 b3) {
        return new UnivariateDerivative2(MathArrays.linearCombination(a1.f0, b1.f0,
                                                                      a2.f0, b2.f0,
                                                                      a3.f0, b3.f0),
                                         MathArrays.linearCombination(new double[] {
                                                                          a1.f0, a1.f1,
                                                                          a2.f0, a2.f1,
                                                                          a3.f0, a3.f1
                                                                      }, new double[] {
                                                                          b1.f1, b1.f0,
                                                                          b2.f1, b2.f0,
                                                                          b3.f1, b3.f0
                                                                      }),
                                         MathArrays.linearCombination(new double[] {
                                                                          a1.f0, 2 * a1.f1, a1.f2,
                                                                          a2.f0, 2 * a2.f1, a2.f2,
                                                                          a3.f0, 2 * a3.f1, a3.f2
                                                                      }, new double[] {
                                                                          b1.f2, b1.f1, b1.f0,
                                                                          b2.f2, b2.f1, b2.f0,
                                                                          b3.f2, b3.f1, b3.f0
                                                                      }));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 linearCombination(final double a1, final UnivariateDerivative2 b1,
                                                   final double a2, final UnivariateDerivative2 b2,
                                                   final double a3, final UnivariateDerivative2 b3) {
        return new UnivariateDerivative2(MathArrays.linearCombination(a1, b1.f0,
                                                                      a2, b2.f0,
                                                                      a3, b3.f0),
                                         MathArrays.linearCombination(a1, b1.f1,
                                                                      a2, b2.f1,
                                                                      a3, b3.f1),
                                         MathArrays.linearCombination(a1, b1.f2,
                                                                      a2, b2.f2,
                                                                      a3, b3.f2));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 linearCombination(final UnivariateDerivative2 a1, final UnivariateDerivative2 b1,
                                                   final UnivariateDerivative2 a2, final UnivariateDerivative2 b2,
                                                   final UnivariateDerivative2 a3, final UnivariateDerivative2 b3,
                                                   final UnivariateDerivative2 a4, final UnivariateDerivative2 b4) {
        return new UnivariateDerivative2(MathArrays.linearCombination(a1.f0, b1.f0,
                                                                      a2.f0, b2.f0,
                                                                      a3.f0, b3.f0,
                                                                      a4.f0, b4.f0),
                                         MathArrays.linearCombination(new double[] {
                                                                          a1.f0, a1.f1,
                                                                          a2.f0, a2.f1,
                                                                          a3.f0, a3.f1,
                                                                          a4.f0, a4.f1
                                                                      }, new double[] {
                                                                          b1.f1, b1.f0,
                                                                          b2.f1, b2.f0,
                                                                          b3.f1, b3.f0,
                                                                          b4.f1, b4.f0
                                                                      }),
                                         MathArrays.linearCombination(new double[] {
                                                                          a1.f0, 2 * a1.f1, a1.f2,
                                                                          a2.f0, 2 * a2.f1, a2.f2,
                                                                          a3.f0, 2 * a3.f1, a3.f2,
                                                                          a4.f0, 2 * a4.f1, a4.f2
                                                                      }, new double[] {
                                                                          b1.f2, b1.f1, b1.f0,
                                                                          b2.f2, b2.f1, b2.f0,
                                                                          b3.f2, b3.f1, b3.f0,
                                                                          b4.f2, b4.f1, b4.f0
                                                                      }));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 linearCombination(final double a1, final UnivariateDerivative2 b1,
                                                   final double a2, final UnivariateDerivative2 b2,
                                                   final double a3, final UnivariateDerivative2 b3,
                                                   final double a4, final UnivariateDerivative2 b4) {
        return new UnivariateDerivative2(MathArrays.linearCombination(a1, b1.f0,
                                                                      a2, b2.f0,
                                                                      a3, b3.f0,
                                                                      a4, b4.f0),
                                         MathArrays.linearCombination(a1, b1.f1,
                                                                      a2, b2.f1,
                                                                      a3, b3.f1,
                                                                      a4, b4.f1),
                                         MathArrays.linearCombination(a1, b1.f2,
                                                                      a2, b2.f2,
                                                                      a3, b3.f2,
                                                                      a4, b4.f2));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative2 getPi() {
        return PI;
    }

    /** Test for the equality of two univariate derivatives.
     * <p>
     * univariate derivatives are considered equal if they have the same derivatives.
     * </p>
     * @param other Object to test for equality to this
     * @return true if two univariate derivatives are equal
     */
    @Override
    public boolean equals(Object other) {

        if (this == other) {
            return true;
        }

        if (other instanceof UnivariateDerivative2) {
            final UnivariateDerivative2 rhs = (UnivariateDerivative2) other;
            return f0 == rhs.f0 && f1 == rhs.f1 && f2 == rhs.f2;
        }

        return false;

    }

    /** Get a hashCode for the univariate derivative.
     * @return a hash code value for this object
     */
    @Override
    public int hashCode() {
        return 317 - 41 * Double.hashCode(f0) + 57 * Double.hashCode(f1) - 103 * Double.hashCode(f2);
    }

    /** {@inheritDoc}
     * <p>
     * Comparison performed considering that derivatives are intrinsically linked to monomials in the corresponding
     * Taylor expansion and that the higher the degree, the smaller the term.
     * </p>
     * @since 3.0
     */
    @Override
    public int compareTo(final UnivariateDerivative2 o) {
        final int cF0 = Double.compare(f0, o.getReal());
        if (cF0 == 0) {
            final int cF1 = Double.compare(f1, o.getFirstDerivative());
            if (cF1 == 0) {
                return Double.compare(f2, o.getSecondDerivative());
            } else {
                return cF1;
            }
        } else {
            return cF0;
        }
    }

}
