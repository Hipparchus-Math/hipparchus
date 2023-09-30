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
 * and {@link DerivativeStructure#getOrder() derivation order} also limited to one.
 * It should have less overhead than {@link DerivativeStructure} in its domain.</p>
 * <p>This class is an implementation of Rall's numbers. Rall's numbers are an
 * extension to the real numbers used throughout mathematical expressions; they hold
 * the derivative together with the value of a function.</p>
 * <p>{@link UnivariateDerivative1} instances can be used directly thanks to
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
 * @see FieldUnivariateDerivative1
 * @see FieldUnivariateDerivative2
 * @see FieldGradient
 * @since 1.7
 */
public class UnivariateDerivative1 extends UnivariateDerivative<UnivariateDerivative1> {

    /** The constant value of π as a {@code UnivariateDerivative1}.
     * @since 2.0
     */
    public static final UnivariateDerivative1 PI = new UnivariateDerivative1(FastMath.PI, 0.0);

    /** Serializable UID. */
    private static final long serialVersionUID = 20200519L;

    /** Value of the function. */
    private final double f0;

    /** First derivative of the function. */
    private final double f1;

    /** Build an instance with values and derivative.
     * @param f0 value of the function
     * @param f1 first derivative of the function
     */
    public UnivariateDerivative1(final double f0, final double f1) {
        this.f0 = f0;
        this.f1 = f1;
    }

    /** Build an instance from a {@link DerivativeStructure}.
     * @param ds derivative structure
     * @exception MathIllegalArgumentException if either {@code ds} parameters
     * is not 1 or {@code ds} order is not 1
     */
    public UnivariateDerivative1(final DerivativeStructure ds) throws MathIllegalArgumentException {
        MathUtils.checkDimension(ds.getFreeParameters(), 1);
        MathUtils.checkDimension(ds.getOrder(), 1);
        this.f0 = ds.getValue();
        this.f1 = ds.getPartialDerivative(1);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 newInstance(final double value) {
        return new UnivariateDerivative1(value, 0.0);
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
            default :
                throw new MathIllegalArgumentException(LocalizedCoreFormats.DERIVATION_ORDER_NOT_ALLOWED, n);
        }
    }

    /** {@inheritDoc} */
    @Override
    public int getOrder() {
        return 1;
    }

    /** Get the first derivative.
     * @return first derivative
     * @see #getValue()
     */
    public double getFirstDerivative() {
        return f1;
    }

    /** {@inheritDoc} */
    @Override
    public DerivativeStructure toDerivativeStructure() {
        return getField().getConversionFactory().build(f0, f1);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 add(final double a) {
        return new UnivariateDerivative1(f0 + a, f1);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 add(final UnivariateDerivative1 a) {
        return new UnivariateDerivative1(f0 + a.f0, f1 + a.f1);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 subtract(final double a) {
        return new UnivariateDerivative1(f0 - a, f1);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 subtract(final UnivariateDerivative1 a) {
        return new UnivariateDerivative1(f0 - a.f0, f1 - a.f1);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 multiply(final int n) {
        return new UnivariateDerivative1(f0 * n, f1 * n);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 multiply(final double a) {
        return new UnivariateDerivative1(f0 * a, f1 * a);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 multiply(final UnivariateDerivative1 a) {
        return new UnivariateDerivative1(f0 * a.f0,
                                         MathArrays.linearCombination(f1, a.f0, f0, a.f1));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 divide(final double a) {
        final double inv1 = 1.0 / a;
        return new UnivariateDerivative1(f0 * inv1, f1 * inv1);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 divide(final UnivariateDerivative1 a) {
        final double inv1 = 1.0 / a.f0;
        final double inv2 = inv1 * inv1;
        return new UnivariateDerivative1(f0 * inv1,
                                         MathArrays.linearCombination(f1, a.f0, -f0, a.f1) * inv2);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 remainder(final double a) {
        return new UnivariateDerivative1(FastMath.IEEEremainder(f0, a), f1);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 remainder(final UnivariateDerivative1 a) {

        // compute k such that lhs % rhs = lhs - k rhs
        final double rem = FastMath.IEEEremainder(f0, a.f0);
        final double k   = FastMath.rint((f0 - rem) / a.f0);

        return new UnivariateDerivative1(rem, f1 - k * a.f1);

    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 negate() {
        return new UnivariateDerivative1(-f0, -f1);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 abs() {
        if (Double.doubleToLongBits(f0) < 0) {
            // we use the bits representation to also handle -0.0
            return negate();
        } else {
            return this;
        }
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 ceil() {
        return new UnivariateDerivative1(FastMath.ceil(f0), 0.0);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 floor() {
        return new UnivariateDerivative1(FastMath.floor(f0), 0.0);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 rint() {
        return new UnivariateDerivative1(FastMath.rint(f0), 0.0);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 sign() {
        return new UnivariateDerivative1(FastMath.signum(f0), 0.0);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 copySign(final UnivariateDerivative1 sign) {
        long m = Double.doubleToLongBits(f0);
        long s = Double.doubleToLongBits(sign.f0);
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
            return this;
        }
        return negate(); // flip sign
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 copySign(final double sign) {
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
    public UnivariateDerivative1 scalb(final int n) {
        return new UnivariateDerivative1(FastMath.scalb(f0, n), FastMath.scalb(f1, n));
    }

    /** {@inheritDoc}
     * <p>
     * The {@code ulp} function is a step function, hence all its derivatives are 0.
     * </p>
     * @since 2.0
     */
    @Override
    public UnivariateDerivative1 ulp() {
        return new UnivariateDerivative1(FastMath.ulp(f0), 0.0);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 hypot(final UnivariateDerivative1 y) {

        if (Double.isInfinite(f0) || Double.isInfinite(y.f0)) {
            return new UnivariateDerivative1(Double.POSITIVE_INFINITY, 0.0);
        } else if (Double.isNaN(f0) || Double.isNaN(y.f0)) {
            return new UnivariateDerivative1(Double.NaN, 0.0);
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
                final UnivariateDerivative1 scaledX = scalb(-middleExp);
                final UnivariateDerivative1 scaledY = y.scalb(-middleExp);

                // compute scaled hypotenuse
                final UnivariateDerivative1 scaledH =
                        scaledX.multiply(scaledX).add(scaledY.multiply(scaledY)).sqrt();

                // remove scaling
                return scaledH.scalb(middleExp);

            }

        }
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 reciprocal() {
        final double inv1 = 1.0 / f0;
        final double inv2 = inv1 * inv1;
        return new UnivariateDerivative1(inv1, -f1 * inv2);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 compose(final double... f) {
        MathUtils.checkDimension(f.length, getOrder() + 1);
        return new UnivariateDerivative1(f[0], f[1] * f1);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 sqrt() {
        final double s = FastMath.sqrt(f0);
        return compose(s, 1 / (2 * s));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 cbrt() {
        final double c = FastMath.cbrt(f0);
        return compose(c, 1 / (3 * c * c));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 rootN(final int n) {
        if (n == 2) {
            return sqrt();
        } else if (n == 3) {
            return cbrt();
        } else {
            final double r = FastMath.pow(f0, 1.0 / n);
            return compose(r, 1 / (n * FastMath.pow(r, n - 1)));
        }
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1Field getField() {
        return UnivariateDerivative1Field.getInstance();
    }

    /** Compute a<sup>x</sup> where a is a double and x a {@link UnivariateDerivative1}
     * @param a number to exponentiate
     * @param x power to apply
     * @return a<sup>x</sup>
     */
    public static UnivariateDerivative1 pow(final double a, final UnivariateDerivative1 x) {
        if (a == 0) {
            return x.getField().getZero();
        } else {
            final double aX = FastMath.pow(a, x.f0);
            return new UnivariateDerivative1(aX, FastMath.log(a) * aX * x.f1);
        }
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 pow(final double p) {
        if (p == 0) {
            return getField().getOne();
        } else {
            final double f0Pm1 = FastMath.pow(f0, p - 1);
            return compose(f0Pm1 * f0, p * f0Pm1);
        }
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 pow(final int n) {
        if (n == 0) {
            return getField().getOne();
        } else {
            final double f0Nm1 = FastMath.pow(f0, n - 1);
            return compose(f0Nm1 * f0, n * f0Nm1);
        }
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 pow(final UnivariateDerivative1 e) {
        return log().multiply(e).exp();
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 exp() {
        final double exp = FastMath.exp(f0);
        return compose(exp, exp);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 expm1() {
        final double exp   = FastMath.exp(f0);
        final double expM1 = FastMath.expm1(f0);
        return compose(expM1, exp);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 log() {
        return compose(FastMath.log(f0), 1 / f0);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 log1p() {
        return compose(FastMath.log1p(f0), 1 / (1 + f0));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 log10() {
        return compose(FastMath.log10(f0), 1 / (f0 * FastMath.log(10.0)));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 cos() {
        final SinCos sinCos = FastMath.sinCos(f0);
        return compose(sinCos.cos(), -sinCos.sin());
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 sin() {
        final SinCos sinCos = FastMath.sinCos(f0);
        return compose(sinCos.sin(), sinCos.cos());
    }

    /** {@inheritDoc} */
    @Override
    public FieldSinCos<UnivariateDerivative1> sinCos() {
        final SinCos sinCos = FastMath.sinCos(f0);
        return new FieldSinCos<>(new UnivariateDerivative1(sinCos.sin(),  f1 * sinCos.cos()),
                                 new UnivariateDerivative1(sinCos.cos(), -f1 * sinCos.sin()));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 tan() {
        final double tan = FastMath.tan(f0);
        return compose(tan, 1 + tan * tan);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 acos() {
        return compose(FastMath.acos(f0), -1 / FastMath.sqrt(1 - f0 * f0));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 asin() {
        return compose(FastMath.asin(f0), 1 / FastMath.sqrt(1 - f0 * f0));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 atan() {
        return compose(FastMath.atan(f0), 1 / (1 + f0 * f0));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 atan2(final UnivariateDerivative1 x) {
        final double inv = 1.0 / (f0 * f0 + x.f0 * x.f0);
        return new UnivariateDerivative1(FastMath.atan2(f0, x.f0),
                                         MathArrays.linearCombination(x.f0, f1, -x.f1, f0) * inv);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 cosh() {
        return compose(FastMath.cosh(f0), FastMath.sinh(f0));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 sinh() {
        return compose(FastMath.sinh(f0), FastMath.cosh(f0));
    }

    /** {@inheritDoc} */
    @Override
    public FieldSinhCosh<UnivariateDerivative1> sinhCosh() {
        final SinhCosh sinhCosh = FastMath.sinhCosh(f0);
        return new FieldSinhCosh<>(new UnivariateDerivative1(sinhCosh.sinh(), f1 * sinhCosh.cosh()),
                                   new UnivariateDerivative1(sinhCosh.cosh(), f1 * sinhCosh.sinh()));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 tanh() {
        final double tanh = FastMath.tanh(f0);
        return compose(tanh, 1 - tanh * tanh);
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 acosh() {
        return compose(FastMath.acosh(f0), 1 / FastMath.sqrt(f0 * f0 - 1));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 asinh() {
        return compose(FastMath.asinh(f0), 1 / FastMath.sqrt(f0 * f0 + 1));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 atanh() {
        return compose(FastMath.atanh(f0), 1 / (1 - f0 * f0));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 toDegrees() {
        return new UnivariateDerivative1(FastMath.toDegrees(f0), FastMath.toDegrees(f1));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 toRadians() {
        return new UnivariateDerivative1(FastMath.toRadians(f0), FastMath.toRadians(f1));
    }

    /** Evaluate Taylor expansion a univariate derivative.
     * @param delta parameter offset Δx
     * @return value of the Taylor expansion at x + Δx
     */
    public double taylor(final double delta) {
        return f0 + delta * f1;
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 linearCombination(final UnivariateDerivative1[] a, final UnivariateDerivative1[] b) {

        // extract values and first derivatives
        final int      n  = a.length;
        final double[] a0 = new double[n];
        final double[] b0 = new double[n];
        final double[] a1 = new double[2 * n];
        final double[] b1 = new double[2 * n];
        for (int i = 0; i < n; ++i) {
            final UnivariateDerivative1 ai = a[i];
            final UnivariateDerivative1 bi = b[i];
            a0[i]         = ai.f0;
            b0[i]         = bi.f0;
            a1[2 * i]     = ai.f0;
            a1[2 * i + 1] = ai.f1;
            b1[2 * i]     = bi.f1;
            b1[2 * i + 1] = bi.f0;
        }

        return new UnivariateDerivative1(MathArrays.linearCombination(a0, b0),
                                         MathArrays.linearCombination(a1, b1));

    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 linearCombination(final double[] a, final UnivariateDerivative1[] b) {

        // extract values and first derivatives
        final int      n  = b.length;
        final double[] b0 = new double[n];
        final double[] b1 = new double[n];
        for (int i = 0; i < n; ++i) {
            b0[i] = b[i].f0;
            b1[i] = b[i].f1;
        }

        return new UnivariateDerivative1(MathArrays.linearCombination(a, b0),
                                         MathArrays.linearCombination(a, b1));

    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 linearCombination(final UnivariateDerivative1 a1, final UnivariateDerivative1 b1,
                                                   final UnivariateDerivative1 a2, final UnivariateDerivative1 b2) {
        return new UnivariateDerivative1(MathArrays.linearCombination(a1.f0, b1.f0,
                                                                      a2.f0, b2.f0),
                                         MathArrays.linearCombination(a1.f0, b1.f1,
                                                                      a1.f1, b1.f0,
                                                                      a2.f0, b2.f1,
                                                                      a2.f1, b2.f0));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 linearCombination(final double a1, final UnivariateDerivative1 b1,
                                                   final double a2, final UnivariateDerivative1 b2) {
        return new UnivariateDerivative1(MathArrays.linearCombination(a1, b1.f0,
                                                                      a2, b2.f0),
                                         MathArrays.linearCombination(a1, b1.f1,
                                                                      a2, b2.f1));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 linearCombination(final UnivariateDerivative1 a1, final UnivariateDerivative1 b1,
                                                   final UnivariateDerivative1 a2, final UnivariateDerivative1 b2,
                                                   final UnivariateDerivative1 a3, final UnivariateDerivative1 b3) {
        return new UnivariateDerivative1(MathArrays.linearCombination(a1.f0, b1.f0,
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
                                                                      }));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 linearCombination(final double a1, final UnivariateDerivative1 b1,
                                                   final double a2, final UnivariateDerivative1 b2,
                                                   final double a3, final UnivariateDerivative1 b3) {
        return new UnivariateDerivative1(MathArrays.linearCombination(a1, b1.f0,
                                                                      a2, b2.f0,
                                                                      a3, b3.f0),
                                         MathArrays.linearCombination(a1, b1.f1,
                                                                      a2, b2.f1,
                                                                      a3, b3.f1));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 linearCombination(final UnivariateDerivative1 a1, final UnivariateDerivative1 b1,
                                                   final UnivariateDerivative1 a2, final UnivariateDerivative1 b2,
                                                   final UnivariateDerivative1 a3, final UnivariateDerivative1 b3,
                                                   final UnivariateDerivative1 a4, final UnivariateDerivative1 b4) {
        return new UnivariateDerivative1(MathArrays.linearCombination(a1.f0, b1.f0,
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
                                                                      }));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 linearCombination(final double a1, final UnivariateDerivative1 b1,
                                                   final double a2, final UnivariateDerivative1 b2,
                                                   final double a3, final UnivariateDerivative1 b3,
                                                   final double a4, final UnivariateDerivative1 b4) {
        return new UnivariateDerivative1(MathArrays.linearCombination(a1, b1.f0,
                                                                      a2, b2.f0,
                                                                      a3, b3.f0,
                                                                      a4, b4.f0),
                                         MathArrays.linearCombination(a1, b1.f1,
                                                                      a2, b2.f1,
                                                                      a3, b3.f1,
                                                                      a4, b4.f1));
    }

    /** {@inheritDoc} */
    @Override
    public UnivariateDerivative1 getPi() {
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

        if (other instanceof UnivariateDerivative1) {
            final UnivariateDerivative1 rhs = (UnivariateDerivative1) other;
            return f0 == rhs.f0 && f1 == rhs.f1;
        }

        return false;

    }

    /** Get a hashCode for the univariate derivative.
     * @return a hash code value for this object
     */
    @Override
    public int hashCode() {
        return 453 - 19 * Double.hashCode(f0) + 37 * Double.hashCode(f1);
    }

    /** {@inheritDoc}
     * <p>
     * Comparison performed considering that derivatives are intrinsically linked to monomials in the corresponding
     * Taylor expansion and that the higher the degree, the smaller the term.
     * </p>
     * @since 3.0
     */
    @Override
    public int compareTo(final UnivariateDerivative1 o) {
        final int cF0 = Double.compare(f0, o.getReal());
        if (cF0 == 0) {
            return Double.compare(f1, o.getFirstDerivative());
        } else {
            return cF0;
        }
    }

}
