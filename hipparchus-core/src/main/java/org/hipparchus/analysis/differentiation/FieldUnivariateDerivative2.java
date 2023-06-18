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

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.FieldSinCos;
import org.hipparchus.util.FieldSinhCosh;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.MathUtils;

/** Class representing both the value and the differentials of a function.
 * <p>This class is a stripped-down version of {@link FieldDerivativeStructure}
 * with only one {@link FieldDerivativeStructure#getFreeParameters() free parameter}
 * and {@link FieldDerivativeStructure#getOrder() derivation order} limited to two.
 * It should have less overhead than {@link FieldDerivativeStructure} in its domain.</p>
 * <p>This class is an implementation of Rall's numbers. Rall's numbers are an
 * extension to the real numbers used throughout mathematical expressions; they hold
 * the derivative together with the value of a function.</p>
 * <p>{@link FieldUnivariateDerivative2} instances can be used directly thanks to
 * the arithmetic operators to the mathematical functions provided as
 * methods by this class (+, -, *, /, %, sin, cos ...).</p>
 * <p>Implementing complex expressions by hand using these classes is
 * a tedious and error-prone task but has the advantage of having no limitation
 * on the derivation order despite not requiring users to compute the derivatives by
 * themselves.</p>
 * <p>Instances of this class are guaranteed to be immutable.</p>
 * @param <T> the type of the function parameters and value
 * @see DerivativeStructure
 * @see UnivariateDerivative1
 * @see UnivariateDerivative2
 * @see Gradient
 * @see FieldDerivativeStructure
 * @see FieldUnivariateDerivative1
 * @see FieldGradient
 * @since 1.7
 */
public class FieldUnivariateDerivative2<T extends CalculusFieldElement<T>>
    extends FieldUnivariateDerivative<T, FieldUnivariateDerivative2<T>> {

    /** Value of the function. */
    private final T f0;

    /** First derivative of the function. */
    private final T f1;

    /** Second derivative of the function. */
    private final T f2;

    /** Build an instance with values and derivative.
     * @param f0 value of the function
     * @param f1 first derivative of the function
     * @param f2 second derivative of the function
     */
    public FieldUnivariateDerivative2(final T f0, final T f1, final T f2) {
        this.f0 = f0;
        this.f1 = f1;
        this.f2 = f2;
    }

    /** Build an instance from a {@link DerivativeStructure}.
     * @param ds derivative structure
     * @exception MathIllegalArgumentException if either {@code ds} parameters
     * is not 1 or {@code ds} order is not 2
     */
    public FieldUnivariateDerivative2(final FieldDerivativeStructure<T> ds) throws MathIllegalArgumentException {
        MathUtils.checkDimension(ds.getFreeParameters(), 1);
        MathUtils.checkDimension(ds.getOrder(), 2);
        this.f0 = ds.getValue();
        this.f1 = ds.getPartialDerivative(1);
        this.f2 = ds.getPartialDerivative(2);
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> newInstance(final double value) {
        final T zero = f0.getField().getZero();
        return new FieldUnivariateDerivative2<>(zero.newInstance(value), zero, zero);
    }

    /** {@inheritDoc} */
    @Override
    public double getReal() {
        return getValue().getReal();
    }

    /** Get the value part of the univariate derivative.
     * @return value part of the univariate derivative
     */
    @Override
    public T getValue() {
        return f0;
    }

    /** Get a derivative from the univariate derivative.
     * @param n derivation order (must be between 0 and {@link #getOrder()}, both inclusive)
     * @return n<sup>th</sup> derivative, or {@code NaN} if n is
     * either negative or strictly larger than {@link #getOrder()}
     */
    @Override
    public T getDerivative(final int n) {
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

    /** Get the derivation order.
     * @return derivation order
     */
    @Override
    public int getOrder() {
        return 2;
    }

    /** Get the first derivative.
     * @return first derivative
     * @see #getValue()
     */
    public T getFirstDerivative() {
        return f1;
    }

    /** Get the second derivative.
     * @return second derivative
     * @see #getValue()
     * @see #getFirstDerivative()
     */
    public T getSecondDerivative() {
        return f2;
    }

    /** Get the {@link Field} the value and parameters of the function belongs to.
     * @return {@link Field} the value and parameters of the function belongs to
     */
    public Field<T> getValueField() {
        return f0.getField();
    }

    /** Convert the instance to a {@link FieldDerivativeStructure}.
     * @return derivative structure with same value and derivative as the instance
     */
    @Override
    public FieldDerivativeStructure<T> toDerivativeStructure() {
        return getField().getConversionFactory().build(f0, f1, f2);
    }

    /** '+' operator.
     * @param a right hand side parameter of the operator
     * @return this+a
     */
    public FieldUnivariateDerivative2<T> add(final T a) {
        return new FieldUnivariateDerivative2<>(f0.add(a), f1, f2);
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> add(final double a) {
        return new FieldUnivariateDerivative2<>(f0.add(a), f1, f2);
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> add(final FieldUnivariateDerivative2<T> a) {
        return new FieldUnivariateDerivative2<>(f0.add(a.f0), f1.add(a.f1), f2.add(a.f2));
    }

    /** '-' operator.
     * @param a right hand side parameter of the operator
     * @return this-a
     */
    public FieldUnivariateDerivative2<T> subtract(final T a) {
        return new FieldUnivariateDerivative2<>(f0.subtract(a), f1, f2);
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> subtract(final double a) {
        return new FieldUnivariateDerivative2<>(f0.subtract(a), f1, f2);
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> subtract(final FieldUnivariateDerivative2<T> a) {
        return new FieldUnivariateDerivative2<>(f0.subtract(a.f0), f1.subtract(a.f1), f2.subtract(a.f2));
    }

    /** '&times;' operator.
     * @param a right hand side parameter of the operator
     * @return this&times;a
     */
    public FieldUnivariateDerivative2<T> multiply(final T a) {
        return new FieldUnivariateDerivative2<>(f0.multiply(a), f1.multiply(a), f2.multiply(a));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> multiply(final int n) {
        return new FieldUnivariateDerivative2<>(f0.multiply(n), f1.multiply(n), f2.multiply(n));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> multiply(final double a) {
        return new FieldUnivariateDerivative2<>(f0.multiply(a), f1.multiply(a), f2.multiply(a));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> multiply(final FieldUnivariateDerivative2<T> a) {
        return new FieldUnivariateDerivative2<>(f0.multiply(a.f0),
                                                a.f0.linearCombination(f1, a.f0, f0, a.f1),
                                                a.f0.linearCombination(f2, a.f0, f1.add(f1), a.f1, f0, a.f2));
    }

    /** '&divide;' operator.
     * @param a right hand side parameter of the operator
     * @return this&divide;a
     */
    public FieldUnivariateDerivative2<T> divide(final T a) {
        final T inv1 = a.reciprocal();
        return new FieldUnivariateDerivative2<>(f0.multiply(inv1), f1.multiply(inv1), f2.multiply(inv1));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> divide(final double a) {
        final double inv1 = 1.0 / a;
        return new FieldUnivariateDerivative2<>(f0.multiply(inv1), f1.multiply(inv1), f2.multiply(inv1));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> divide(final FieldUnivariateDerivative2<T> a) {
        final T inv1 = a.f0.reciprocal();
        final T inv2 = inv1.multiply(inv1);
        final T inv3 = inv1.multiply(inv2);
        return new FieldUnivariateDerivative2<>(f0.multiply(inv1),
                                                a.f0.linearCombination(f1, a.f0, f0.negate(), a.f1).multiply(inv2),
                                                a.f0.linearCombination(f2, a.f0.multiply(a.f0),
                                                                       f1.multiply(-2), a.f0.multiply(a.f1),
                                                                       f0.add(f0), a.f1.multiply(a.f1),
                                                                       f0.negate(), a.f0.multiply(a.f2)).multiply(inv3));
    }

    /** IEEE remainder operator.
     * @param a right hand side parameter of the operator
     * @return this - n &times; a where n is the closest integer to this/a
     * (the even integer is chosen for n if this/a is halfway between two integers)
     */
    public FieldUnivariateDerivative2<T> remainder(final T a) {
        return new FieldUnivariateDerivative2<>(FastMath.IEEEremainder(f0, a), f1, f2);
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> remainder(final double a) {
        return new FieldUnivariateDerivative2<>(FastMath.IEEEremainder(f0, a), f1, f2);
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> remainder(final FieldUnivariateDerivative2<T> a) {

        // compute k such that lhs % rhs = lhs - k rhs
        final T rem = FastMath.IEEEremainder(f0, a.f0);
        final T k   = FastMath.rint(f0.subtract(rem).divide(a.f0));

        return new FieldUnivariateDerivative2<>(rem,
                                                f1.subtract(k.multiply(a.f1)),
                                                f2.subtract(k.multiply(a.f2)));

    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> negate() {
        return new FieldUnivariateDerivative2<>(f0.negate(), f1.negate(), f2.negate());
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> abs() {
        if (Double.doubleToLongBits(f0.getReal()) < 0) {
            // we use the bits representation to also handle -0.0
            return negate();
        } else {
            return this;
        }
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> ceil() {
        final T zero = f0.getField().getZero();
        return new FieldUnivariateDerivative2<>(FastMath.ceil(f0), zero, zero);
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> floor() {
        final T zero = f0.getField().getZero();
        return new FieldUnivariateDerivative2<>(FastMath.floor(f0), zero, zero);
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> rint() {
        final T zero = f0.getField().getZero();
        return new FieldUnivariateDerivative2<>(FastMath.rint(f0), zero, zero);
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> sign() {
        final T zero = f0.getField().getZero();
        return new FieldUnivariateDerivative2<>(FastMath.sign(f0), zero, zero);
    }

    /**
     * Returns the instance with the sign of the argument.
     * A NaN {@code sign} argument is treated as positive.
     *
     * @param sign the sign for the returned value
     * @return the instance with the same sign as the {@code sign} argument
     */
    public FieldUnivariateDerivative2<T> copySign(final T sign) {
        long m = Double.doubleToLongBits(f0.getReal());
        long s = Double.doubleToLongBits(sign.getReal());
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
            return this;
        }
        return negate(); // flip sign
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> copySign(final FieldUnivariateDerivative2<T> sign) {
        long m = Double.doubleToLongBits(f0.getReal());
        long s = Double.doubleToLongBits(sign.f0.getReal());
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
            return this;
        }
        return negate(); // flip sign
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> copySign(final double sign) {
        long m = Double.doubleToLongBits(f0.getReal());
        long s = Double.doubleToLongBits(sign);
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
            return this;
        }
        return negate(); // flip sign
    }

    /** {@inheritDoc} */
    @Override
    public int getExponent() {
        return FastMath.getExponent(f0.getReal());
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> scalb(final int n) {
        return new FieldUnivariateDerivative2<>(FastMath.scalb(f0, n),
                                                FastMath.scalb(f1, n),
                                                FastMath.scalb(f2, n));
    }

    /** {@inheritDoc}
     * <p>
     * The {@code ulp} function is a step function, hence all its derivatives are 0.
     * </p>
     * @since 2.0
     */
    @Override
    public FieldUnivariateDerivative2<T> ulp() {
        final T zero = getValueField().getZero();
        return new FieldUnivariateDerivative2<>(FastMath.ulp(f0), zero, zero);
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> hypot(final FieldUnivariateDerivative2<T> y) {

        if (Double.isInfinite(f0.getReal()) || Double.isInfinite(y.f0.getReal())) {
            final T zero = f0.getField().getZero();
            return new FieldUnivariateDerivative2<>(f0.newInstance(Double.POSITIVE_INFINITY),
                                                    zero, zero);
        } else if (Double.isNaN(f0.getReal()) || Double.isNaN(y.f0.getReal())) {
            final T zero = f0.getField().getZero();
            return new FieldUnivariateDerivative2<>(f0.newInstance(Double.NaN),
                                                    zero, zero);
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
                final FieldUnivariateDerivative2<T> scaledX = scalb(-middleExp);
                final FieldUnivariateDerivative2<T> scaledY = y.scalb(-middleExp);

                // compute scaled hypotenuse
                final FieldUnivariateDerivative2<T> scaledH =
                        scaledX.multiply(scaledX).add(scaledY.multiply(scaledY)).sqrt();

                // remove scaling
                return scaledH.scalb(middleExp);

            }

        }
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> reciprocal() {
        final T inv1 = f0.reciprocal();
        final T inv2 = inv1.multiply(inv1);
        final T inv3 = inv1.multiply(inv2);
        return new FieldUnivariateDerivative2<>(inv1,
                                                f1.negate().multiply(inv2),
                                                f0.linearCombination(f1.add(f1), f1, f0.negate(), f2).multiply(inv3));
    }

    /** Compute composition of the instance by a function.
     * @param g0 value of the function at the current point (i.e. at {@code g(getValue())})
     * @param g1 first derivative of the function at the current point (i.e. at {@code g'(getValue())})
     * @param g2 second derivative of the function at the current point (i.e. at {@code g''(getValue())})
     * @return g(this)
     */
    public FieldUnivariateDerivative2<T> compose(final T g0, final T g1, final T g2) {
        return new FieldUnivariateDerivative2<>(g0,
                                                g1.multiply(f1),
                                                f0.linearCombination(g1, f2, g2, f1.multiply(f1)));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> sqrt() {
        final T s0 = FastMath.sqrt(f0);
        final T s0twice = s0.multiply(2);
        final T s1 = f1.divide(s0twice);
        final T s2 = (f2.subtract(s1.multiply(s1).multiply(2))).divide(s0twice);
        return new FieldUnivariateDerivative2<>(s0, s1, s2);
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> cbrt() {
        final T c  = FastMath.cbrt(f0);
        final T c2 = c.multiply(c);
        return compose(c, c2.multiply(3).reciprocal(), c2.multiply(-4.5).multiply(f0).reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> rootN(final int n) {
        if (n == 2) {
            return sqrt();
        } else if (n == 3) {
            return cbrt();
        } else {
            final T r = FastMath.pow(f0, 1.0 / n);
            final T z = FastMath.pow(r, n - 1).multiply(n);
            return compose(r, z.reciprocal(), z.multiply(z).multiply(r).reciprocal().multiply(1 -n));
        }
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2Field<T> getField() {
        return FieldUnivariateDerivative2Field.getUnivariateDerivative2Field(f0.getField());
    }

    /** Compute a<sup>x</sup> where a is a double and x a {@link FieldUnivariateDerivative2}
     * @param a number to exponentiate
     * @param x power to apply
     * @param <T> the type of the function parameters and value
     * @return a<sup>x</sup>
     */
    public static <T extends CalculusFieldElement<T>> FieldUnivariateDerivative2<T> pow(final double a, final FieldUnivariateDerivative2<T> x) {
        if (a == 0) {
            return x.getField().getZero();
        } else {
            final T      aX    = FastMath.pow(x.f0.newInstance(a), x.f0);
            final double lnA   = FastMath.log(a);
            final T      aXlnA = aX.multiply(lnA);
            return new FieldUnivariateDerivative2<>(aX,
                                                    aXlnA.multiply(x.f1),
                                                    aXlnA.multiply(x.f1.multiply(x.f1).multiply(lnA).add(x.f2)));
        }
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> pow(final double p) {
        if (p == 0) {
            return getField().getOne();
        } else {
            final T f0Pm2 = FastMath.pow(f0, p - 2);
            final T f0Pm1 = f0Pm2.multiply(f0);
            final T f0P   = f0Pm1.multiply(f0);
            return compose(f0P, f0Pm1.multiply(p), f0Pm2.multiply(p * (p - 1)));
        }
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> pow(final int n) {
        if (n == 0) {
            return getField().getOne();
        } else {
            final T f0Nm2 = FastMath.pow(f0, n - 2);
            final T f0Nm1 = f0Nm2.multiply(f0);
            final T f0N   = f0Nm1.multiply(f0);
            return compose(f0N, f0Nm1.multiply(n), f0Nm2.multiply(n * (n - 1)));
        }
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> pow(final FieldUnivariateDerivative2<T> e) {
        return log().multiply(e).exp();
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> exp() {
        final T exp = FastMath.exp(f0);
        return compose(exp, exp, exp);
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> expm1() {
        final T exp   = FastMath.exp(f0);
        final T expM1 = FastMath.expm1(f0);
        return compose(expM1, exp, exp);
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> log() {
        final T inv = f0.reciprocal();
        return compose(FastMath.log(f0), inv, inv.multiply(inv).negate());
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> log1p() {
        final T inv = f0.add(1).reciprocal();
        return compose(FastMath.log1p(f0), inv, inv.multiply(inv).negate());
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> log10() {
        final T invF0 = f0.reciprocal();
        final T inv = invF0.divide(FastMath.log(10.0));
        return compose(FastMath.log10(f0), inv, inv.multiply(invF0).negate());
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> cos() {
        final FieldSinCos<T> sinCos = FastMath.sinCos(f0);
        return compose(sinCos.cos(), sinCos.sin().negate(), sinCos.cos().negate());
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> sin() {
        final FieldSinCos<T> sinCos = FastMath.sinCos(f0);
        return compose(sinCos.sin(), sinCos.cos(), sinCos.sin().negate());
    }

    /** {@inheritDoc} */
    @Override
    public FieldSinCos<FieldUnivariateDerivative2<T>> sinCos() {
        final FieldSinCos<T> sinCos = FastMath.sinCos(f0);
        final T mSin = sinCos.sin().negate();
        final T mCos = sinCos.cos().negate();
        return new FieldSinCos<>(compose(sinCos.sin(), sinCos.cos(), mSin),
                                 compose(sinCos.cos(), mSin, mCos));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> tan() {
        final T tan  = FastMath.tan(f0);
        final T sec2 = tan.multiply(tan).add(1);
        return compose(tan, sec2, sec2.add(sec2).multiply(tan));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> acos() {
        final T inv = f0.multiply(f0).negate().add(1).reciprocal();
        final T mS  = inv.sqrt().negate();
        return compose(FastMath.acos(f0), mS, mS.multiply(f0).multiply(inv));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> asin() {
        final T inv = f0.multiply(f0).negate().add(1).reciprocal();
        final T s   = inv.sqrt();
        return compose(FastMath.asin(f0), s, s.multiply(f0).multiply(inv));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> atan() {
        final T inv = f0.multiply(f0).add(1).reciprocal();
        return compose(FastMath.atan(f0), inv, f0.multiply(-2).multiply(inv).multiply(inv));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> atan2(final FieldUnivariateDerivative2<T> x) {
        final T x2    = x.f0.multiply(x.f0);
        final T f02   = f0.add(f0);
        final T inv   = f0.multiply(f0).add(x2).reciprocal();
        final T atan0 = FastMath.atan2(f0, x.f0);
        final T atan1 = f0.linearCombination(x.f0, f1, x.f1.negate(), f0).multiply(inv);
        final T c     = f0.linearCombination(f2, x2,
                                             f1.multiply(-2), x.f0.multiply(x.f1),
                                             f02, x.f1.multiply(x.f1),
                                             f0.negate(), x.f0.multiply(x.f2)).multiply(inv);
        return new FieldUnivariateDerivative2<>(atan0,
                                                atan1,
                                                c.subtract(f02.multiply(atan1).multiply(atan1)).divide(x.f0));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> cosh() {
        final T c = FastMath.cosh(f0);
        final T s = FastMath.sinh(f0);
        return compose(c, s, c);
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> sinh() {
        final T c = FastMath.cosh(f0);
        final T s = FastMath.sinh(f0);
        return compose(s, c, s);
    }

    /** {@inheritDoc} */
    @Override
    public FieldSinhCosh<FieldUnivariateDerivative2<T>> sinhCosh() {
        final FieldSinhCosh<T> sinhCosh = FastMath.sinhCosh(f0);
        return new FieldSinhCosh<>(compose(sinhCosh.sinh(), sinhCosh.cosh(), sinhCosh.sinh()),
                                   compose(sinhCosh.cosh(), sinhCosh.sinh(), sinhCosh.cosh()));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> tanh() {
        final T tanh  = FastMath.tanh(f0);
        final T sech2 = tanh.multiply(tanh).negate().add(1);
        return compose(tanh, sech2, sech2.multiply(-2).multiply(tanh));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> acosh() {
        final T inv = f0.multiply(f0).subtract(1).reciprocal();
        final T s   = inv.sqrt();
        return compose(FastMath.acosh(f0), s, f0.negate().multiply(s).multiply(inv));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> asinh() {
        final T inv = f0.multiply(f0).add(1).reciprocal();
        final T s   = inv.sqrt();
        return compose(FastMath.asinh(f0), s, f0.negate().multiply(s).multiply(inv));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> atanh() {
        final T inv = f0.multiply(f0).negate().add(1).reciprocal();
        return compose(FastMath.atanh(f0), inv, f0.add(f0).multiply(inv).multiply(inv));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> toDegrees() {
        return new FieldUnivariateDerivative2<>(FastMath.toDegrees(f0), FastMath.toDegrees(f1), FastMath.toDegrees(f2));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> toRadians() {
        return new FieldUnivariateDerivative2<>(FastMath.toRadians(f0), FastMath.toRadians(f1), FastMath.toRadians(f2));
    }

    /** Evaluate Taylor expansion a univariate derivative.
     * @param delta parameter offset Δx
     * @return value of the Taylor expansion at x + Δx
     */
    public T taylor(final double delta) {
        return f0.add(f1.add(f2.multiply(0.5 * delta)).multiply(delta));
    }

    /** Evaluate Taylor expansion a univariate derivative.
     * @param delta parameter offset Δx
     * @return value of the Taylor expansion at x + Δx
     */
    public T taylor(final T delta) {
        return f0.add(f1.add(f2.multiply(delta.multiply(0.5))).multiply(delta));
    }

    /**
     * Compute a linear combination.
     * @param a Factors.
     * @param b Factors.
     * @return <code>&Sigma;<sub>i</sub> a<sub>i</sub> b<sub>i</sub></code>.
     * @throws MathIllegalArgumentException if arrays dimensions don't match
     */
    public FieldUnivariateDerivative2<T> linearCombination(final T[] a, final FieldUnivariateDerivative2<T>[] b) {

        // extract values and derivatives
        final Field<T> field = b[0].f0.getField();
        final int      n  = b.length;
        final T[] b0 = MathArrays.buildArray(field, n);
        final T[] b1 = MathArrays.buildArray(field, n);
        final T[] b2 = MathArrays.buildArray(field, n);
        for (int i = 0; i < n; ++i) {
            b0[i] = b[i].f0;
            b1[i] = b[i].f1;
            b2[i] = b[i].f2;
        }

        return new FieldUnivariateDerivative2<>(b[0].f0.linearCombination(a, b0),
                                                b[0].f0.linearCombination(a, b1),
                                                b[0].f0.linearCombination(a, b2));

    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> linearCombination(final FieldUnivariateDerivative2<T>[] a,
                                                           final FieldUnivariateDerivative2<T>[] b) {

        // extract values and derivatives
        final Field<T> field = a[0].f0.getField();
        final int n  = a.length;
        final T[] a0 = MathArrays.buildArray(field, n);
        final T[] b0 = MathArrays.buildArray(field, n);
        final T[] a1 = MathArrays.buildArray(field, 2 * n);
        final T[] b1 = MathArrays.buildArray(field, 2 * n);
        final T[] a2 = MathArrays.buildArray(field, 3 * n);
        final T[] b2 = MathArrays.buildArray(field, 3 * n);
        for (int i = 0; i < n; ++i) {
            final FieldUnivariateDerivative2<T> ai = a[i];
            final FieldUnivariateDerivative2<T> bi = b[i];
            a0[i]         = ai.f0;
            b0[i]         = bi.f0;
            a1[2 * i]     = ai.f0;
            a1[2 * i + 1] = ai.f1;
            b1[2 * i]     = bi.f1;
            b1[2 * i + 1] = bi.f0;
            a2[3 * i]     = ai.f0;
            a2[3 * i + 1] = ai.f1.add(ai.f1);
            a2[3 * i + 2] = ai.f2;
            b2[3 * i]     = bi.f2;
            b2[3 * i + 1] = bi.f1;
            b2[3 * i + 2] = bi.f0;
        }

        return new FieldUnivariateDerivative2<>(a[0].f0.linearCombination(a0, b0),
                                                a[0].f0.linearCombination(a1, b1),
                                                a[0].f0.linearCombination(a2, b2));

    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> linearCombination(final double[] a, final FieldUnivariateDerivative2<T>[] b) {

        // extract values and derivatives
        final Field<T> field = b[0].f0.getField();
        final int      n  = b.length;
        final T[] b0 = MathArrays.buildArray(field, n);
        final T[] b1 = MathArrays.buildArray(field, n);
        final T[] b2 = MathArrays.buildArray(field, n);
        for (int i = 0; i < n; ++i) {
            b0[i] = b[i].f0;
            b1[i] = b[i].f1;
            b2[i] = b[i].f2;
        }

        return new FieldUnivariateDerivative2<>(b[0].f0.linearCombination(a, b0),
                                                b[0].f0.linearCombination(a, b1),
                                                b[0].f0.linearCombination(a, b2));

    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> linearCombination(final FieldUnivariateDerivative2<T> a1, final FieldUnivariateDerivative2<T> b1,
                                                           final FieldUnivariateDerivative2<T> a2, final FieldUnivariateDerivative2<T> b2) {
        final Field<T> field = a1.f0.getField();
        final T[]      u2    = MathArrays.buildArray(field, 6);
        final T[]      v2    = MathArrays.buildArray(field, 6);
        u2[0] = a1.f0;
        u2[1] = a1.f1.add(a1.f1);
        u2[2] = a1.f2;
        u2[3] = a2.f0;
        u2[4] = a2.f1.add(a2.f1);
        u2[5] = a2.f2;
        v2[0] = b1.f2;
        v2[1] = b1.f1;
        v2[2] = b1.f0;
        v2[3] = b2.f2;
        v2[4] = b2.f1;
        v2[5] = b2.f0;
        return new FieldUnivariateDerivative2<>(a1.f0.linearCombination(a1.f0, b1.f0,
                                                                        a2.f0, b2.f0),
                                                a1.f0.linearCombination(a1.f0, b1.f1,
                                                                        a1.f1, b1.f0,
                                                                        a2.f0, b2.f1,
                                                                        a2.f1, b2.f0),
                                                a1.f0.linearCombination(u2, v2));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> linearCombination(final double a1, final FieldUnivariateDerivative2<T> b1,
                                                           final double a2, final FieldUnivariateDerivative2<T> b2) {
        return new FieldUnivariateDerivative2<>(b1.f0.linearCombination(a1, b1.f0,
                                                                        a2, b2.f0),
                                                b1.f0.linearCombination(a1, b1.f1,
                                                                        a2, b2.f1),
                                                b1.f0.linearCombination(a1, b1.f2,
                                                                        a2, b2.f2));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> linearCombination(final FieldUnivariateDerivative2<T> a1, final FieldUnivariateDerivative2<T> b1,
                                                           final FieldUnivariateDerivative2<T> a2, final FieldUnivariateDerivative2<T> b2,
                                                           final FieldUnivariateDerivative2<T> a3, final FieldUnivariateDerivative2<T> b3) {
        final Field<T> field = a1.f0.getField();
        final T[]      u1     = MathArrays.buildArray(field, 6);
        final T[]      v1     = MathArrays.buildArray(field, 6);
        u1[0] = a1.f0;
        u1[1] = a1.f1;
        u1[2] = a2.f0;
        u1[3] = a2.f1;
        u1[4] = a3.f0;
        u1[5] = a3.f1;
        v1[0] = b1.f1;
        v1[1] = b1.f0;
        v1[2] = b2.f1;
        v1[3] = b2.f0;
        v1[4] = b3.f1;
        v1[5] = b3.f0;
        final T[]      u2     = MathArrays.buildArray(field, 9);
        final T[]      v2     = MathArrays.buildArray(field, 9);
        u2[0] = a1.f0;
        u2[1] = a1.f1.add(a1.f1);
        u2[2] = a1.f2;
        u2[3] = a2.f0;
        u2[4] = a2.f1.add(a2.f1);
        u2[5] = a2.f2;
        u2[6] = a3.f0;
        u2[7] = a3.f1.add(a3.f1);
        u2[8] = a3.f2;
        v2[0] = b1.f2;
        v2[1] = b1.f1;
        v2[2] = b1.f0;
        v2[3] = b2.f2;
        v2[4] = b2.f1;
        v2[5] = b2.f0;
        v2[6] = b3.f2;
        v2[7] = b3.f1;
        v2[8] = b3.f0;
        return new FieldUnivariateDerivative2<>(a1.f0.linearCombination(a1.f0, b1.f0,
                                                                        a2.f0, b2.f0,
                                                                        a3.f0, b3.f0),
                                                a1.f0.linearCombination(u1, v1),
                                                a1.f0.linearCombination(u2, v2));
    }

    /**
     * Compute a linear combination.
     * @param a1 first factor of the first term
     * @param b1 second factor of the first term
     * @param a2 first factor of the second term
     * @param b2 second factor of the second term
     * @param a3 first factor of the third term
     * @param b3 second factor of the third term
     * @return a<sub>1</sub>&times;b<sub>1</sub> +
     * a<sub>2</sub>&times;b<sub>2</sub> + a<sub>3</sub>&times;b<sub>3</sub>
     * @see #linearCombination(double, FieldUnivariateDerivative2, double, FieldUnivariateDerivative2)
     * @see #linearCombination(double, FieldUnivariateDerivative2, double, FieldUnivariateDerivative2, double, FieldUnivariateDerivative2, double, FieldUnivariateDerivative2)
     * @exception MathIllegalArgumentException if number of free parameters or orders are inconsistent
     */
    public FieldUnivariateDerivative2<T> linearCombination(final T a1, final FieldUnivariateDerivative2<T> b1,
                                                           final T a2, final FieldUnivariateDerivative2<T> b2,
                                                           final T a3, final FieldUnivariateDerivative2<T> b3) {
        return new FieldUnivariateDerivative2<>(b1.f0.linearCombination(a1, b1.f0,
                                                                        a2, b2.f0,
                                                                        a3, b3.f0),
                                                b1.f0.linearCombination(a1, b1.f1,
                                                                        a2, b2.f1,
                                                                        a3, b3.f1),
                                                b1.f0.linearCombination(a1, b1.f2,
                                                                        a2, b2.f2,
                                                                        a3, b3.f2));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> linearCombination(final double a1, final FieldUnivariateDerivative2<T> b1,
                                                           final double a2, final FieldUnivariateDerivative2<T> b2,
                                                           final double a3, final FieldUnivariateDerivative2<T> b3) {
        return new FieldUnivariateDerivative2<>(b1.f0.linearCombination(a1, b1.f0,
                                                                        a2, b2.f0,
                                                                        a3, b3.f0),
                                                b1.f0.linearCombination(a1, b1.f1,
                                                                        a2, b2.f1,
                                                                        a3, b3.f1),
                                                b1.f0.linearCombination(a1, b1.f2,
                                                                        a2, b2.f2,
                                                                        a3, b3.f2));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> linearCombination(final FieldUnivariateDerivative2<T> a1, final FieldUnivariateDerivative2<T> b1,
                                                           final FieldUnivariateDerivative2<T> a2, final FieldUnivariateDerivative2<T> b2,
                                                           final FieldUnivariateDerivative2<T> a3, final FieldUnivariateDerivative2<T> b3,
                                                           final FieldUnivariateDerivative2<T> a4, final FieldUnivariateDerivative2<T> b4) {
        final Field<T> field = a1.f0.getField();
        final T[] u1 = MathArrays.buildArray(field, 8);
        final T[] v1 = MathArrays.buildArray(field, 8);
        u1[0] = a1.f0;
        u1[1] = a1.f1;
        u1[2] = a2.f0;
        u1[3] = a2.f1;
        u1[4] = a3.f0;
        u1[5] = a3.f1;
        u1[6] = a4.f0;
        u1[7] = a4.f1;
        v1[0] = b1.f1;
        v1[1] = b1.f0;
        v1[2] = b2.f1;
        v1[3] = b2.f0;
        v1[4] = b3.f1;
        v1[5] = b3.f0;
        v1[6] = b4.f1;
        v1[7] = b4.f0;
        final T[] u2 = MathArrays.buildArray(field, 12);
        final T[] v2 = MathArrays.buildArray(field, 12);
        u2[ 0] = a1.f0;
        u2[ 1] = a1.f1.add(a1.f1);
        u2[ 2] = a1.f2;
        u2[ 3] = a2.f0;
        u2[ 4] = a2.f1.add(a2.f1);
        u2[ 5] = a2.f2;
        u2[ 6] = a3.f0;
        u2[ 7] = a3.f1.add(a3.f1);
        u2[ 8] = a3.f2;
        u2[ 9] = a4.f0;
        u2[10] = a4.f1.add(a4.f1);
        u2[11] = a4.f2;
        v2[ 0] = b1.f2;
        v2[ 1] = b1.f1;
        v2[ 2] = b1.f0;
        v2[ 3] = b2.f2;
        v2[ 4] = b2.f1;
        v2[ 5] = b2.f0;
        v2[ 6] = b3.f2;
        v2[ 7] = b3.f1;
        v2[ 8] = b3.f0;
        v2[ 9] = b4.f2;
        v2[10] = b4.f1;
        v2[11] = b4.f0;
        return new FieldUnivariateDerivative2<>(a1.f0.linearCombination(a1.f0, b1.f0,
                                                                        a2.f0, b2.f0,
                                                                        a3.f0, b3.f0,
                                                                        a4.f0, b4.f0),
                                                a1.f0.linearCombination(u1, v1),
                                                a1.f0.linearCombination(u2, v2));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> linearCombination(final double a1, final FieldUnivariateDerivative2<T> b1,
                                                           final double a2, final FieldUnivariateDerivative2<T> b2,
                                                           final double a3, final FieldUnivariateDerivative2<T> b3,
                                                           final double a4, final FieldUnivariateDerivative2<T> b4) {
        return new FieldUnivariateDerivative2<>(b1.f0.linearCombination(a1, b1.f0,
                                                                        a2, b2.f0,
                                                                        a3, b3.f0,
                                                                        a4, b4.f0),
                                                b1.f0.linearCombination(a1, b1.f1,
                                                                        a2, b2.f1,
                                                                        a3, b3.f1,
                                                                        a4, b4.f1),
                                                b1.f0.linearCombination(a1, b1.f2,
                                                                        a2, b2.f2,
                                                                        a3, b3.f2,
                                                                        a4, b4.f2));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative2<T> getPi() {
        final T zero = getValueField().getZero();
        return new FieldUnivariateDerivative2<>(zero.getPi(), zero, zero);
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

        if (other instanceof FieldUnivariateDerivative2) {
            @SuppressWarnings("unchecked")
            final FieldUnivariateDerivative2<T> rhs = (FieldUnivariateDerivative2<T>) other;
            return f0.equals(rhs.f0) && f1.equals(rhs.f1) && f2.equals(rhs.f2);
        }

        return false;

    }

    /** Get a hashCode for the univariate derivative.
     * @return a hash code value for this object
     */
    @Override
    public int hashCode() {
        return 317 - 41 * f0.hashCode() + 57 * f1.hashCode() - 103 * f2.hashCode();
    }

}
