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
 * and {@link FieldDerivativeStructure#getOrder() derivation order} also limited to one.
 * It should have less overhead than {@link FieldDerivativeStructure} in its domain.</p>
 * <p>This class is an implementation of Rall's numbers. Rall's numbers are an
 * extension to the real numbers used throughout mathematical expressions; they hold
 * the derivative together with the value of a function.</p>
 * <p>{@link FieldUnivariateDerivative1} instances can be used directly thanks to
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
 * @see FieldUnivariateDerivative2
 * @see FieldGradient
 * @since 1.7
 */
public class FieldUnivariateDerivative1<T extends CalculusFieldElement<T>>
    extends FieldUnivariateDerivative<T, FieldUnivariateDerivative1<T>> {

    /** Value of the function. */
    private final T f0;

    /** First derivative of the function. */
    private final T f1;

    /** Build an instance with values and derivative.
     * @param f0 value of the function
     * @param f1 first derivative of the function
     */
    public FieldUnivariateDerivative1(final T f0, final T f1) {
        this.f0 = f0;
        this.f1 = f1;
    }

    /** Build an instance from a {@link DerivativeStructure}.
     * @param ds derivative structure
     * @exception MathIllegalArgumentException if either {@code ds} parameters
     * is not 1 or {@code ds} order is not 1
     */
    public FieldUnivariateDerivative1(final FieldDerivativeStructure<T> ds) throws MathIllegalArgumentException {
        MathUtils.checkDimension(ds.getFreeParameters(), 1);
        MathUtils.checkDimension(ds.getOrder(), 1);
        this.f0 = ds.getValue();
        this.f1 = ds.getPartialDerivative(1);
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> newInstance(final double value) {
        final T zero = f0.getField().getZero();
        return new FieldUnivariateDerivative1<>(zero.newInstance(value), zero);
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
            default :
                throw new MathIllegalArgumentException(LocalizedCoreFormats.DERIVATION_ORDER_NOT_ALLOWED, n);
        }
    }

    /** Get the derivation order.
     * @return derivation order
     */
    @Override
    public int getOrder() {
        return 1;
    }

    /** Get the first derivative.
     * @return first derivative
     * @see #getValue()
     */
    public T getFirstDerivative() {
        return f1;
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
        return getField().getConversionFactory().build(f0, f1);
    }

    /** '+' operator.
     * @param a right hand side parameter of the operator
     * @return this+a
     */
    public FieldUnivariateDerivative1<T> add(final T a) {
        return new FieldUnivariateDerivative1<>(f0.add(a), f1);
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> add(final double a) {
        return new FieldUnivariateDerivative1<>(f0.add(a), f1);
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> add(final FieldUnivariateDerivative1<T> a) {
        return new FieldUnivariateDerivative1<>(f0.add(a.f0), f1.add(a.f1));
    }

    /** '-' operator.
     * @param a right hand side parameter of the operator
     * @return this-a
     */
    public FieldUnivariateDerivative1<T> subtract(final T a) {
        return new FieldUnivariateDerivative1<>(f0.subtract(a), f1);
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> subtract(final double a) {
        return new FieldUnivariateDerivative1<>(f0.subtract(a), f1);
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> subtract(final FieldUnivariateDerivative1<T> a) {
        return new FieldUnivariateDerivative1<>(f0.subtract(a.f0), f1.subtract(a.f1));
    }

    /** '&times;' operator.
     * @param a right hand side parameter of the operator
     * @return this&times;a
     */
    public FieldUnivariateDerivative1<T> multiply(final T a) {
        return new FieldUnivariateDerivative1<>(f0.multiply(a), f1.multiply(a));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> multiply(final int n) {
        return new FieldUnivariateDerivative1<>(f0.multiply(n), f1.multiply(n));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> multiply(final double a) {
        return new FieldUnivariateDerivative1<>(f0.multiply(a), f1.multiply(a));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> multiply(final FieldUnivariateDerivative1<T> a) {
        return new FieldUnivariateDerivative1<>(f0.multiply(a.f0),
                                                a.f0.linearCombination(f1, a.f0, f0, a.f1));
    }

    /** '&divide;' operator.
     * @param a right hand side parameter of the operator
     * @return this&divide;a
     */
    public FieldUnivariateDerivative1<T> divide(final T a) {
        final T inv1 = a.reciprocal();
        return new FieldUnivariateDerivative1<>(f0.multiply(inv1), f1.multiply(inv1));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> divide(final double a) {
        final double inv1 = 1.0 / a;
        return new FieldUnivariateDerivative1<>(f0.multiply(inv1), f1.multiply(inv1));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> divide(final FieldUnivariateDerivative1<T> a) {
        final T inv1 = a.f0.reciprocal();
        final T inv2 = inv1.multiply(inv1);
        return new FieldUnivariateDerivative1<>(f0.multiply(inv1),
                                                a.f0.linearCombination(f1, a.f0, f0.negate(), a.f1).multiply(inv2));
    }

    /** IEEE remainder operator.
     * @param a right hand side parameter of the operator
     * @return this - n &times; a where n is the closest integer to this/a
     * (the even integer is chosen for n if this/a is halfway between two integers)
     */
    public FieldUnivariateDerivative1<T> remainder(final T a) {
        return new FieldUnivariateDerivative1<>(FastMath.IEEEremainder(f0, a), f1);
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> remainder(final double a) {
        return new FieldUnivariateDerivative1<>(FastMath.IEEEremainder(f0, a), f1);
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> remainder(final FieldUnivariateDerivative1<T> a) {

        // compute k such that lhs % rhs = lhs - k rhs
        final T rem = FastMath.IEEEremainder(f0, a.f0);
        final T k   = FastMath.rint(f0.subtract(rem).divide(a.f0));

        return new FieldUnivariateDerivative1<>(rem, f1.subtract(k.multiply(a.f1)));

    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> negate() {
        return new FieldUnivariateDerivative1<>(f0.negate(), f1.negate());
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> abs() {
        if (Double.doubleToLongBits(f0.getReal()) < 0) {
            // we use the bits representation to also handle -0.0
            return negate();
        } else {
            return this;
        }
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> ceil() {
        return new FieldUnivariateDerivative1<>(FastMath.ceil(f0), f0.getField().getZero());
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> floor() {
        return new FieldUnivariateDerivative1<>(FastMath.floor(f0), f0.getField().getZero());
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> rint() {
        return new FieldUnivariateDerivative1<>(FastMath.rint(f0), f0.getField().getZero());
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> sign() {
        return new FieldUnivariateDerivative1<>(FastMath.sign(f0), f0.getField().getZero());
    }

    /**
     * Returns the instance with the sign of the argument.
     * A NaN {@code sign} argument is treated as positive.
     *
     * @param sign the sign for the returned value
     * @return the instance with the same sign as the {@code sign} argument
     */
    public FieldUnivariateDerivative1<T> copySign(final T sign) {
        long m = Double.doubleToLongBits(f0.getReal());
        long s = Double.doubleToLongBits(sign.getReal());
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
            return this;
        }
        return negate(); // flip sign
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> copySign(final FieldUnivariateDerivative1<T> sign) {
        long m = Double.doubleToLongBits(f0.getReal());
        long s = Double.doubleToLongBits(sign.f0.getReal());
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
            return this;
        }
        return negate(); // flip sign
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> copySign(final double sign) {
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
    public FieldUnivariateDerivative1<T> scalb(final int n) {
        return new FieldUnivariateDerivative1<>(FastMath.scalb(f0, n), FastMath.scalb(f1, n));
    }

    /** {@inheritDoc}
     * <p>
     * The {@code ulp} function is a step function, hence all its derivatives are 0.
     * </p>
     * @since 2.0
     */
    @Override
    public FieldUnivariateDerivative1<T> ulp() {
        return new FieldUnivariateDerivative1<>(FastMath.ulp(f0), getValueField().getZero());
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> hypot(final FieldUnivariateDerivative1<T> y) {

        if (Double.isInfinite(f0.getReal()) || Double.isInfinite(y.f0.getReal())) {
            return new FieldUnivariateDerivative1<>(f0.newInstance(Double.POSITIVE_INFINITY),
                                                    f0.getField().getZero());
        } else if (Double.isNaN(f0.getReal()) || Double.isNaN(y.f0.getReal())) {
            return new FieldUnivariateDerivative1<>(f0.newInstance(Double.NaN),
                                                    f0.getField().getZero());
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
                final FieldUnivariateDerivative1<T> scaledX = scalb(-middleExp);
                final FieldUnivariateDerivative1<T> scaledY = y.scalb(-middleExp);

                // compute scaled hypotenuse
                final FieldUnivariateDerivative1<T> scaledH =
                        scaledX.multiply(scaledX).add(scaledY.multiply(scaledY)).sqrt();

                // remove scaling
                return scaledH.scalb(middleExp);

            }

        }
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> reciprocal() {
        final T inv1 = f0.reciprocal();
        final T inv2 = inv1.multiply(inv1);
        return new FieldUnivariateDerivative1<>(inv1, f1.negate().multiply(inv2));
    }

    /** Compute composition of the instance by a function.
     * @param g0 value of the function at the current point (i.e. at {@code g(getValue())})
     * @param g1 first derivative of the function at the current point (i.e. at {@code g'(getValue())})
     * @return g(this)
     */
    public FieldUnivariateDerivative1<T> compose(final T g0, final T g1) {
        return new FieldUnivariateDerivative1<>(g0, g1.multiply(f1));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> sqrt() {
        final T s = FastMath.sqrt(f0);
        return compose(s, s.add(s).reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> cbrt() {
        final T c = FastMath.cbrt(f0);
        return compose(c, c.multiply(c).multiply(3).reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> rootN(final int n) {
        if (n == 2) {
            return sqrt();
        } else if (n == 3) {
            return cbrt();
        } else {
            final T r = FastMath.pow(f0, 1.0 / n);
            return compose(r, FastMath.pow(r, n - 1).multiply(n).reciprocal());
        }
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1Field<T> getField() {
        return FieldUnivariateDerivative1Field.getUnivariateDerivative1Field(f0.getField());
    }

    /** Compute a<sup>x</sup> where a is a double and x a {@link FieldUnivariateDerivative1}
     * @param a number to exponentiate
     * @param x power to apply
     * @param <T> the type of the function parameters and value
     * @return a<sup>x</sup>
     */
    public static <T extends CalculusFieldElement<T>> FieldUnivariateDerivative1<T> pow(final double a, final FieldUnivariateDerivative1<T> x) {
        if (a == 0) {
            return x.getField().getZero();
        } else {
            final T aX = FastMath.pow(x.f0.newInstance(a), x.f0);
            return new FieldUnivariateDerivative1<>(aX, aX.multiply(FastMath.log(a)).multiply(x.f1));
        }
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> pow(final double p) {
        if (p == 0) {
            return getField().getOne();
        } else {
            final T f0Pm1 = FastMath.pow(f0, p - 1);
            return compose(f0Pm1.multiply(f0), f0Pm1.multiply(p));
        }
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> pow(final int n) {
        if (n == 0) {
            return getField().getOne();
        } else {
            final T f0Nm1 = FastMath.pow(f0, n - 1);
            return compose(f0Nm1.multiply(f0), f0Nm1.multiply(n));
        }
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> pow(final FieldUnivariateDerivative1<T> e) {
        return log().multiply(e).exp();
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> exp() {
        final T exp = FastMath.exp(f0);
        return compose(exp, exp);
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> expm1() {
        final T exp   = FastMath.exp(f0);
        final T expM1 = FastMath.expm1(f0);
        return compose(expM1, exp);
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> log() {
        return compose(FastMath.log(f0), f0.reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> log1p() {
        return compose(FastMath.log1p(f0), f0.add(1).reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> log10() {
        return compose(FastMath.log10(f0), f0.multiply(FastMath.log(10.0)).reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> cos() {
        final FieldSinCos<T> sinCos = FastMath.sinCos(f0);
        return compose(sinCos.cos(), sinCos.sin().negate());
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> sin() {
        final FieldSinCos<T> sinCos = FastMath.sinCos(f0);
        return compose(sinCos.sin(), sinCos.cos());
    }

    /** {@inheritDoc} */
    @Override
    public FieldSinCos<FieldUnivariateDerivative1<T>> sinCos() {
        final FieldSinCos<T> sinCos = FastMath.sinCos(f0);
        return new FieldSinCos<>(new FieldUnivariateDerivative1<>(sinCos.sin(), f1.multiply(sinCos.cos())),
                                 new FieldUnivariateDerivative1<>(sinCos.cos(), f1.multiply(sinCos.sin()).negate()));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> tan() {
        final T tan = FastMath.tan(f0);
        return compose(tan, tan.multiply(tan).add(1));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> acos() {
        return compose(FastMath.acos(f0), f0.multiply(f0).negate().add(1).sqrt().reciprocal().negate());
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> asin() {
        return compose(FastMath.asin(f0), f0.multiply(f0).negate().add(1).sqrt().reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> atan() {
        return compose(FastMath.atan(f0), f0.multiply(f0).add(1).reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> atan2(final FieldUnivariateDerivative1<T> x) {
        final T inv = f0.multiply(f0).add(x.f0.multiply(x.f0)).reciprocal();
        return new FieldUnivariateDerivative1<>(FastMath.atan2(f0, x.f0),
                                                f0.linearCombination(x.f0, f1, x.f1.negate(), f0).multiply(inv));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> cosh() {
        return compose(FastMath.cosh(f0), FastMath.sinh(f0));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> sinh() {
        return compose(FastMath.sinh(f0), FastMath.cosh(f0));
    }

    /** {@inheritDoc} */
    @Override
    public FieldSinhCosh<FieldUnivariateDerivative1<T>> sinhCosh() {
        final FieldSinhCosh<T> sinhCosh = FastMath.sinhCosh(f0);
        return new FieldSinhCosh<>(new FieldUnivariateDerivative1<>(sinhCosh.sinh(), f1.multiply(sinhCosh.cosh())),
                                   new FieldUnivariateDerivative1<>(sinhCosh.cosh(), f1.multiply(sinhCosh.sinh())));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> tanh() {
        final T tanh = FastMath.tanh(f0);
        return compose(tanh, tanh.multiply(tanh).negate().add(1));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> acosh() {
        return compose(FastMath.acosh(f0), f0.multiply(f0).subtract(1).sqrt().reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> asinh() {
        return compose(FastMath.asinh(f0), f0.multiply(f0).add(1).sqrt().reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> atanh() {
        return compose(FastMath.atanh(f0), f0.multiply(f0).negate().add(1).reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> toDegrees() {
        return new FieldUnivariateDerivative1<>(FastMath.toDegrees(f0), FastMath.toDegrees(f1));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> toRadians() {
        return new FieldUnivariateDerivative1<>(FastMath.toRadians(f0), FastMath.toRadians(f1));
    }

    /** Evaluate Taylor expansion of a univariate derivative.
     * @param delta parameter offset Δx
     * @return value of the Taylor expansion at x + Δx
     */
    public T taylor(final double delta) {
        return f0.add(f1.multiply(delta));
    }

    /** Evaluate Taylor expansion of a univariate derivative.
     * @param delta parameter offset Δx
     * @return value of the Taylor expansion at x + Δx
     */
    public T taylor(final T delta) {
        return f0.add(f1.multiply(delta));
    }

    /**
     * Compute a linear combination.
     * @param a Factors.
     * @param b Factors.
     * @return <code>&Sigma;<sub>i</sub> a<sub>i</sub> b<sub>i</sub></code>.
     * @throws MathIllegalArgumentException if arrays dimensions don't match
     */
    public FieldUnivariateDerivative1<T> linearCombination(final T[] a, final FieldUnivariateDerivative1<T>[] b) {

        // extract values and first derivatives
        final Field<T> field = b[0].f0.getField();
        final int      n  = b.length;
        final T[] b0 = MathArrays.buildArray(field, n);
        final T[] b1 = MathArrays.buildArray(field, n);
        for (int i = 0; i < n; ++i) {
            b0[i] = b[i].f0;
            b1[i] = b[i].f1;
        }

        return new FieldUnivariateDerivative1<>(b[0].f0.linearCombination(a, b0),
                                                b[0].f0.linearCombination(a, b1));

    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> linearCombination(final FieldUnivariateDerivative1<T>[] a,
                                                           final FieldUnivariateDerivative1<T>[] b) {

        // extract values and first derivatives
        final Field<T> field = a[0].f0.getField();
        final int n  = a.length;
        final T[] a0 = MathArrays.buildArray(field, n);
        final T[] b0 = MathArrays.buildArray(field, n);
        final T[] a1 = MathArrays.buildArray(field, 2 * n);
        final T[] b1 = MathArrays.buildArray(field, 2 * n);
        for (int i = 0; i < n; ++i) {
            final FieldUnivariateDerivative1<T> ai = a[i];
            final FieldUnivariateDerivative1<T> bi = b[i];
            a0[i]         = ai.f0;
            b0[i]         = bi.f0;
            a1[2 * i]     = ai.f0;
            a1[2 * i + 1] = ai.f1;
            b1[2 * i]     = bi.f1;
            b1[2 * i + 1] = bi.f0;
        }

        return new FieldUnivariateDerivative1<>(a[0].f0.linearCombination(a0, b0),
                                                a[0].f0.linearCombination(a1, b1));

    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> linearCombination(final double[] a, final FieldUnivariateDerivative1<T>[] b) {

        // extract values and first derivatives
        final Field<T> field = b[0].f0.getField();
        final int      n  = b.length;
        final T[] b0 = MathArrays.buildArray(field, n);
        final T[] b1 = MathArrays.buildArray(field, n);
        for (int i = 0; i < n; ++i) {
            b0[i] = b[i].f0;
            b1[i] = b[i].f1;
        }

        return new FieldUnivariateDerivative1<>(b[0].f0.linearCombination(a, b0),
                                                b[0].f0.linearCombination(a, b1));

    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> linearCombination(final FieldUnivariateDerivative1<T> a1, final FieldUnivariateDerivative1<T> b1,
                                                           final FieldUnivariateDerivative1<T> a2, final FieldUnivariateDerivative1<T> b2) {
        return new FieldUnivariateDerivative1<>(a1.f0.linearCombination(a1.f0, b1.f0,
                                                                        a2.f0, b2.f0),
                                                a1.f0.linearCombination(a1.f0, b1.f1,
                                                                        a1.f1, b1.f0,
                                                                        a2.f0, b2.f1,
                                                                        a2.f1, b2.f0));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> linearCombination(final double a1, final FieldUnivariateDerivative1<T> b1,
                                                           final double a2, final FieldUnivariateDerivative1<T> b2) {
        return new FieldUnivariateDerivative1<>(b1.f0.linearCombination(a1, b1.f0,
                                                                        a2, b2.f0),
                                                b1.f0.linearCombination(a1, b1.f1,
                                                                        a2, b2.f1));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> linearCombination(final FieldUnivariateDerivative1<T> a1, final FieldUnivariateDerivative1<T> b1,
                                                           final FieldUnivariateDerivative1<T> a2, final FieldUnivariateDerivative1<T> b2,
                                                           final FieldUnivariateDerivative1<T> a3, final FieldUnivariateDerivative1<T> b3) {
        final Field<T> field = a1.f0.getField();
        final T[] a = MathArrays.buildArray(field, 6);
        final T[] b = MathArrays.buildArray(field, 6);
        a[0] = a1.f0;
        a[1] = a1.f1;
        a[2] = a2.f0;
        a[3] = a2.f1;
        a[4] = a3.f0;
        a[5] = a3.f1;
        b[0] = b1.f1;
        b[1] = b1.f0;
        b[2] = b2.f1;
        b[3] = b2.f0;
        b[4] = b3.f1;
        b[5] = b3.f0;
        return new FieldUnivariateDerivative1<>(a1.f0.linearCombination(a1.f0, b1.f0,
                                                                        a2.f0, b2.f0,
                                                                        a3.f0, b3.f0),
                                                a1.f0.linearCombination(a, b));
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
     * @see #linearCombination(double, FieldUnivariateDerivative1, double, FieldUnivariateDerivative1)
     * @see #linearCombination(double, FieldUnivariateDerivative1, double, FieldUnivariateDerivative1, double, FieldUnivariateDerivative1, double, FieldUnivariateDerivative1)
     * @exception MathIllegalArgumentException if number of free parameters or orders are inconsistent
     */
    public FieldUnivariateDerivative1<T> linearCombination(final T a1, final FieldUnivariateDerivative1<T> b1,
                                                           final T a2, final FieldUnivariateDerivative1<T> b2,
                                                           final T a3, final FieldUnivariateDerivative1<T> b3) {
        return new FieldUnivariateDerivative1<>(b1.f0.linearCombination(a1, b1.f0,
                                                                        a2, b2.f0,
                                                                        a3, b3.f0),
                                                b1.f0.linearCombination(a1, b1.f1,
                                                                        a2, b2.f1,
                                                                        a3, b3.f1));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> linearCombination(final double a1, final FieldUnivariateDerivative1<T> b1,
                                                           final double a2, final FieldUnivariateDerivative1<T> b2,
                                                           final double a3, final FieldUnivariateDerivative1<T> b3) {
        return new FieldUnivariateDerivative1<>(b1.f0.linearCombination(a1, b1.f0,
                                                                        a2, b2.f0,
                                                                        a3, b3.f0),
                                                b1.f0.linearCombination(a1, b1.f1,
                                                                        a2, b2.f1,
                                                                        a3, b3.f1));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> linearCombination(final FieldUnivariateDerivative1<T> a1, final FieldUnivariateDerivative1<T> b1,
                                                           final FieldUnivariateDerivative1<T> a2, final FieldUnivariateDerivative1<T> b2,
                                                           final FieldUnivariateDerivative1<T> a3, final FieldUnivariateDerivative1<T> b3,
                                                           final FieldUnivariateDerivative1<T> a4, final FieldUnivariateDerivative1<T> b4) {
        final Field<T> field = a1.f0.getField();
        final T[] a = MathArrays.buildArray(field, 8);
        final T[] b = MathArrays.buildArray(field, 8);
        a[0] = a1.f0;
        a[1] = a1.f1;
        a[2] = a2.f0;
        a[3] = a2.f1;
        a[4] = a3.f0;
        a[5] = a3.f1;
        a[6] = a4.f0;
        a[7] = a4.f1;
        b[0] = b1.f1;
        b[1] = b1.f0;
        b[2] = b2.f1;
        b[3] = b2.f0;
        b[4] = b3.f1;
        b[5] = b3.f0;
        b[6] = b4.f1;
        b[7] = b4.f0;
        return new FieldUnivariateDerivative1<>(a1.f0.linearCombination(a1.f0, b1.f0,
                                                                        a2.f0, b2.f0,
                                                                        a3.f0, b3.f0,
                                                                        a4.f0, b4.f0),
                                                a1.f0.linearCombination(a, b));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> linearCombination(final double a1, final FieldUnivariateDerivative1<T> b1,
                                                           final double a2, final FieldUnivariateDerivative1<T> b2,
                                                           final double a3, final FieldUnivariateDerivative1<T> b3,
                                                           final double a4, final FieldUnivariateDerivative1<T> b4) {
        return new FieldUnivariateDerivative1<>(b1.f0.linearCombination(a1, b1.f0,
                                                                        a2, b2.f0,
                                                                        a3, b3.f0,
                                                                        a4, b4.f0),
                                                b1.f0.linearCombination(a1, b1.f1,
                                                                        a2, b2.f1,
                                                                        a3, b3.f1,
                                                                        a4, b4.f1));
    }

    /** {@inheritDoc} */
    @Override
    public FieldUnivariateDerivative1<T> getPi() {
        final T zero = getValueField().getZero();
        return new FieldUnivariateDerivative1<>(zero.getPi(), zero);
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

        if (other instanceof FieldUnivariateDerivative1) {
            @SuppressWarnings("unchecked")
            final FieldUnivariateDerivative1<T> rhs = (FieldUnivariateDerivative1<T>) other;
            return f0.equals(rhs.f0) && f1.equals(rhs.f1);
        }

        return false;

    }

    /** Get a hashCode for the univariate derivative.
     * @return a hash code value for this object
     */
    @Override
    public int hashCode() {
        return 453 - 19 * f0.hashCode() + 37 * f1.hashCode();
    }

}
