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

import java.util.Arrays;

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
 * with {@link FieldDerivativeStructure#getOrder() derivation order} limited to one.
 * It should have less overhead than {@link FieldDerivativeStructure} in its domain.</p>
 * <p>This class is an implementation of Rall's numbers. Rall's numbers are an
 * extension to the real numbers used throughout mathematical expressions; they hold
 * the derivative together with the value of a function.</p>
 * <p>{@link FieldGradient} instances can be used directly thanks to
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
 * @see FieldUnivariateDerivative2
 * @since 1.7
 */
public class FieldGradient<T extends CalculusFieldElement<T>> implements FieldDerivative<T, FieldGradient<T>> {

    /** Value of the function. */
    private final T value;

    /** Gradient of the function. */
    private final T[] grad;

    /** Build an instance with values and unitialized derivatives array.
     * @param value value of the function
     * @param freeParameters number of free parameters
     */
    private FieldGradient(final T value, int freeParameters) {
        this.value = value;
        this.grad  = MathArrays.buildArray(value.getField(), freeParameters);
    }

    /** Build an instance with values and derivative.
     * @param value value of the function
     * @param gradient gradient of the function
     */
    @SafeVarargs
    public FieldGradient(final T value, final T... gradient) {
        this(value, gradient.length);
        System.arraycopy(gradient, 0, grad, 0, grad.length);
    }

    /** Build an instance from a {@link DerivativeStructure}.
     * @param ds derivative structure
     * @exception MathIllegalArgumentException if {@code ds} order
     * is not 1
     */
    public FieldGradient(final FieldDerivativeStructure<T> ds) throws MathIllegalArgumentException {
        this(ds.getValue(), ds.getFreeParameters());
        MathUtils.checkDimension(ds.getOrder(), 1);
        System.arraycopy(ds.getAllDerivatives(), 1, grad, 0, grad.length);
    }

    /** Build an instance corresponding to a constant value.
     * @param freeParameters number of free parameters (i.e. dimension of the gradient)
     * @param value constant value of the function
     * @param <T> the type of the function parameters and value
     * @return a {@code FieldGradient} with a constant value and all derivatives set to 0.0
     */
    public static <T extends CalculusFieldElement<T>> FieldGradient<T> constant(final int freeParameters, final T value) {
        final FieldGradient<T> g = new FieldGradient<>(value, freeParameters);
        Arrays.fill(g.grad, value.getField().getZero());
        return g;
    }

    /** Build a {@code Gradient} representing a variable.
     * <p>Instances built using this method are considered
     * to be the free variables with respect to which differentials
     * are computed. As such, their differential with respect to
     * themselves is +1.</p>
     * @param freeParameters number of free parameters (i.e. dimension of the gradient)
     * @param index index of the variable (from 0 to {@link #getFreeParameters() getFreeParameters()} - 1)
     * @param value value of the variable
     * @param <T> the type of the function parameters and value
     * @return a {@code FieldGradient} with a constant value and all derivatives set to 0.0 except the
     * one at {@code index} which will be set to 1.0
     */
    public static <T extends CalculusFieldElement<T>> FieldGradient<T> variable(final int freeParameters,
                                                                                final int index, final T value) {
        final FieldGradient<T> g = new FieldGradient<>(value, freeParameters);
        final Field<T> field = value.getField();
        Arrays.fill(g.grad, field.getZero());
        g.grad[index] = field.getOne();
        return g;
    }

    /** Get the {@link Field} the value and parameters of the function belongs to.
     * @return {@link Field} the value and parameters of the function belongs to
     */
    public Field<T> getValueField() {
        return value.getField();
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> newInstance(final double c) {
        return newInstance(getValueField().getZero().newInstance(c));
    }

    /** Create an instance corresponding to a constant real value.
     * <p>
     * The default implementation creates the instance by adding
     * the value to {@code getField().getZero()}. This is not optimal
     * and does not work when called with a negative zero as the
     * sign of zero is lost with the addition. The default implementation
     * should therefore be overridden in concrete classes. The default
     * implementation will be removed at the next major version.
     * </p>
     * @param c constant real value
     * @return instance corresponding to a constant real value
     */
    public FieldGradient<T> newInstance(final T c) {
        return new FieldGradient<>(c, MathArrays.buildArray(value.getField(), grad.length));
    }

    /** {@inheritDoc} */
    @Override
    public double getReal() {
        return getValue().getReal();
    }

    /** Get the value part of the function.
     * @return value part of the value of the function
     */
    @Override
    public T getValue() {
        return value;
    }

    /** Get the gradient part of the function.
     * @return gradient part of the value of the function
     */
    public T[] getGradient() {
        return grad.clone();
    }

    /** Get the number of free parameters.
     * @return number of free parameters
     */
    @Override
    public int getFreeParameters() {
        return grad.length;
    }

    /** {@inheritDoc} */
    @Override
    public int getOrder() {
        return 1;
    }

    /** {@inheritDoc} */
    @Override
    public T getPartialDerivative(final int ... orders)
        throws MathIllegalArgumentException {

        // check the number of components
        if (orders.length != grad.length) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   orders.length, grad.length);
        }

        // check that either all derivation orders are set to 0,
        // or that only one is set to 1 and all other ones are set to 0
        int selected = -1;
        for (int i = 0; i < orders.length; ++i) {
            if (orders[i] != 0) {
                if (selected >= 0 || orders[i] != 1) {
                     throw new MathIllegalArgumentException(LocalizedCoreFormats.DERIVATION_ORDER_NOT_ALLOWED,
                                                           orders[i]);
                }
                // found the component set to derivation order 1
                selected = i;
            }
        }

        return (selected < 0) ? value : grad[selected];

    }

    /** Get the partial derivative with respect to one parameter.
     * @param n index of the parameter (counting from 0)
     * @return partial derivative with respect to the n<sup>th</sup> parameter
     * @exception MathIllegalArgumentException if n is either negative or larger
     * or equal to {@link #getFreeParameters()}
     */
    public T getPartialDerivative(final int n) throws MathIllegalArgumentException {
        if (n < 0 || n >= grad.length) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.OUT_OF_RANGE_SIMPLE, n, 0, grad.length - 1);
        }
        return grad[n];
    }

    /** Convert the instance to a {@link FieldDerivativeStructure}.
     * @return derivative structure with same value and derivative as the instance
     */
    public FieldDerivativeStructure<T> toDerivativeStructure() {
        final T[] derivatives = MathArrays.buildArray(getValueField(), 1 + grad.length);
        derivatives[0] = value;
        System.arraycopy(grad, 0, derivatives, 1, grad.length);
        return getField().getConversionFactory().build(derivatives);
    }

    /** '+' operator.
     * @param a right hand side parameter of the operator
     * @return this+a
     */
    public FieldGradient<T> add(final T a) {
        return new FieldGradient<>(value.add(a), grad);
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> add(final double a) {
        return new FieldGradient<>(value.add(a), grad);
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> add(final FieldGradient<T> a) {
        final FieldGradient<T> result = newInstance(value.add(a.value));
        for (int i = 0; i < grad.length; ++i) {
            result.grad[i] = grad[i].add(a.grad[i]);
        }
        return result;
    }

    /** '-' operator.
     * @param a right hand side parameter of the operator
     * @return this-a
     */
    public FieldGradient<T> subtract(final T a) {
        return new FieldGradient<>(value.subtract(a), grad);
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> subtract(final double a) {
        return new FieldGradient<>(value.subtract(a), grad);
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> subtract(final FieldGradient<T> a) {
        final FieldGradient<T> result = newInstance(value.subtract(a.value));
        for (int i = 0; i < grad.length; ++i) {
            result.grad[i] = grad[i].subtract(a.grad[i]);
        }
        return result;
    }

    /** '&times;' operator.
     * @param n right hand side parameter of the operator
     * @return this&times;n
     */
    public FieldGradient<T> multiply(final T n) {
        final FieldGradient<T> result = newInstance(value.multiply(n));
        for (int i = 0; i < grad.length; ++i) {
            result.grad[i] = grad[i].multiply(n);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> multiply(final int n) {
        final FieldGradient<T> result = newInstance(value.multiply(n));
        for (int i = 0; i < grad.length; ++i) {
            result.grad[i] = grad[i].multiply(n);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> multiply(final double a) {
        final FieldGradient<T> result = newInstance(value.multiply(a));
        for (int i = 0; i < grad.length; ++i) {
            result.grad[i] = grad[i].multiply(a);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> multiply(final FieldGradient<T> a) {
        final FieldGradient<T> result = newInstance(value.multiply(a.value));
        for (int i = 0; i < grad.length; ++i) {
            result.grad[i] = grad[i].multiply(a.value).add(value.multiply(a.grad[i]));
        }
        return result;
    }

    /** '&divide;' operator.
     * @param a right hand side parameter of the operator
     * @return this&divide;a
     */
    public FieldGradient<T> divide(final T a) {
        final FieldGradient<T> result = newInstance(value.divide(a));
        for (int i = 0; i < grad.length; ++i) {
            result.grad[i] = grad[i].divide(a);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> divide(final double a) {
        final FieldGradient<T> result = newInstance(value.divide(a));
        for (int i = 0; i < grad.length; ++i) {
            result.grad[i] = grad[i].divide(a);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> divide(final FieldGradient<T> a) {
        final T inv1 = a.value.reciprocal();
        final T inv2 = inv1.multiply(inv1);
        final FieldGradient<T> result = newInstance(value.multiply(inv1));
        for (int i = 0; i < grad.length; ++i) {
            result.grad[i] = grad[i].multiply(a.value).subtract(value.multiply(a.grad[i])).multiply(inv2);
        }
        return result;
    }

    /** IEEE remainder operator.
     * @param a right hand side parameter of the operator
     * @return this - n &times; a where n is the closest integer to this/a
     * (the even integer is chosen for n if this/a is halfway between two integers)
     */
    public FieldGradient<T> remainder(final T a) {
        return new FieldGradient<>(FastMath.IEEEremainder(value, a), grad);
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> remainder(final double a) {
        return new FieldGradient<>(FastMath.IEEEremainder(value, a), grad);
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> remainder(final FieldGradient<T> a) {

        // compute k such that lhs % rhs = lhs - k rhs
        final T rem = FastMath.IEEEremainder(value, a.value);
        final T k   = FastMath.rint(value.subtract(rem).divide(a.value));

        final FieldGradient<T> result = newInstance(rem);
        for (int i = 0; i < grad.length; ++i) {
            result.grad[i] = grad[i].subtract(k.multiply(a.grad[i]));
        }
        return result;

    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> negate() {
        final FieldGradient<T> result = newInstance(value.negate());
        for (int i = 0; i < grad.length; ++i) {
            result.grad[i] = grad[i].negate();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> abs() {
        if (Double.doubleToLongBits(value.getReal()) < 0) {
            // we use the bits representation to also handle -0.0
            return negate();
        } else {
            return this;
        }
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> ceil() {
        return newInstance(FastMath.ceil(value));
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> floor() {
        return newInstance(FastMath.floor(value));
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> rint() {
        return newInstance(FastMath.rint(value));
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> sign() {
        return newInstance(FastMath.sign(value));
    }

    /**
     * Returns the instance with the sign of the argument.
     * A NaN {@code sign} argument is treated as positive.
     *
     * @param sign the sign for the returned value
     * @return the instance with the same sign as the {@code sign} argument
     */
    public FieldGradient<T> copySign(final T sign) {
        long m = Double.doubleToLongBits(value.getReal());
        long s = Double.doubleToLongBits(sign.getReal());
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
            return this;
        }
        return negate(); // flip sign
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> copySign(final FieldGradient<T> sign) {
        long m = Double.doubleToLongBits(value.getReal());
        long s = Double.doubleToLongBits(sign.value.getReal());
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
            return this;
        }
        return negate(); // flip sign
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> copySign(final double sign) {
        long m = Double.doubleToLongBits(value.getReal());
        long s = Double.doubleToLongBits(sign);
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
            return this;
        }
        return negate(); // flip sign
    }

    /** {@inheritDoc} */
    @Override
    public int getExponent() {
        return FastMath.getExponent(value.getReal());
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> scalb(final int n) {
        final FieldGradient<T> result = newInstance(FastMath.scalb(value, n));
        for (int i = 0; i < grad.length; ++i) {
            result.grad[i] = FastMath.scalb(grad[i], n);
        }
        return result;
    }

    /** {@inheritDoc}
     * <p>
     * The {@code ulp} function is a step function, hence all its derivatives are 0.
     * </p>
     * @since 2.0
     */
    @Override
    public FieldGradient<T> ulp() {
        return newInstance(FastMath.ulp(value));
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> hypot(final FieldGradient<T> y) {

        if (Double.isInfinite(value.getReal()) || Double.isInfinite(y.value.getReal())) {
            return newInstance(Double.POSITIVE_INFINITY);
        } else if (Double.isNaN(value.getReal()) || Double.isNaN(y.value.getReal())) {
            return newInstance(Double.NaN);
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
                final FieldGradient<T> scaledX = scalb(-middleExp);
                final FieldGradient<T> scaledY = y.scalb(-middleExp);

                // compute scaled hypotenuse
                final FieldGradient<T> scaledH =
                        scaledX.multiply(scaledX).add(scaledY.multiply(scaledY)).sqrt();

                // remove scaling
                return scaledH.scalb(middleExp);

            }

        }
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> reciprocal() {
        final T inv1  = value.reciprocal();
        final T mInv2 = inv1.multiply(inv1).negate();
        final FieldGradient<T> result = newInstance(inv1);
        for (int i = 0; i < grad.length; ++i) {
            result.grad[i] = grad[i].multiply(mInv2);
        }
        return result;
    }

    /** Compute composition of the instance by a function.
     * @param g0 value of the function at the current point (i.e. at {@code g(getValue())})
     * @param g1 first derivative of the function at the current point (i.e. at {@code g'(getValue())})
     * @return g(this)
     */
    public FieldGradient<T> compose(final T g0, final T g1) {
        final FieldGradient<T> result = newInstance(g0);
        for (int i = 0; i < grad.length; ++i) {
            result.grad[i] = g1.multiply(grad[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> sqrt() {
        final T s = FastMath.sqrt(value);
        return compose(s, s.add(s).reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> cbrt() {
        final T c = FastMath.cbrt(value);
        return compose(c, c.multiply(c).multiply(3).reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> rootN(final int n) {
        if (n == 2) {
            return sqrt();
        } else if (n == 3) {
            return cbrt();
        } else {
            final T r = FastMath.pow(value, 1.0 / n);
            return compose(r, FastMath.pow(r, n - 1).multiply(n).reciprocal());
        }
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradientField<T> getField() {
        return FieldGradientField.getField(getValueField(), getFreeParameters());
    }

    /** Compute a<sup>x</sup> where a is a double and x a {@link FieldGradient}
     * @param a number to exponentiate
     * @param x power to apply
     * @param <T> the type of the function parameters and value
     * @return a<sup>x</sup>
     */
    public static <T extends CalculusFieldElement<T>> FieldGradient<T> pow(final double a, final FieldGradient<T> x) {
        if (a == 0) {
            return x.getField().getZero();
        } else {
            final T aX    = FastMath.pow(x.value.newInstance(a), x.value);
            final T aXlnA = aX.multiply(FastMath.log(a));
            final FieldGradient<T> result = x.newInstance(aX);
            for (int i = 0; i < x.grad.length; ++i) {
                result.grad[i] =  aXlnA.multiply(x.grad[i]);
            }
            return result;
        }
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> pow(final double p) {
        if (p == 0) {
            return getField().getOne();
        } else {
            final T f0Pm1 = FastMath.pow(value, p - 1);
            return compose(f0Pm1.multiply(value), f0Pm1.multiply(p));
        }
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> pow(final int n) {
        if (n == 0) {
            return getField().getOne();
        } else {
            final T f0Nm1 = FastMath.pow(value, n - 1);
            return compose(f0Nm1.multiply(value), f0Nm1.multiply(n));
        }
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> pow(final FieldGradient<T> e) {
        return log().multiply(e).exp();
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> exp() {
        final T exp = FastMath.exp(value);
        return compose(exp, exp);
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> expm1() {
        final T exp   = FastMath.exp(value);
        final T expM1 = FastMath.expm1(value);
        return compose(expM1, exp);
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> log() {
        return compose(FastMath.log(value), value.reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> log1p() {
        return compose(FastMath.log1p(value), value.add(1).reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> log10() {
        return compose(FastMath.log10(value), value.multiply(FastMath.log(10.0)).reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> cos() {
        final FieldSinCos<T> sinCos = FastMath.sinCos(value);
        return compose(sinCos.cos(), sinCos.sin().negate());
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> sin() {
        final FieldSinCos<T> sinCos = FastMath.sinCos(value);
        return compose(sinCos.sin(), sinCos.cos());
    }

    /** {@inheritDoc} */
    @Override
    public FieldSinCos<FieldGradient<T>> sinCos() {
        final FieldSinCos<T> sinCos = FastMath.sinCos(value);
        final FieldGradient<T> sin = newInstance(sinCos.sin());
        final FieldGradient<T> cos = newInstance(sinCos.cos());
        final T mSin = sinCos.sin().negate();
        for (int i = 0; i < grad.length; ++i) {
            sin.grad[i] =  grad[i].multiply(sinCos.cos());
            cos.grad[i] =  grad[i].multiply(mSin);
        }
        return new FieldSinCos<>(sin, cos);
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> tan() {
        final T tan = FastMath.tan(value);
        return compose(tan, tan.multiply(tan).add(1));
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> acos() {
        return compose(FastMath.acos(value), value.multiply(value).negate().add(1).sqrt().reciprocal().negate());
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> asin() {
        return compose(FastMath.asin(value), value.multiply(value).negate().add(1).sqrt().reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> atan() {
        return compose(FastMath.atan(value), value.multiply(value).add(1).reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> atan2(final FieldGradient<T> x) {
        final T inv = value.multiply(value).add(x.value.multiply(x.value)).reciprocal();
        final FieldGradient<T> result = newInstance(FastMath.atan2(value, x.value));
        final T xValueInv = x.value.multiply(inv);
        final T mValueInv = value.negate().multiply(inv);
        for (int i = 0; i < grad.length; ++i) {
            result.grad[i] = xValueInv.multiply(grad[i]).add(x.grad[i].multiply(mValueInv));
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> cosh() {
        return compose(FastMath.cosh(value), FastMath.sinh(value));
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> sinh() {
        return compose(FastMath.sinh(value), FastMath.cosh(value));
    }

    /** {@inheritDoc} */
    @Override
    public FieldSinhCosh<FieldGradient<T>> sinhCosh() {
        final FieldSinhCosh<T> sinhCosh = FastMath.sinhCosh(value);
        final FieldGradient<T> sinh = newInstance(sinhCosh.sinh());
        final FieldGradient<T> cosh = newInstance(sinhCosh.cosh());
        for (int i = 0; i < grad.length; ++i) {
            sinh.grad[i] = grad[i].multiply(sinhCosh.cosh());
            cosh.grad[i] = grad[i].multiply(sinhCosh.sinh());
        }
        return new FieldSinhCosh<>(sinh, cosh);
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> tanh() {
        final T tanh = FastMath.tanh(value);
        return compose(tanh, tanh.multiply(tanh).negate().add(1));
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> acosh() {
        return compose(FastMath.acosh(value), value.multiply(value).subtract(1).sqrt().reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> asinh() {
        return compose(FastMath.asinh(value), value.multiply(value).add(1).sqrt().reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> atanh() {
        return compose(FastMath.atanh(value), value.multiply(value).negate().add(1).reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> toDegrees() {
        final FieldGradient<T> result = newInstance(FastMath.toDegrees(value));
        for (int i = 0; i < grad.length; ++i) {
            result.grad[i] = FastMath.toDegrees(grad[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> toRadians() {
        final FieldGradient<T> result = newInstance(FastMath.toRadians(value));
        for (int i = 0; i < grad.length; ++i) {
            result.grad[i] = FastMath.toRadians(grad[i]);
        }
        return result;
    }

    /** Evaluate Taylor expansion of a gradient.
     * @param delta parameters offsets (&Delta;x, &Delta;y, ...)
     * @return value of the Taylor expansion at x + &Delta;x, y + &Delta;y, ...
     */
    public T taylor(final double... delta) {
        T result = value;
        for (int i = 0; i < grad.length; ++i) {
            result = result.add(grad[i].multiply(delta[i]));
        }
        return result;
    }

    /** Evaluate Taylor expansion of a gradient.
     * @param delta parameters offsets (&Delta;x, &Delta;y, ...)
     * @return value of the Taylor expansion at x + &Delta;x, y + &Delta;y, ...
     */
    public T taylor(@SuppressWarnings("unchecked") final T... delta) {
        T result = value;
        for (int i = 0; i < grad.length; ++i) {
            result = result.add(grad[i].multiply(delta[i]));
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> linearCombination(final FieldGradient<T>[] a, final FieldGradient<T>[] b) {

        // extract values and first derivatives
        final Field<T> field = a[0].value.getField();
        final int n  = a.length;
        final T[] a0 = MathArrays.buildArray(field, n);
        final T[] b0 = MathArrays.buildArray(field, n);
        final T[] a1 = MathArrays.buildArray(field, 2 * n);
        final T[] b1 = MathArrays.buildArray(field, 2 * n);
        for (int i = 0; i < n; ++i) {
            final FieldGradient<T> ai = a[i];
            final FieldGradient<T> bi = b[i];
            a0[i]         = ai.value;
            b0[i]         = bi.value;
            a1[2 * i]     = ai.value;
            b1[2 * i + 1] = bi.value;
        }

        final FieldGradient<T> result = newInstance(a[0].value.linearCombination(a0, b0));
        for (int k = 0; k < grad.length; ++k) {
            for (int i = 0; i < n; ++i) {
                a1[2 * i + 1] = a[i].grad[k];
                b1[2 * i]     = b[i].grad[k];
            }
            result.grad[k] = a[0].value.linearCombination(a1, b1);
        }
        return result;

    }

    /**
     * Compute a linear combination.
     * @param a Factors.
     * @param b Factors.
     * @return <code>&Sigma;<sub>i</sub> a<sub>i</sub> b<sub>i</sub></code>.
     * @throws MathIllegalArgumentException if arrays dimensions don't match
     */
    public FieldGradient<T> linearCombination(final T[] a, final FieldGradient<T>[] b) {

        // extract values and first derivatives
        final Field<T> field = b[0].value.getField();
        final int      n  = b.length;
        final T[] b0 = MathArrays.buildArray(field, n);
        final T[] b1 = MathArrays.buildArray(field, n);
        for (int i = 0; i < n; ++i) {
            b0[i] = b[i].value;
        }

        final FieldGradient<T> result = newInstance(b[0].value.linearCombination(a, b0));
        for (int k = 0; k < grad.length; ++k) {
            for (int i = 0; i < n; ++i) {
                b1[i] = b[i].grad[k];
            }
            result.grad[k] = b[0].value.linearCombination(a, b1);
        }
        return result;

    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> linearCombination(final double[] a, final FieldGradient<T>[] b) {

        // extract values and first derivatives
        final Field<T> field = b[0].value.getField();
        final int      n  = b.length;
        final T[] b0 = MathArrays.buildArray(field, n);
        final T[] b1 = MathArrays.buildArray(field, n);
        for (int i = 0; i < n; ++i) {
            b0[i] = b[i].value;
        }

        final FieldGradient<T> result = newInstance(b[0].value.linearCombination(a, b0));
        for (int k = 0; k < grad.length; ++k) {
            for (int i = 0; i < n; ++i) {
                b1[i] = b[i].grad[k];
            }
            result.grad[k] = b[0].value.linearCombination(a, b1);
        }
        return result;

    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> linearCombination(final FieldGradient<T> a1, final FieldGradient<T> b1,
                                              final FieldGradient<T> a2, final FieldGradient<T> b2) {
        final FieldGradient<T> result = newInstance(a1.value.linearCombination(a1.value, b1.value,
                                                                               a2.value, b2.value));
        for (int i = 0; i < b1.grad.length; ++i) {
            result.grad[i] = a1.value.linearCombination(a1.value,       b1.grad[i],
                                                            a1.grad[i], b1.value,
                                                            a2.value,       b2.grad[i],
                                                            a2.grad[i], b2.value);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> linearCombination(final double a1, final FieldGradient<T> b1,
                                              final double a2, final FieldGradient<T> b2) {
        final FieldGradient<T> result = newInstance(b1.value.linearCombination(a1, b1.value,
                                                                               a2, b2.value));
        for (int i = 0; i < b1.grad.length; ++i) {
            result.grad[i] = b1.value.linearCombination(a1, b1.grad[i],
                                                            a2, b2.grad[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> linearCombination(final FieldGradient<T> a1, final FieldGradient<T> b1,
                                              final FieldGradient<T> a2, final FieldGradient<T> b2,
                                              final FieldGradient<T> a3, final FieldGradient<T> b3) {
        final Field<T> field = a1.value.getField();
        final T[] a = MathArrays.buildArray(field, 6);
        final T[] b = MathArrays.buildArray(field, 6);
        a[0] = a1.value;
        a[2] = a2.value;
        a[4] = a3.value;
        b[1] = b1.value;
        b[3] = b2.value;
        b[5] = b3.value;
        final FieldGradient<T> result = newInstance(a1.value.linearCombination(a1.value, b1.value,
                                                                               a2.value, b2.value,
                                                                               a3.value, b3.value));
        for (int i = 0; i < b1.grad.length; ++i) {
            a[1] = a1.grad[i];
            a[3] = a2.grad[i];
            a[5] = a3.grad[i];
            b[0] = b1.grad[i];
            b[2] = b2.grad[i];
            b[4] = b3.grad[i];
            result.grad[i] = a1.value.linearCombination(a, b);
        }
        return result;
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
     * @see #linearCombination(double, FieldGradient, double, FieldGradient)
     * @see #linearCombination(double, FieldGradient, double, FieldGradient, double, FieldGradient, double, FieldGradient)
     * @exception MathIllegalArgumentException if number of free parameters or orders are inconsistent
     */
    public FieldGradient<T> linearCombination(final T a1, final FieldGradient<T> b1,
                                              final T a2, final FieldGradient<T> b2,
                                              final T a3, final FieldGradient<T> b3) {
        final FieldGradient<T> result = newInstance(b1.value.linearCombination(a1, b1.value,
                                                                               a2, b2.value,
                                                                               a3, b3.value));
        for (int i = 0; i < b1.grad.length; ++i) {
            result.grad[i] = b1.value.linearCombination(a1, b1.grad[i],
                                                        a2, b2.grad[i],
                                                        a3, b3.grad[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> linearCombination(final double a1, final FieldGradient<T> b1,
                                              final double a2, final FieldGradient<T> b2,
                                              final double a3, final FieldGradient<T> b3) {
        final FieldGradient<T> result = newInstance(b1.value.linearCombination(a1, b1.value,
                                                                               a2, b2.value,
                                                                               a3, b3.value));
        for (int i = 0; i < b1.grad.length; ++i) {
            result.grad[i] = b1.value.linearCombination(a1, b1.grad[i],
                                                            a2, b2.grad[i],
                                                            a3, b3.grad[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> linearCombination(final FieldGradient<T> a1, final FieldGradient<T> b1,
                                              final FieldGradient<T> a2, final FieldGradient<T> b2,
                                              final FieldGradient<T> a3, final FieldGradient<T> b3,
                                              final FieldGradient<T> a4, final FieldGradient<T> b4) {
        final Field<T> field = a1.value.getField();
        final T[] a = MathArrays.buildArray(field, 8);
        final T[] b = MathArrays.buildArray(field, 8);
        a[0] = a1.value;
        a[2] = a2.value;
        a[4] = a3.value;
        a[6] = a4.value;
        b[1] = b1.value;
        b[3] = b2.value;
        b[5] = b3.value;
        b[7] = b4.value;
        final FieldGradient<T> result = newInstance(a1.value.linearCombination(a1.value, b1.value,
                                                                               a2.value, b2.value,
                                                                               a3.value, b3.value,
                                                                               a4.value, b4.value));
        for (int i = 0; i < b1.grad.length; ++i) {
            a[1] = a1.grad[i];
            a[3] = a2.grad[i];
            a[5] = a3.grad[i];
            a[7] = a4.grad[i];
            b[0] = b1.grad[i];
            b[2] = b2.grad[i];
            b[4] = b3.grad[i];
            b[6] = b4.grad[i];
            result.grad[i] = a1.value.linearCombination(a, b);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> linearCombination(final double a1, final FieldGradient<T> b1,
                                              final double a2, final FieldGradient<T> b2,
                                              final double a3, final FieldGradient<T> b3,
                                              final double a4, final FieldGradient<T> b4) {
        final FieldGradient<T> result = newInstance(b1.value.linearCombination(a1, b1.value,
                                                                               a2, b2.value,
                                                                               a3, b3.value,
                                                                               a4, b4.value));
        for (int i = 0; i < b1.grad.length; ++i) {
            result.grad[i] = b1.value.linearCombination(a1, b1.grad[i],
                                                            a2, b2.grad[i],
                                                            a3, b3.grad[i],
                                                            a4, b4.grad[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldGradient<T> getPi() {
        return new FieldGradient<>(getValueField().getZero().getPi(), getFreeParameters());
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

        if (other instanceof FieldGradient) {
            @SuppressWarnings("unchecked")
            final FieldGradient<T> rhs = (FieldGradient<T>) other;
            if (!value.equals(rhs.value) || grad.length != rhs.grad.length) {
                return false;
            }
            for (int i = 0; i < grad.length; ++i) {
                if (!grad[i].equals(rhs.grad[i])) {
                    return false;
                }
            }
            return true;
        }

        return false;

    }

    /** Get a hashCode for the univariate derivative.
     * @return a hash code value for this object
     */
    @Override
    public int hashCode() {
        return 129 + 7 *value.hashCode() - 15 * Arrays.hashCode(grad);
    }

}
