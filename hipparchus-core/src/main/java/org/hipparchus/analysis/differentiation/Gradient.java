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

import java.io.Serializable;
import java.util.Arrays;

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
 * with {@link DerivativeStructure#getOrder() derivation order} limited to one.
 * It should have less overhead than {@link DerivativeStructure} in its domain.</p>
 * <p>This class is an implementation of Rall's numbers. Rall's numbers are an
 * extension to the real numbers used throughout mathematical expressions; they hold
 * the derivative together with the value of a function.</p>
 * <p>{@link Gradient} instances can be used directly thanks to
 * the arithmetic operators to the mathematical functions provided as
 * methods by this class (+, -, *, /, %, sin, cos ...).</p>
 * <p>Implementing complex expressions by hand using {@link Derivative}-based
 * classes (or in fact any {@link org.hipparchus.CalculusFieldElement} class) is
 * a tedious and error-prone task but has the advantage of not requiring users
 * to compute the derivatives by themselves and allowing to switch for one
 * derivative implementation to another as they all share the same filed API.</p>
 * <p>Instances of this class are guaranteed to be immutable.</p>
 * @see DerivativeStructure
 * @see UnivariateDerivative1
 * @see UnivariateDerivative2
 * @see FieldDerivativeStructure
 * @see FieldUnivariateDerivative1
 * @see FieldUnivariateDerivative2
 * @see FieldGradient
 * @since 1.7
 */
public class Gradient implements Derivative1<Gradient>, Serializable {

    /** Serializable UID. */
    private static final long serialVersionUID = 20200520L;

    /** Value of the function. */
    private final double value;

    /** Gradient of the function. */
    private final double[] grad;

    /** Build an instance with values and uninitialized derivatives array.
     * @param value value of the function
     * @param freeParameters number of free parameters
     */
    private Gradient(final double value, int freeParameters) {
        this.value = value;
        this.grad  = new double[freeParameters];
    }

    /** Build an instance with value and derivatives array, used for performance internally (no copy).
     * @param value value of the function
     * @param gradient derivatives
     * @since 3.1
     */
    private Gradient(final double[] gradient, final double value) {
        this.value = value;
        this.grad  = gradient;
    }

    /** Build an instance with values and derivative.
     * @param value value of the function
     * @param gradient gradient of the function
     */
    public Gradient(final double value, final double... gradient) {
        this(value, gradient.length);
        System.arraycopy(gradient, 0, grad, 0, grad.length);
    }

    /** Build an instance from a {@link DerivativeStructure}.
     * @param ds derivative structure
     * @exception MathIllegalArgumentException if {@code ds} order
     * is not 1
     */
    public Gradient(final DerivativeStructure ds) throws MathIllegalArgumentException {
        this(ds.getValue(), ds.getFreeParameters());
        MathUtils.checkDimension(ds.getOrder(), 1);
        System.arraycopy(ds.getAllDerivatives(), 1, grad, 0, grad.length);
    }

    /** Build an instance corresponding to a constant value.
     * @param freeParameters number of free parameters (i.e. dimension of the gradient)
     * @param value constant value of the function
     * @return a {@code Gradient} with a constant value and all derivatives set to 0.0
     */
    public static Gradient constant(final int freeParameters, final double value) {
        return new Gradient(value, freeParameters);
    }

    /** Build a {@code Gradient} representing a variable.
     * <p>Instances built using this method are considered
     * to be the free variables with respect to which differentials
     * are computed. As such, their differential with respect to
     * themselves is +1.</p>
     * @param freeParameters number of free parameters (i.e. dimension of the gradient)
     * @param index index of the variable (from 0 to {@link #getFreeParameters() getFreeParameters()} - 1)
     * @param value value of the variable
     * @return a {@code Gradient} with a constant value and all derivatives set to 0.0 except the
     * one at {@code index} which will be set to 1.0
     */
    public static Gradient variable(final int freeParameters, final int index, final double value) {
        final Gradient g = new Gradient(value, freeParameters);
        g.grad[index] = 1.0;
        return g;
    }

    /** {@inheritDoc} */
    @Override
    public Gradient newInstance(final double c) {
        return new Gradient(c, new double[grad.length]);
    }

    /** {@inheritDoc} */
    @Override
    public Gradient withValue(final double v) {
        return new Gradient(v, grad);
    }

    /** {@inheritDoc} */
    @Override
    public Gradient getAddendum() {
        return new Gradient(0, grad);
    }

    /** Get the value part of the function.
     * @return value part of the value of the function
     */
    @Override
    public double getValue() {
        return value;
    }

    /** Get the gradient part of the function.
     * @return gradient part of the value of the function
     * @see #getPartialDerivative(int)
     */
    public double[] getGradient() {
        return grad.clone();
    }

    /** {@inheritDoc} */
    @Override
    public int getFreeParameters() {
        return grad.length;
    }

    /** {@inheritDoc} */
    @Override
    public double getPartialDerivative(final int ... orders)
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
    public double getPartialDerivative(final int n) throws MathIllegalArgumentException {
        if (n < 0 || n >= grad.length) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.OUT_OF_RANGE_SIMPLE, n, 0, grad.length - 1);
        }
        return grad[n];
    }

    /** Convert the instance to a {@link DerivativeStructure}.
     * @return derivative structure with same value and derivative as the instance
     */
    public DerivativeStructure toDerivativeStructure() {
        final double[] derivatives = new double[1 + grad.length];
        derivatives[0] = value;
        System.arraycopy(grad, 0, derivatives, 1, grad.length);
        return getField().getConversionFactory().build(derivatives);
    }

    /** {@inheritDoc} */
    @Override
    public Gradient add(final Gradient a) {
        final double[] gradient = new double[grad.length];
        for (int i = 0; i < grad.length; ++i) {
            gradient[i] = grad[i] + a.grad[i];
        }
        return new Gradient(gradient, value + a.value);
    }

    /** {@inheritDoc} */
    @Override
    public Gradient subtract(final Gradient a) {
        final double[] gradient = new double[grad.length];
        for (int i = 0; i < grad.length; ++i) {
            gradient[i] = grad[i] - a.grad[i];
        }
        return new Gradient(gradient, value - a.value);
    }

    /** {@inheritDoc} */
    @Override
    public Gradient multiply(final int n) {
        final double[] gradient = new double[grad.length];
        for (int i = 0; i < grad.length; ++i) {
            gradient[i] = grad[i] * n;
        }
        return new Gradient(gradient, value * n);
    }

    /** {@inheritDoc} */
    @Override
    public Gradient multiply(final double a) {
        final double[] gradient = new double[grad.length];
        for (int i = 0; i < grad.length; ++i) {
            gradient[i] = grad[i] * a;
        }
        return new Gradient(gradient, value * a);
    }

    /** {@inheritDoc} */
    @Override
    public Gradient multiply(final Gradient a) {
        final double[] gradient = new double[grad.length];
        for (int i = 0; i < grad.length; ++i) {
            gradient[i] = grad[i] * a.value + value * a.grad[i];
        }
        return new Gradient(gradient, value * a.value);
    }

    /** {@inheritDoc} */
    @Override
    public Gradient divide(final double a) {
        final double[] gradient = new double[grad.length];
        for (int i = 0; i < grad.length; ++i) {
            gradient[i] = grad[i] / a;
        }
        return new Gradient(gradient, value / a);
    }

    /** {@inheritDoc} */
    @Override
    public Gradient divide(final Gradient a) {
        final double inv1 = 1.0 / a.value;
        final double inv2 = inv1 * inv1;
        final double[] gradient = new double[grad.length];
        for (int i = 0; i < grad.length; ++i) {
            gradient[i] = (grad[i] * a.value - value * a.grad[i]) * inv2;
        }
        return new Gradient(gradient, value * inv1);
    }

    /** {@inheritDoc} */
    @Override
    public Gradient remainder(final Gradient a) {

        // compute k such that lhs % rhs = lhs - k rhs
        final double rem = FastMath.IEEEremainder(value, a.value);
        final double k   = FastMath.rint((value - rem) / a.value);

        final double[] gradient = new double[grad.length];
        for (int i = 0; i < grad.length; ++i) {
            gradient[i] = grad[i] - k * a.grad[i];
        }
        return new Gradient(gradient, rem);

    }

    /** {@inheritDoc} */
    @Override
    public Gradient negate() {
        final double[] gradient = new double[grad.length];
        for (int i = 0; i < grad.length; ++i) {
            gradient[i] = -grad[i];
        }
        return new Gradient(gradient, -value);
    }

    /** {@inheritDoc} */
    @Override
    public Gradient abs() {
        if (Double.doubleToLongBits(value) < 0) {
            // we use the bits representation to also handle -0.0
            return negate();
        } else {
            return this;
        }
    }

    /** {@inheritDoc} */
    @Override
    public Gradient copySign(final Gradient sign) {
        long m = Double.doubleToLongBits(value);
        long s = Double.doubleToLongBits(sign.value);
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
            return this;
        }
        return negate(); // flip sign
    }

    /** {@inheritDoc} */
    @Override
    public Gradient copySign(final double sign) {
        long m = Double.doubleToLongBits(value);
        long s = Double.doubleToLongBits(sign);
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
            return this;
        }
        return negate(); // flip sign
    }

    /** {@inheritDoc} */
    @Override
    public Gradient scalb(final int n) {
        final double[] gradient = new double[grad.length];
        for (int i = 0; i < grad.length; ++i) {
            gradient[i] = FastMath.scalb(grad[i], n);
        }
        return new Gradient(gradient, FastMath.scalb(value, n));
    }

    /** {@inheritDoc} */
    @Override
    public Gradient hypot(final Gradient y) {

        if (Double.isInfinite(value) || Double.isInfinite(y.value)) {
            return newInstance(Double.POSITIVE_INFINITY);
        } else if (Double.isNaN(value) || Double.isNaN(y.value)) {
            return newInstance(Double.NaN);
        } else {

            final int expX = getExponent();
            final int expY = y.getExponent();
            if (expX > expY + 27) {
                // y is negligible with respect to x
                return abs();
            } else if (expY > expX + 27) {
                // x is negligible with respect to y
                return y.abs();
            } else {

                // find an intermediate scale to avoid both overflow and underflow
                final int middleExp = (expX + expY) / 2;

                // scale parameters without losing precision
                final Gradient scaledX = scalb(-middleExp);
                final Gradient scaledY = y.scalb(-middleExp);

                // compute scaled hypotenuse
                final Gradient scaledH =
                        scaledX.multiply(scaledX).add(scaledY.multiply(scaledY)).sqrt();

                // remove scaling
                return scaledH.scalb(middleExp);

            }

        }
    }

    /** {@inheritDoc} */
    @Override
    public Gradient compose(final double... f) {
        MathUtils.checkDimension(f.length, getOrder() + 1);
        return compose(f[0], f[1]);
    }

    /** {@inheritDoc} */
    @Override
    public Gradient compose(final double f0, final double f1) {
        final double[] gradient = new double[grad.length];
        for (int i = 0; i < grad.length; ++i) {
            gradient[i] = f1 * grad[i];
        }
        return new Gradient(gradient, f0);
    }

    /** {@inheritDoc} */
    @Override
    public GradientField getField() {
        return GradientField.getField(getFreeParameters());
    }

    /** Compute a<sup>x</sup> where a is a double and x a {@link Gradient}
     * @param a number to exponentiate
     * @param x power to apply
     * @return a<sup>x</sup>
     */
    public static Gradient pow(final double a, final Gradient x) {
        if (a == 0) {
            return x.getField().getZero();
        } else {
            final double aX = FastMath.pow(a, x.value);
            final double aXlnA = aX * FastMath.log(a);
            final double[] gradient = new double[x.getFreeParameters()];
            for (int i = 0; i < gradient.length; ++i) {
                gradient[i] =  aXlnA * x.grad[i];
            }
            return new Gradient(gradient, aX);
        }
    }

    /** {@inheritDoc} */
    @Override
    public Gradient pow(final double p) {
        if (p == 0) {
            return getField().getOne();
        } else {
            final double valuePm1 = FastMath.pow(value, p - 1);
            return compose(valuePm1 * value, p * valuePm1);
        }
    }

    /** {@inheritDoc} */
    @Override
    public Gradient pow(final int n) {
        if (n == 0) {
            return getField().getOne();
        } else {
            final double valueNm1 = FastMath.pow(value, n - 1);
            return compose(valueNm1 * value, n * valueNm1);
        }
    }

    /** {@inheritDoc} */
    @Override
    public FieldSinCos<Gradient> sinCos() {
        final SinCos sinCos = FastMath.sinCos(value);
        final double[] gradSin = new double[grad.length];
        final double[] gradCos = new double[grad.length];
        for (int i = 0; i < grad.length; ++i) {
            gradSin[i] =  grad[i] * sinCos.cos();
            gradCos[i] =  -grad[i] * sinCos.sin();
        }
        final Gradient sin = new Gradient(gradSin, sinCos.sin());
        final Gradient cos = new Gradient(gradCos, sinCos.cos());
        return new FieldSinCos<>(sin, cos);
    }

    /** {@inheritDoc} */
    @Override
    public Gradient atan2(final Gradient x) {
        final double inv = 1.0 / (value * value + x.value * x.value);
        final double[] gradient = new double[grad.length];
        for (int i = 0; i < grad.length; ++i) {
            gradient[i] = (x.value * grad[i] - x.grad[i] * value) * inv;
        }
        return new Gradient(gradient, FastMath.atan2(value, x.value));
    }

    /** {@inheritDoc} */
    @Override
    public FieldSinhCosh<Gradient> sinhCosh() {
        final SinhCosh sinhCosh = FastMath.sinhCosh(value);
        final double[] gradSinh = new double[grad.length];
        final double[] gradCosh = new double[grad.length];
        for (int i = 0; i < grad.length; ++i) {
            gradSinh[i] = grad[i] * sinhCosh.cosh();
            gradCosh[i] = grad[i] * sinhCosh.sinh();
        }
        final Gradient sinh = new Gradient(gradSinh, sinhCosh.sinh());
        final Gradient cosh = new Gradient(gradCosh, sinhCosh.cosh());
        return new FieldSinhCosh<>(sinh, cosh);
    }

    /** {@inheritDoc} */
    @Override
    public Gradient toDegrees() {
        final double[] gradient = new double[grad.length];
        for (int i = 0; i < grad.length; ++i) {
            gradient[i] = FastMath.toDegrees(grad[i]);
        }
        return new Gradient(gradient, FastMath.toDegrees(value));
    }

    /** {@inheritDoc} */
    @Override
    public Gradient toRadians() {
        final double[] gradient = new double[grad.length];
        for (int i = 0; i < grad.length; ++i) {
            gradient[i] = FastMath.toRadians(grad[i]);
        }
        return new Gradient(gradient, FastMath.toRadians(value));
    }

    /** Evaluate Taylor expansion a derivative structure.
     * @param delta parameters offsets (&Delta;x, &Delta;y, ...)
     * @return value of the Taylor expansion at x + &Delta;x, y + &Delta;y, ...
     */
    public double taylor(final double... delta) {
        double result = value;
        for (int i = 0; i < grad.length; ++i) {
            result += delta[i] * grad[i];
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Gradient linearCombination(final Gradient[] a, final Gradient[] b) {

        // extract values and first derivatives
        final int      n  = a.length;
        final double[] a0 = new double[n];
        final double[] b0 = new double[n];
        final double[] a1 = new double[2 * n];
        final double[] b1 = new double[2 * n];
        for (int i = 0; i < n; ++i) {
            final Gradient ai = a[i];
            final Gradient bi = b[i];
            a0[i]         = ai.value;
            b0[i]         = bi.value;
            a1[2 * i]     = ai.value;
            b1[2 * i + 1] = bi.value;
        }

        final Gradient result = newInstance(MathArrays.linearCombination(a0, b0));
        for (int k = 0; k < grad.length; ++k) {
            for (int i = 0; i < n; ++i) {
                a1[2 * i + 1] = a[i].grad[k];
                b1[2 * i]     = b[i].grad[k];
            }
            result.grad[k] = MathArrays.linearCombination(a1, b1);
        }
        return result;

    }

    /** {@inheritDoc} */
    @Override
    public Gradient linearCombination(final double[] a, final Gradient[] b) {

        // extract values and first derivatives
        final int      n  = b.length;
        final double[] b0 = new double[n];
        final double[] b1 = new double[n];
        for (int i = 0; i < n; ++i) {
            b0[i] = b[i].value;
        }

        final Gradient result = newInstance(MathArrays.linearCombination(a, b0));
        for (int k = 0; k < grad.length; ++k) {
            for (int i = 0; i < n; ++i) {
                b1[i] = b[i].grad[k];
            }
            result.grad[k] = MathArrays.linearCombination(a, b1);
        }
        return result;

    }

    /** {@inheritDoc} */
    @Override
    public Gradient linearCombination(final Gradient a1, final Gradient b1,
                                      final Gradient a2, final Gradient b2) {
        final Gradient result = newInstance(MathArrays.linearCombination(a1.value, b1.value,
                a2.value, b2.value));
        for (int i = 0; i < b1.grad.length; ++i) {
            result.grad[i] = MathArrays.linearCombination(a1.value,       b1.grad[i],
                    a1.grad[i], b1.value,
                    a2.value,       b2.grad[i],
                    a2.grad[i], b2.value);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Gradient linearCombination(final double a1, final Gradient b1,
                                      final double a2, final Gradient b2) {
        final Gradient result = newInstance(MathArrays.linearCombination(a1, b1.value,
                a2, b2.value));
        for (int i = 0; i < b1.grad.length; ++i) {
            result.grad[i] = MathArrays.linearCombination(a1, b1.grad[i],
                    a2, b2.grad[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Gradient linearCombination(final Gradient a1, final Gradient b1,
                                      final Gradient a2, final Gradient b2,
                                      final Gradient a3, final Gradient b3) {
        final double[] a = {
                a1.value, 0, a2.value, 0, a3.value, 0
        };
        final double[] b = {
                0, b1.value, 0, b2.value, 0, b3.value
        };
        final Gradient result = newInstance(MathArrays.linearCombination(a1.value, b1.value,
                a2.value, b2.value,
                a3.value, b3.value));
        for (int i = 0; i < b1.grad.length; ++i) {
            a[1] = a1.grad[i];
            a[3] = a2.grad[i];
            a[5] = a3.grad[i];
            b[0] = b1.grad[i];
            b[2] = b2.grad[i];
            b[4] = b3.grad[i];
            result.grad[i] = MathArrays.linearCombination(a, b);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Gradient linearCombination(final double a1, final Gradient b1,
                                      final double a2, final Gradient b2,
                                      final double a3, final Gradient b3) {
        final Gradient result = newInstance(MathArrays.linearCombination(a1, b1.value,
                a2, b2.value,
                a3, b3.value));
        for (int i = 0; i < b1.grad.length; ++i) {
            result.grad[i] = MathArrays.linearCombination(a1, b1.grad[i],
                    a2, b2.grad[i],
                    a3, b3.grad[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Gradient linearCombination(final Gradient a1, final Gradient b1,
                                      final Gradient a2, final Gradient b2,
                                      final Gradient a3, final Gradient b3,
                                      final Gradient a4, final Gradient b4) {
        final double[] a = {
                a1.value, 0, a2.value, 0, a3.value, 0, a4.value, 0
        };
        final double[] b = {
                0, b1.value, 0, b2.value, 0, b3.value, 0, b4.value
        };
        final Gradient result = newInstance(MathArrays.linearCombination(a1.value, b1.value,
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
            result.grad[i] = MathArrays.linearCombination(a, b);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Gradient linearCombination(final double a1, final Gradient b1,
                                      final double a2, final Gradient b2,
                                      final double a3, final Gradient b3,
                                      final double a4, final Gradient b4) {
        final Gradient result = newInstance(MathArrays.linearCombination(a1, b1.value,
                a2, b2.value,
                a3, b3.value,
                a4, b4.value));
        for (int i = 0; i < b1.grad.length; ++i) {
            result.grad[i] = MathArrays.linearCombination(a1, b1.grad[i],
                    a2, b2.grad[i],
                    a3, b3.grad[i],
                    a4, b4.grad[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public double norm() {
        double sum = FastMath.abs(getValue());
        for (int i = 0; i < getFreeParameters(); i++) {
            sum += FastMath.abs(getPartialDerivative(i));
        }
        return sum;
    }

    /**
     * Add an independent variable to the Taylor expansion.
     * @return object with one more variable
     * @since 4.0
     */
    public Gradient stackVariable() {
        final double[] gradient = new double[this.getFreeParameters() + 1];
        System.arraycopy(this.grad, 0, gradient, 0, this.getFreeParameters());
        return new Gradient(gradient, this.getValue());
    }

    /** {@inheritDoc} */
    @Override
    public boolean isSmall(final Gradient base, final double relativeThreshold) {
        if (grad.length != base.grad.length) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   grad.length, base.grad.length);
        }
        boolean ok = FastMath.abs(value) <= FastMath.abs(base.value * relativeThreshold);
        for (int i = 0; i < grad.length && ok; ++i) {
            ok = FastMath.abs(grad[i]) <= FastMath.abs(base.grad[i] * relativeThreshold);
        }
        return ok;
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

        if (other instanceof Gradient) {
            final Gradient rhs = (Gradient) other;
            return value == rhs.value && MathArrays.equals(grad, rhs.grad);
        }

        return false;

    }

    /** Get a hashCode for the univariate derivative.
     * @return a hash code value for this object
     */
    @Override
    public int hashCode() {
        return 129 + 7 * Double.hashCode(value) - 15 * Arrays.hashCode(grad);
    }

}
