/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */
package org.hipparchus.analysis.differentiation;

import java.io.Serializable;

import org.hipparchus.Field;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.FieldSinCos;
import org.hipparchus.util.FieldSinhCosh;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.MathUtils;

/** Class representing both the value and the differentials of a function.
 * <p>This class is the workhorse of the differentiation package.</p>
 * <p>This class is an implementation of the extension to Rall's
 * numbers described in Dan Kalman's paper <a
 * href="http://www.dankalman.net/AUhome/pdffiles/mmgautodiff.pdf">Doubly
 * Recursive Multivariate Automatic Differentiation</a>, Mathematics Magazine, vol. 75,
 * no. 3, June 2002. Rall's numbers are an extension to the real numbers used
 * throughout mathematical expressions; they hold the derivative together with the
 * value of a function. Dan Kalman's derivative structures hold all partial derivatives
 * up to any specified order, with respect to any number of free parameters. Rall's
 * numbers therefore can be seen as derivative structures for order one derivative and
 * one free parameter, and real numbers can be seen as derivative structures with zero
 * order derivative and no free parameters.</p>
 * <p>{@link DerivativeStructure} instances can be used directly thanks to
 * the arithmetic operators to the mathematical functions provided as
 * methods by this class (+, -, *, /, %, sin, cos ...).</p>
 * <p>Implementing complex expressions by hand using these classes is
 * a tedious and error-prone task but has the advantage of having no limitation
 * on the derivation order despite not requiring users to compute the derivatives by
 * themselves. Implementing complex expression can also be done by developing computation
 * code using standard primitive double values and to use {@link
 * UnivariateFunctionDifferentiator differentiators} to create the {@link
 * DerivativeStructure}-based instances. This method is simpler but may be limited in
 * the accuracy and derivation orders and may be computationally intensive (this is
 * typically the case for {@link FiniteDifferencesDifferentiator finite differences
 * differentiator}.</p>
 * <p>Instances of this class are guaranteed to be immutable.</p>
 * @see DSCompiler
 * @see FieldDerivativeStructure
 */
public class DerivativeStructure implements Derivative<DerivativeStructure>, Serializable {

    /** Serializable UID. */
    private static final long serialVersionUID = 20161220L;

    /** Factory that built the instance. */
    private final DSFactory factory;

    /** Combined array holding all values. */
    private final double[] data;

    /** Build an instance with all values and derivatives set to 0.
     * @param factory factory that built the instance
     * @param data combined array holding all values
     */
    DerivativeStructure(final DSFactory factory, final double[] data) {
        this.factory = factory;
        this.data    = data.clone();
    }

    /** Build an instance with all values and derivatives set to 0.
     * @param factory factory that built the instance
     * @since 1.4
     */
    DerivativeStructure(final DSFactory factory) {
        this.factory = factory;
        this.data    = new double[factory.getCompiler().getSize()];
    }

    /** {@inheritDoc} */
    @Override
    public DerivativeStructure newInstance(final double value) {
        return factory.constant(value);
    }

    /** Get the factory that built the instance.
     * @return factory that built the instance
     */
    public DSFactory getFactory() {
        return factory;
    }

    @Override
    /** {@inheritDoc} */
    public int getFreeParameters() {
        return getFactory().getCompiler().getFreeParameters();
    }

    @Override
    /** {@inheritDoc} */
    public int getOrder() {
        return getFactory().getCompiler().getOrder();
    }

    /** Set a derivative component.
     * <p>
     * This method is package-private (no modifier specified), as it is intended
     * to be used only by Hipparchus classes since it relied on the ordering of
     * derivatives within the class. This allows avoiding checks on the index,
     * for performance reasons.
     * </p>
     * @param index index of the derivative
     * @param value of the derivative to set
     * @since 1.4
     */
    void setDerivativeComponent(final int index, final double value) {
        data[index] = value;
    }

    /** Get a derivative component.
     * <p>
     * This method is package-private (no modifier specified), as it is intended
     * to be used only by Hipparchus classes since it relied on the ordering of
     * derivatives within the class. This allows avoiding checks on the index,
     * for performance reasons.
     * </p>
     * @param index index of the derivative
     * @return value of the derivative
     * @since 2.2
     */
    double getDerivativeComponent(final int index) {
        return data[index];
    }

    /** {@inheritDoc}
     */
    @Override
    public double getReal() {
        return data[0];
    }

    /** Get the value part of the derivative structure.
     * @return value part of the derivative structure
     * @see #getPartialDerivative(int...)
     */
    @Override
    public double getValue() {
        return data[0];
    }

    /** {@inheritDoc} */
    @Override
    public double getPartialDerivative(final int ... orders)
        throws MathIllegalArgumentException {
        return data[getFactory().getCompiler().getPartialDerivativeIndex(orders)];
    }

    /** Get all partial derivatives.
     * @return a fresh copy of partial derivatives, in an array sorted according to
     * {@link DSCompiler#getPartialDerivativeIndex(int...)}
     */
    public double[] getAllDerivatives() {
        return data.clone();
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure add(final double a) {
        final DerivativeStructure ds = factory.build();
        System.arraycopy(data, 0, ds.data, 0, data.length);
        ds.data[0] += a;
        return ds;
    }

    /** {@inheritDoc}
     * @exception MathIllegalArgumentException if number of free parameters
     * or orders do not match
     */
    @Override
    public DerivativeStructure add(final DerivativeStructure a)
        throws MathIllegalArgumentException {
        factory.checkCompatibility(a.factory);
        final DerivativeStructure ds = factory.build();
        factory.getCompiler().add(data, 0, a.data, 0, ds.data, 0);
        return ds;
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure subtract(final double a) {
        return add(-a);
    }

    /** {@inheritDoc}
     * @exception MathIllegalArgumentException if number of free parameters
     * or orders do not match
     */
    @Override
    public DerivativeStructure subtract(final DerivativeStructure a)
        throws MathIllegalArgumentException {
        factory.checkCompatibility(a.factory);
        final DerivativeStructure ds = factory.build();
        factory.getCompiler().subtract(data, 0, a.data, 0, ds.data, 0);
        return ds;
    }

    /** {@inheritDoc} */
    @Override
    public DerivativeStructure multiply(final int n) {
        return multiply((double) n);
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure multiply(final double a) {
        final DerivativeStructure ds = factory.build();
        for (int i = 0; i < ds.data.length; ++i) {
            ds.data[i] = data[i] * a;
        }
        return ds;
    }

    /** {@inheritDoc}
     * @exception MathIllegalArgumentException if number of free parameters
     * or orders do not match
     */
    @Override
    public DerivativeStructure multiply(final DerivativeStructure a)
        throws MathIllegalArgumentException {
        factory.checkCompatibility(a.factory);
        final DerivativeStructure result = factory.build();
        factory.getCompiler().multiply(data, 0, a.data, 0, result.data, 0);
        return result;
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure divide(final double a) {
        final DerivativeStructure ds = factory.build();
        final double inv = 1.0 / a;
        for (int i = 0; i < ds.data.length; ++i) {
            ds.data[i] = data[i] * inv;
        }
        return ds;
    }

    /** {@inheritDoc}
     * @exception MathIllegalArgumentException if number of free parameters
     * or orders do not match
     */
    @Override
    public DerivativeStructure divide(final DerivativeStructure a)
        throws MathIllegalArgumentException {
        factory.checkCompatibility(a.factory);
        final DerivativeStructure result = factory.build();
        factory.getCompiler().divide(data, 0, a.data, 0, result.data, 0);
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public DerivativeStructure remainder(final double a) {
        final DerivativeStructure ds = factory.build();
        System.arraycopy(data, 0, ds.data, 0, data.length);
        ds.data[0] = FastMath.IEEEremainder(ds.data[0], a);
        return ds;
    }

    /** {@inheritDoc}
     * @exception MathIllegalArgumentException if number of free parameters
     * or orders do not match
     */
    @Override
    public DerivativeStructure remainder(final DerivativeStructure a)
        throws MathIllegalArgumentException {
        factory.checkCompatibility(a.factory);
        final DerivativeStructure result = factory.build();
        factory.getCompiler().remainder(data, 0, a.data, 0, result.data, 0);
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public DerivativeStructure negate() {
        final DerivativeStructure ds = factory.build();
        for (int i = 0; i < ds.data.length; ++i) {
            ds.data[i] = -data[i];
        }
        return ds;
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure abs() {
        if (Double.doubleToLongBits(data[0]) < 0) {
            // we use the bits representation to also handle -0.0
            return negate();
        } else {
            return this;
        }
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure ceil() {
        return factory.constant(FastMath.ceil(data[0]));
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure floor() {
        return factory.constant(FastMath.floor(data[0]));
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure rint() {
        return factory.constant(FastMath.rint(data[0]));
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure sign() {
        return factory.constant(FastMath.signum(data[0]));
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure copySign(final DerivativeStructure sign) {
        long m = Double.doubleToLongBits(data[0]);
        long s = Double.doubleToLongBits(sign.data[0]);
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
            return this;
        }
        return negate(); // flip sign
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure copySign(final double sign) {
        long m = Double.doubleToLongBits(data[0]);
        long s = Double.doubleToLongBits(sign);
        if ((m >= 0 && s >= 0) || (m < 0 && s < 0)) { // Sign is currently OK
            return this;
        }
        return negate(); // flip sign
    }

    /**
     * Return the exponent of the instance value, removing the bias.
     * <p>
     * For double numbers of the form 2<sup>x</sup>, the unbiased
     * exponent is exactly x.
     * </p>
     * @return exponent for instance in IEEE754 representation, without bias
     */
    @Override
    public int getExponent() {
        return FastMath.getExponent(data[0]);
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure scalb(final int n) {
        final DerivativeStructure ds = factory.build();
        for (int i = 0; i < ds.data.length; ++i) {
            ds.data[i] = FastMath.scalb(data[i], n);
        }
        return ds;
    }

    /** {@inheritDoc}
     * <p>
     * The {@code ulp} function is a step function, hence all its derivatives are 0.
     * </p>
     * @since 2.0
     */
    @Override
    public DerivativeStructure ulp() {
        final DerivativeStructure ds = factory.build();
        ds.data[0] = FastMath.ulp(data[0]);
        return ds;
    }

    /** {@inheritDoc}
     * @exception MathIllegalArgumentException if number of free parameters
     * or orders do not match
     */
    @Override
    public DerivativeStructure hypot(final DerivativeStructure y)
        throws MathIllegalArgumentException {

        factory.checkCompatibility(y.factory);

        if (Double.isInfinite(data[0]) || Double.isInfinite(y.data[0])) {
            return factory.constant(Double.POSITIVE_INFINITY);
        } else if (Double.isNaN(data[0]) || Double.isNaN(y.data[0])) {
            return factory.constant(Double.NaN);
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
                final DerivativeStructure scaledX = scalb(-middleExp);
                final DerivativeStructure scaledY = y.scalb(-middleExp);

                // compute scaled hypotenuse
                final DerivativeStructure scaledH =
                        scaledX.multiply(scaledX).add(scaledY.multiply(scaledY)).sqrt();

                // remove scaling
                return scaledH.scalb(middleExp);

            }

        }
    }

    /**
     * Returns the hypotenuse of a triangle with sides {@code x} and {@code y}
     * - sqrt(<i>x</i><sup>2</sup>&nbsp;+<i>y</i><sup>2</sup>)
     * avoiding intermediate overflow or underflow.
     *
     * <ul>
     * <li> If either argument is infinite, then the result is positive infinity.</li>
     * <li> else, if either argument is NaN then the result is NaN.</li>
     * </ul>
     *
     * @param x a value
     * @param y a value
     * @return sqrt(<i>x</i><sup>2</sup>&nbsp;+<i>y</i><sup>2</sup>)
     * @exception MathIllegalArgumentException if number of free parameters
     * or orders do not match
     */
    public static DerivativeStructure hypot(final DerivativeStructure x, final DerivativeStructure y)
        throws MathIllegalArgumentException {
        return x.hypot(y);
    }

    /** Compute composition of the instance by a univariate function.
     * @param f array of value and derivatives of the function at
     * the current point (i.e. [f({@link #getValue()}),
     * f'({@link #getValue()}), f''({@link #getValue()})...]).
     * @return f(this)
     * @exception MathIllegalArgumentException if the number of derivatives
     * in the array is not equal to {@link #getOrder() order} + 1
     */
    @Override
    public DerivativeStructure compose(final double ... f)
        throws MathIllegalArgumentException {

        MathUtils.checkDimension(f.length, getOrder() + 1);
        final DerivativeStructure result = factory.build();
        factory.getCompiler().compose(data, 0, f, result.data, 0);
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public DerivativeStructure reciprocal() {
        final DerivativeStructure result = factory.build();
        factory.getCompiler().reciprocal(data, 0, result.data, 0);
        return result;
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure sqrt() {
        final DerivativeStructure result = factory.build();
        factory.getCompiler().sqrt(data, 0, result.data, 0);
        return result;
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure cbrt() {
        return rootN(3);
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure rootN(final int n) {
        final DerivativeStructure result = factory.build();
        factory.getCompiler().rootN(data, 0, n, result.data, 0);
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Field<DerivativeStructure> getField() {
        return factory.getDerivativeField();
    }

    /** Compute a<sup>x</sup> where a is a double and x a {@link DerivativeStructure}
     * @param a number to exponentiate
     * @param x power to apply
     * @return a<sup>x</sup>
     */
    public static DerivativeStructure pow(final double a, final DerivativeStructure x) {
        final DerivativeStructure result = x.factory.build();
        x.factory.getCompiler().pow(a, x.data, 0, result.data, 0);
        return result;
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure pow(final double p) {
        final DerivativeStructure result = factory.build();
        factory.getCompiler().pow(data, 0, p, result.data, 0);
        return result;
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure pow(final int n) {
        final DerivativeStructure result = factory.build();
        factory.getCompiler().pow(data, 0, n, result.data, 0);
        return result;
    }

    /** {@inheritDoc}
     * @exception MathIllegalArgumentException if number of free parameters
     * or orders do not match
     */
    @Override
    public DerivativeStructure pow(final DerivativeStructure e)
        throws MathIllegalArgumentException {
        factory.checkCompatibility(e.factory);
        final DerivativeStructure result = factory.build();
        factory.getCompiler().pow(data, 0, e.data, 0, result.data, 0);
        return result;
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure exp() {
        final DerivativeStructure result = factory.build();
        factory.getCompiler().exp(data, 0, result.data, 0);
        return result;
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure expm1() {
        final DerivativeStructure result = factory.build();
        factory.getCompiler().expm1(data, 0, result.data, 0);
        return result;
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure log() {
        final DerivativeStructure result = factory.build();
        factory.getCompiler().log(data, 0, result.data, 0);
        return result;
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure log1p() {
        final DerivativeStructure result = factory.build();
        factory.getCompiler().log1p(data, 0, result.data, 0);
        return result;
    }

    /** Base 10 logarithm.
     * @return base 10 logarithm of the instance
     */
    @Override
    public DerivativeStructure log10() {
        final DerivativeStructure result = factory.build();
        factory.getCompiler().log10(data, 0, result.data, 0);
        return result;
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure cos() {
        final DerivativeStructure result = factory.build();
        factory.getCompiler().cos(data, 0, result.data, 0);
        return result;
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure sin() {
        final DerivativeStructure result = factory.build();
        factory.getCompiler().sin(data, 0, result.data, 0);
        return result;
    }

    /** {@inheritDoc}
     */
    @Override
    public FieldSinCos<DerivativeStructure> sinCos() {
        final DerivativeStructure sin = factory.build();
        final DerivativeStructure cos = factory.build();
        factory.getCompiler().sinCos(data, 0, sin.data, 0, cos.data, 0);
        return new FieldSinCos<>(sin, cos);
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure tan() {
        final DerivativeStructure result = factory.build();
        factory.getCompiler().tan(data, 0, result.data, 0);
        return result;
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure acos() {
        final DerivativeStructure result = factory.build();
        factory.getCompiler().acos(data, 0, result.data, 0);
        return result;
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure asin() {
        final DerivativeStructure result = factory.build();
        factory.getCompiler().asin(data, 0, result.data, 0);
        return result;
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure atan() {
        final DerivativeStructure result = factory.build();
        factory.getCompiler().atan(data, 0, result.data, 0);
        return result;
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure atan2(final DerivativeStructure x)
        throws MathIllegalArgumentException {
        factory.checkCompatibility(x.factory);
        final DerivativeStructure result = factory.build();
        factory.getCompiler().atan2(data, 0, x.data, 0, result.data, 0);
        return result;
    }

    /** Two arguments arc tangent operation.
     * @param y first argument of the arc tangent
     * @param x second argument of the arc tangent
     * @return atan2(y, x)
     * @exception MathIllegalArgumentException if number of free parameters
     * or orders do not match
     */
    public static DerivativeStructure atan2(final DerivativeStructure y, final DerivativeStructure x)
        throws MathIllegalArgumentException {
        return y.atan2(x);
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure cosh() {
        final DerivativeStructure result = factory.build();
        factory.getCompiler().cosh(data, 0, result.data, 0);
        return result;
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure sinh() {
        final DerivativeStructure result = factory.build();
        factory.getCompiler().sinh(data, 0, result.data, 0);
        return result;
    }

    /** {@inheritDoc}
     */
    @Override
    public FieldSinhCosh<DerivativeStructure> sinhCosh() {
        final DerivativeStructure sinh = factory.build();
        final DerivativeStructure cosh = factory.build();
        factory.getCompiler().sinhCosh(data, 0, sinh.data, 0, cosh.data, 0);
        return new FieldSinhCosh<>(sinh, cosh);
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure tanh() {
        final DerivativeStructure result = factory.build();
        factory.getCompiler().tanh(data, 0, result.data, 0);
        return result;
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure acosh() {
        final DerivativeStructure result = factory.build();
        factory.getCompiler().acosh(data, 0, result.data, 0);
        return result;
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure asinh() {
        final DerivativeStructure result = factory.build();
        factory.getCompiler().asinh(data, 0, result.data, 0);
        return result;
    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure atanh() {
        final DerivativeStructure result = factory.build();
        factory.getCompiler().atanh(data, 0, result.data, 0);
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public DerivativeStructure toDegrees() {
        final DerivativeStructure ds = factory.build();
        for (int i = 0; i < ds.data.length; ++i) {
            ds.data[i] = FastMath.toDegrees(data[i]);
        }
        return ds;
    }

    /** {@inheritDoc} */
    @Override
    public DerivativeStructure toRadians() {
        final DerivativeStructure ds = factory.build();
        for (int i = 0; i < ds.data.length; ++i) {
            ds.data[i] = FastMath.toRadians(data[i]);
        }
        return ds;
    }

    /** Integrate w.r.t. one independent variable.
     * <p>
     * Rigorously, if the derivatives of a function are known up to
     * order N, the ones of its M-th integral w.r.t. a given variable
     * (seen as a function itself) are actually known up to order N+M.
     * However, this method still casts the output as a DerivativeStructure
     * of order N. The integration constants are systematically set to zero.
     * </p>
     * @param varIndex Index of independent variable w.r.t. which integration is done.
     * @param integrationOrder Number of times the integration operator must be applied. If non-positive, call the
     *                         differentiation operator.
     * @return DerivativeStructure on which integration operator has been applied a certain number of times.
     * @since 2.2
     */
    public DerivativeStructure integrate(final int varIndex, final int integrationOrder) {

        // Deal first with trivial case
        if (integrationOrder > getOrder()) {
            return factory.constant(0.);
        } else if (integrationOrder == 0) {
            return factory.build(data);
        }

        // Call 'inverse' (not rigorously) operation if necessary
        if (integrationOrder < 0) {
            return differentiate(varIndex, -integrationOrder);
        }

        final double[] newData = new double[data.length];
        final DSCompiler dsCompiler = factory.getCompiler();
        for (int i = 0; i < newData.length; i++) {
            if (data[i] != 0.) {
                final int[] orders = dsCompiler.getPartialDerivativeOrders(i);
                int sum = 0;
                for (int order : orders) {
                    sum += order;
                }
                if (sum + integrationOrder <= getOrder()) {
                    final int saved = orders[varIndex];
                    orders[varIndex] += integrationOrder;
                    final int index = dsCompiler.getPartialDerivativeIndex(orders);
                    orders[varIndex] = saved;
                    newData[index] = data[i];
                }
            }
        }

        return factory.build(newData);
    }

    /** Differentiate w.r.t. one independent variable.
     * <p>
     * Rigorously, if the derivatives of a function are known up to
     * order N, the ones of its M-th derivative w.r.t. a given variable
     * (seen as a function itself) are only known up to order N-M.
     * However, this method still casts the output as a DerivativeStructure
     * of order N with zeroes for the higher order terms.
     * </p>
     * @param varIndex Index of independent variable w.r.t. which differentiation is done.
     * @param differentiationOrder Number of times the differentiation operator must be applied. If non-positive, call
     *                             the integration operator instead.
     * @return DerivativeStructure on which differentiation operator has been applied a certain number of times
     * @since 2.2
     */
    public DerivativeStructure differentiate(final int varIndex, final int differentiationOrder) {

        // Deal first with trivial case
        if (differentiationOrder > getOrder()) {
            return factory.constant(0.);
        } else if (differentiationOrder == 0) {
            return factory.build(data);
        }

        // Call 'inverse' (not rigorously) operation if necessary
        if (differentiationOrder < 0) {
            return integrate(varIndex, -differentiationOrder);
        }

        final double[] newData = new double[data.length];
        final DSCompiler dsCompiler = factory.getCompiler();
        for (int i = 0; i < newData.length; i++) {
            if (data[i] != 0.) {
                final int[] orders = dsCompiler.getPartialDerivativeOrders(i);
                if (orders[varIndex] - differentiationOrder >= 0) {
                    final int saved = orders[varIndex];
                    orders[varIndex] -= differentiationOrder;
                    final int index = dsCompiler.getPartialDerivativeIndex(orders);
                    orders[varIndex] = saved;
                    newData[index] = data[i];
                }
            }
        }

        return factory.build(newData);
    }

    /** Evaluate Taylor expansion a derivative structure.
     * @param delta parameters offsets (&Delta;x, &Delta;y, ...)
     * @return value of the Taylor expansion at x + &Delta;x, y + &Delta;y, ...
     * @throws MathRuntimeException if factorials becomes too large
     */
    public double taylor(final double ... delta) throws MathRuntimeException {
        return factory.getCompiler().taylor(data, 0, delta);
    }

    /** Rebase instance with respect to low level parameter functions.
     * <p>
     * The instance is considered to be a function of {@link #getFreeParameters()
     * n free parameters} up to order {@link #getOrder() o} \(f(p_0, p_1, \ldots p_{n-1})\).
     * Its {@link #getPartialDerivative(int...) partial derivatives} are therefore
     * \(f, \frac{\partial f}{\partial p_0}, \frac{\partial f}{\partial p_1}, \ldots
     * \frac{\partial^2 f}{\partial p_0^2}, \frac{\partial^2 f}{\partial p_0 p_1},
     * \ldots \frac{\partial^o f}{\partial p_{n-1}^o}\). The free parameters
     * \(p_0, p_1, \ldots p_{n-1}\) are considered to be functions of \(m\) lower
     * level other parameters \(q_0, q_1, \ldots q_{m-1}\).
     * </p>
     * \( \begin{align}
     * p_0 &amp; = p_0(q_0, q_1, \ldots q_{m-1})\\
     * p_1 &amp; = p_1(q_0, q_1, \ldots q_{m-1})\\
     * p_{n-1} &amp; = p_{n-1}(q_0, q_1, \ldots q_{m-1})
     * \end{align}\)
     * <p>
     * This method compute the composition of the partial derivatives of \(f\)
     * and the partial derivatives of \(p_0, p_1, \ldots p_{n-1}\), i.e. the
     * {@link #getPartialDerivative(int...) partial derivatives} of the value
     * returned will be
     * \(f, \frac{\partial f}{\partial q_0}, \frac{\partial f}{\partial q_1}, \ldots
     * \frac{\partial^2 f}{\partial q_0^2}, \frac{\partial^2 f}{\partial q_0 q_1},
     * \ldots \frac{\partial^o f}{\partial q_{m-1}^o}\).
     * </p>
     * <p>
     * The number of parameters must match {@link #getFreeParameters()} and the
     * derivation orders of the instance and parameters must also match.
     * </p>
     * @param p base parameters with respect to which partial derivatives
     * were computed in the instance
     * @return derivative structure with partial derivatives computed
     * with respect to the lower level parameters used in the \(p_i\)
     * @since 2.2
     */
    public DerivativeStructure rebase(final DerivativeStructure... p) {

        MathUtils.checkDimension(getFreeParameters(), p.length);

        // handle special case of no variables at all
        if (p.length == 0) {
            return this;
        }

        final int pSize = p[0].getFactory().getCompiler().getSize();
        final double[] pData = new double[p.length * pSize];
        for (int i = 0; i < p.length; ++i) {
            MathUtils.checkDimension(getOrder(), p[i].getOrder());
            MathUtils.checkDimension(p[0].getFreeParameters(), p[i].getFreeParameters());
            System.arraycopy(p[i].data, 0, pData, i * pSize, pSize);
        }

        final DerivativeStructure result = p[0].factory.build();
        factory.getCompiler().rebase(data, 0, p[0].factory.getCompiler(), pData, result.data, 0);
        return result;

    }

    /** {@inheritDoc}
     * @exception MathIllegalArgumentException if number of free parameters
     * or orders do not match
     */
    @Override
    public DerivativeStructure linearCombination(final DerivativeStructure[] a, final DerivativeStructure[] b)
        throws MathIllegalArgumentException {

        // compute an accurate value, taking care of cancellations
        final double[] aDouble = new double[a.length];
        for (int i = 0; i < a.length; ++i) {
            aDouble[i] = a[i].getValue();
        }
        final double[] bDouble = new double[b.length];
        for (int i = 0; i < b.length; ++i) {
            bDouble[i] = b[i].getValue();
        }
        final double accurateValue = MathArrays.linearCombination(aDouble, bDouble);

        // compute a simple value, with all partial derivatives
        DerivativeStructure simpleValue = a[0].getField().getZero();
        for (int i = 0; i < a.length; ++i) {
            simpleValue = simpleValue.add(a[i].multiply(b[i]));
        }

        // create a result with accurate value and all derivatives (not necessarily as accurate as the value)
        final double[] all = simpleValue.getAllDerivatives();
        all[0] = accurateValue;
        return factory.build(all);

    }

    /** {@inheritDoc}
     * @exception MathIllegalArgumentException if number of free parameters
     * or orders do not match
     */
    @Override
    public DerivativeStructure linearCombination(final double[] a, final DerivativeStructure[] b)
        throws MathIllegalArgumentException {

        // compute an accurate value, taking care of cancellations
        final double[] bDouble = new double[b.length];
        for (int i = 0; i < b.length; ++i) {
            bDouble[i] = b[i].getValue();
        }
        final double accurateValue = MathArrays.linearCombination(a, bDouble);

        // compute a simple value, with all partial derivatives
        DerivativeStructure simpleValue = b[0].getField().getZero();
        for (int i = 0; i < a.length; ++i) {
            simpleValue = simpleValue.add(b[i].multiply(a[i]));
        }

        // create a result with accurate value and all derivatives (not necessarily as accurate as the value)
        final double[] all = simpleValue.getAllDerivatives();
        all[0] = accurateValue;
        return factory.build(all);

    }

    /** {@inheritDoc}
     * @exception MathIllegalArgumentException if number of free parameters
     * or orders do not match
     */
    @Override
    public DerivativeStructure linearCombination(final DerivativeStructure a1, final DerivativeStructure b1,
                                                 final DerivativeStructure a2, final DerivativeStructure b2)
        throws MathIllegalArgumentException {

        // compute an accurate value, taking care of cancellations
        final double accurateValue = MathArrays.linearCombination(a1.getValue(), b1.getValue(),
                                                                  a2.getValue(), b2.getValue());

        // compute a simple value, with all partial derivatives
        final DerivativeStructure simpleValue = a1.multiply(b1).add(a2.multiply(b2));

        // create a result with accurate value and all derivatives (not necessarily as accurate as the value)
        final double[] all = simpleValue.getAllDerivatives();
        all[0] = accurateValue;
        return factory.build(all);

    }

    /** {@inheritDoc}
     * @exception MathIllegalArgumentException if number of free parameters
     * or orders do not match
     */
    @Override
    public DerivativeStructure linearCombination(final double a1, final DerivativeStructure b1,
                                                 final double a2, final DerivativeStructure b2)
        throws MathIllegalArgumentException {

        factory.checkCompatibility(b1.factory);
        factory.checkCompatibility(b2.factory);

        final DerivativeStructure ds = factory.build();
        factory.getCompiler().linearCombination(a1, b1.data, 0,
                                                a2, b2.data, 0,
                                                ds.data, 0);

        return ds;

    }

    /** {@inheritDoc}
     * @exception MathIllegalArgumentException if number of free parameters
     * or orders do not match
     */
    @Override
    public DerivativeStructure linearCombination(final DerivativeStructure a1, final DerivativeStructure b1,
                                                 final DerivativeStructure a2, final DerivativeStructure b2,
                                                 final DerivativeStructure a3, final DerivativeStructure b3)
        throws MathIllegalArgumentException {

        // compute an accurate value, taking care of cancellations
        final double accurateValue = MathArrays.linearCombination(a1.getValue(), b1.getValue(),
                                                                  a2.getValue(), b2.getValue(),
                                                                  a3.getValue(), b3.getValue());

        // compute a simple value, with all partial derivatives
        final DerivativeStructure simpleValue = a1.multiply(b1).add(a2.multiply(b2)).add(a3.multiply(b3));

        // create a result with accurate value and all derivatives (not necessarily as accurate as the value)
        final double[] all = simpleValue.getAllDerivatives();
        all[0] = accurateValue;
        return factory.build(all);

    }

    /** {@inheritDoc}
     * @exception MathIllegalArgumentException if number of free parameters
     * or orders do not match
     */
    @Override
    public DerivativeStructure linearCombination(final double a1, final DerivativeStructure b1,
                                                 final double a2, final DerivativeStructure b2,
                                                 final double a3, final DerivativeStructure b3)
        throws MathIllegalArgumentException {

        factory.checkCompatibility(b1.factory);
        factory.checkCompatibility(b2.factory);
        factory.checkCompatibility(b3.factory);

        final DerivativeStructure ds = factory.build();
        factory.getCompiler().linearCombination(a1, b1.data, 0,
                                                a2, b2.data, 0,
                                                a3, b3.data, 0,
                                                ds.data, 0);

        return ds;

    }

    /** {@inheritDoc}
     * @exception MathIllegalArgumentException if number of free parameters
     * or orders do not match
     */
    @Override
    public DerivativeStructure linearCombination(final DerivativeStructure a1, final DerivativeStructure b1,
                                                 final DerivativeStructure a2, final DerivativeStructure b2,
                                                 final DerivativeStructure a3, final DerivativeStructure b3,
                                                 final DerivativeStructure a4, final DerivativeStructure b4)
        throws MathIllegalArgumentException {

        // compute an accurate value, taking care of cancellations
        final double accurateValue = MathArrays.linearCombination(a1.getValue(), b1.getValue(),
                                                                  a2.getValue(), b2.getValue(),
                                                                  a3.getValue(), b3.getValue(),
                                                                  a4.getValue(), b4.getValue());

        // compute a simple value, with all partial derivatives
        final DerivativeStructure simpleValue = a1.multiply(b1).add(a2.multiply(b2)).add(a3.multiply(b3)).add(a4.multiply(b4));

        // create a result with accurate value and all derivatives (not necessarily as accurate as the value)
        final double[] all = simpleValue.getAllDerivatives();
        all[0] = accurateValue;
        return factory.build(all);

    }

    /** {@inheritDoc}
     * @exception MathIllegalArgumentException if number of free parameters
     * or orders do not match
     */
    @Override
    public DerivativeStructure linearCombination(final double a1, final DerivativeStructure b1,
                                                 final double a2, final DerivativeStructure b2,
                                                 final double a3, final DerivativeStructure b3,
                                                 final double a4, final DerivativeStructure b4)
        throws MathIllegalArgumentException {

        factory.checkCompatibility(b1.factory);
        factory.checkCompatibility(b2.factory);
        factory.checkCompatibility(b3.factory);
        factory.checkCompatibility(b4.factory);

        final DerivativeStructure ds = factory.build();
        factory.getCompiler().linearCombination(a1, b1.data, 0,
                                                a2, b2.data, 0,
                                                a3, b3.data, 0,
                                                a4, b4.data, 0,
                                                ds.data, 0);

        return ds;

    }

    /** {@inheritDoc}
     */
    @Override
    public DerivativeStructure getPi() {
        return factory.getDerivativeField().getPi();
    }

    /**
     * Test for the equality of two derivative structures.
     * <p>
     * Derivative structures are considered equal if they have the same number
     * of free parameters, the same derivation order, and the same derivatives.
     * </p>
     * @param other Object to test for equality to this
     * @return true if two derivative structures are equal
     */
    @Override
    public boolean equals(Object other) {

        if (this == other) {
            return true;
        }

        if (other instanceof DerivativeStructure) {
            final DerivativeStructure rhs = (DerivativeStructure)other;
            return (getFreeParameters() == rhs.getFreeParameters()) &&
                   (getOrder() == rhs.getOrder()) &&
                   MathArrays.equals(data, rhs.data);
        }

        return false;

    }

    /**
     * Get a hashCode for the derivative structure.
     * @return a hash code value for this object
     */
    @Override
    public int hashCode() {
        return 227 + 229 * getFreeParameters() + 233 * getOrder() + 239 * MathUtils.hash(data);
    }

}
