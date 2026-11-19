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
import org.hipparchus.util.FastMath;
import org.hipparchus.util.FieldSinCos;
import org.hipparchus.util.FieldSinhCosh;

/** Interface representing a Field object holding partial derivatives up to first order.
 * @param <S> the type of the field elements
 * @param <T> the type of the function derivative
 * @see FieldDerivative
 * @see FieldUnivariateDerivative1
 * @see FieldGradient
 * @see Derivative1
 * @since 3.1
 */
public interface FieldDerivative1<S extends CalculusFieldElement<S>, T extends FieldDerivative<S, T>>
        extends FieldDerivative<S, T> {

    /** {@inheritDoc} */
    @Override
    default int getOrder() {
        return 1;
    }

    /** Compute composition of the instance by a univariate function differentiable at order 1.
     * @param f0 value of function
     * @param f1 first-order derivative
     * @return f(this)
     */
    T compose(S f0, S f1);

    /** {@inheritDoc} */
    @Override
    default T square() {
        final S f0 = getValue();
        return compose(f0.square(), f0.twice());
    }

    /** {@inheritDoc} */
    @Override
    default T reciprocal() {
        final S inv1 = getValue().reciprocal();
        final S inv2 = inv1.square().negate();
        return compose(inv1, inv2);
    }

    /** {@inheritDoc} */
    @Override
    default T exp() {
        final S exp = getValue().exp();
        return compose(exp, exp);
    }

    /** {@inheritDoc} */
    @Override
    default T sqrt() {
        final S s = getValue().sqrt();
        return compose(s, s.add(s).reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    default T cbrt() {
        final S c = getValue().cbrt();
        return compose(c, c.square().multiply(3).reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    default T expm1() {
        final S exp   = FastMath.exp(getValue());
        final S expM1 = FastMath.expm1(getValue());
        return compose(expM1, exp);
    }

    /** {@inheritDoc} */
    @Override
    default T log() {
        return compose(getValue().log(), getValue().reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    default T log1p() {
        return compose(getValue().log1p(), getValue().add(1).reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    default T log10() {
        return compose(getValue().log10(), getValue().multiply(FastMath.log(10.0)).reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    default T cos() {
        final FieldSinCos<S> sinCos = getValue().sinCos();
        return compose(sinCos.cos(), sinCos.sin().negate());
    }

    /** {@inheritDoc} */
    @Override
    default T sin() {
        final FieldSinCos<S> sinCos = getValue().sinCos();
        return compose(sinCos.sin(), sinCos.cos());
    }

    /** {@inheritDoc} */
    @Override
    default FieldSinCos<T> sinCos() {
        final FieldSinCos<S> sinCos = getValue().sinCos();
        return new FieldSinCos<>(compose(sinCos.sin(), sinCos.cos()),
                compose(sinCos.cos(), sinCos.sin().negate()));
    }

    /** {@inheritDoc} */
    @Override
    default T tan() {
        final S tan = getValue().tan();
        return compose(tan, tan.multiply(tan).add(1));
    }

    /** {@inheritDoc} */
    @Override
    default T acos() {
        return compose(getValue().acos(), getValue().square().negate().add(1).sqrt().reciprocal().negate());
    }

    /** {@inheritDoc} */
    @Override
    default T asin() {
        return compose(getValue().asin(), getValue().square().negate().add(1).sqrt().reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    default T atan() {
        return compose(getValue().atan(), getValue().square().add(1).reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    default T cosh() {
        final FieldSinhCosh<S> sinhCosh = getValue().sinhCosh();
        return compose(sinhCosh.cosh(), sinhCosh.sinh());
    }

    /** {@inheritDoc} */
    @Override
    default T sinh() {
        final FieldSinhCosh<S> sinhCosh = getValue().sinhCosh();
        return compose(sinhCosh.sinh(), sinhCosh.cosh());
    }

    /** {@inheritDoc} */
    @Override
    default FieldSinhCosh<T> sinhCosh() {
        final FieldSinhCosh<S> sinhCosh = getValue().sinhCosh();
        return new FieldSinhCosh<>(compose(sinhCosh.sinh(), sinhCosh.cosh()),
                compose(sinhCosh.cosh(), sinhCosh.sinh()));
    }

    /** {@inheritDoc} */
    @Override
    default T tanh() {
        final S tanh = getValue().tanh();
        return compose(tanh, tanh.multiply(tanh).negate().add(1));
    }

    /** {@inheritDoc} */
    @Override
    default T acosh() {
        return compose(getValue().acosh(), getValue().square().subtract(1).sqrt().reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    default T asinh() {
        return compose(getValue().asinh(), getValue().square().add(1).sqrt().reciprocal());
    }

    /** {@inheritDoc} */
    @Override
    default T atanh() {
        return compose(getValue().atanh(), getValue().square().negate().add(1).reciprocal());
    }

}
