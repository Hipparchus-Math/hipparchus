package org.hipparchus.analysis.differentiation;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.SinCos;
import org.hipparchus.util.SinhCosh;
import org.hipparchus.util.FieldSinCos;
import org.hipparchus.util.FieldSinhCosh;

/** Interface representing an object holding partial derivatives up to first order.
 * @param <T> the type of the field elements
 * @see Derivative
 * @see UnivariateDerivative1
 * @see Gradient
 * @see SparseGradient
 * @since 3.1
 */
public interface Derivative1<T extends CalculusFieldElement<T>> extends Derivative<T> {

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
    T compose(double f0, double f1);

    /** {@inheritDoc} */
    @Override
    default T square() {
        final double f0 = getValue();
        return compose(f0 * f0, 2 * f0);
    }

    /** {@inheritDoc} */
    @Override
    default T reciprocal () {
        final double inv1 = 1.0 / getValue();
        final double inv2 = -inv1 * inv1;
        return compose(inv1, inv2);
    }

    /** {@inheritDoc} */
    @Override
    default T sqrt() {
        final double s = FastMath.sqrt(getValue());
        return compose(s, 1 / (2 * s));
    }

    /** {@inheritDoc} */
    @Override
    default T cbrt() {
        final double c = FastMath.cbrt(getValue());
        return compose(c, 1 / (3 * c * c));
    }

    /** {@inheritDoc} */
    @Override
    default T rootN(int n) {
        if (n == 2) {
            return sqrt();
        } else if (n == 3) {
            return cbrt();
        } else {
            final double r = FastMath.pow(getValue(), 1.0 / n);
            return compose(r, 1 / (n * FastMath.pow(r, n - 1)));
        }
    }

    /** {@inheritDoc} */
    @Override
    default T exp() {
        final double exp = FastMath.exp(getValue());
        return compose(exp, exp);
    }

    /** {@inheritDoc} */
    @Override
    default T expm1() {
        final double exp   = FastMath.exp(getValue());
        final double expM1 = FastMath.expm1(getValue());
        return compose(expM1, exp);
    }

    /** {@inheritDoc} */
    @Override
    default T log() {
        return compose(FastMath.log(getValue()), 1 / getValue());
    }

    /** {@inheritDoc} */
    @Override
    default T log1p() {
        return compose(FastMath.log1p(getValue()), 1 / (1 + getValue()));
    }

    /** {@inheritDoc} */
    @Override
    default T log10() {
        return compose(FastMath.log10(getValue()), 1 / (getValue() * FastMath.log(10.0)));
    }

    /** {@inheritDoc} */
    @Override
    default T cos() {
        final SinCos sinCos = FastMath.sinCos(getValue());
        return compose(sinCos.cos(), -sinCos.sin());
    }

    /** {@inheritDoc} */
    @Override
    default T sin() {
        final SinCos sinCos = FastMath.sinCos(getValue());
        return compose(sinCos.sin(), sinCos.cos());
    }

    /** {@inheritDoc} */
    @Override
    default FieldSinCos<T> sinCos() {
        final SinCos sinCos = FastMath.sinCos(getValue());
        return new FieldSinCos<>(compose(sinCos.sin(), sinCos.cos()),
                compose(sinCos.cos(), -sinCos.sin()));
    }

    /** {@inheritDoc} */
    @Override
    default T tan() {
        final double tan = FastMath.tan(getValue());
        return compose(tan, 1 + tan * tan);
    }

    /** {@inheritDoc} */
    @Override
    default T acos() {
        final double f0 = getValue();
        return compose(FastMath.acos(f0), -1 / FastMath.sqrt(1 - f0 * f0));
    }

    /** {@inheritDoc} */
    @Override
    default T asin() {
        final double f0 = getValue();
        return compose(FastMath.asin(f0), 1 / FastMath.sqrt(1 - f0 * f0));
    }

    /** {@inheritDoc} */
    @Override
    default T atan() {
        final double f0 = getValue();
        return compose(FastMath.atan(f0), 1 / (1 + f0 * f0));
    }

    /** {@inheritDoc} */
    @Override
    default T cosh() {
        return compose(FastMath.cosh(getValue()), FastMath.sinh(getValue()));
    }

    /** {@inheritDoc} */
    @Override
    default T sinh() {
        return compose(FastMath.sinh(getValue()), FastMath.cosh(getValue()));
    }

    /** {@inheritDoc} */
    @Override
    default FieldSinhCosh<T> sinhCosh() {
        final SinhCosh sinhCosh = FastMath.sinhCosh(getValue());
        return new FieldSinhCosh<>(compose(sinhCosh.sinh(), sinhCosh.cosh()),
                compose(sinhCosh.cosh(), sinhCosh.sinh()));
    }

    /** {@inheritDoc} */
    @Override
    default T tanh() {
        final double tanh = FastMath.tanh(getValue());
        return compose(tanh, 1 - tanh * tanh);
    }

    /** {@inheritDoc} */
    @Override
    default T acosh() {
        final double f0 = getValue();
        return compose(FastMath.acosh(f0), 1 / FastMath.sqrt(f0 * f0 - 1));
    }

    /** {@inheritDoc} */
    @Override
    default T asinh() {
        final double f0 = getValue();
        return compose(FastMath.asinh(f0), 1 / FastMath.sqrt(f0 * f0 + 1));
    }

    /** {@inheritDoc} */
    @Override
    default T atanh() {
        final double f0 = getValue();
        return compose(FastMath.atanh(f0), 1 / (1 - f0 * f0));
    }

}
