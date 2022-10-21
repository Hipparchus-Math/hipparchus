package org.hipparchus.analysis.polynomials;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.util.CombinatoricsUtils;
import org.hipparchus.util.FastMath;

/**
 * Smoothstep function factory.
 * <p>
 * It allows for quick creation of common/generic smoothstep functions.
 *
 * @author Vincent Cucchietti
 */
public class SmoothStepFactory {

    /**
     * Get the {@link SmoothStepFunction clamping smoothstep function}.
     * @return clamping smoothstep function
     */
    public static SmoothStepFunction getClamp() {
        return getGeneralOrder(0);
    }

    /**
     * Get the {@link SmoothStepFunction cubic smoothstep function}.
     * @return cubic smoothstep function
     */
    public static SmoothStepFunction getCubic() {
        return getGeneralOrder(1);
    }

    /**
     * Get the {@link SmoothStepFunction quintic smoothstep function}.
     * @return quintic smoothstep function
     */
    public static SmoothStepFunction getQuintic() {
        return getGeneralOrder(2);
    }

    /**
     * Create a {@link SmoothStepFunction smoothstep function} of order <b>2N + 1</b>.
     * <p>
     * It uses the general smoothstep equation presented <a href="https://en.wikipedia.org/wiki/Smoothstep">here</a>.
     *
     * @param N determines the order of the output smoothstep function (=2N + 1)
     * @return smoothstep function of order <b>2N + 1</b>
     */
    public static SmoothStepFunction getGeneralOrder(final int N) {

        final int twoNPlusOne = 2 * N + 1;

        final double[] coefficients = new double[twoNPlusOne + 1];

        int n = N;
        for (int i = twoNPlusOne; i > N; i--) {

            coefficients[i] = FastMath.pow(-1, n)
                    * CombinatoricsUtils.binomialCoefficient(N + n, n)
                    * CombinatoricsUtils.binomialCoefficient(twoNPlusOne, N - n);

            n--;
        }

        return new SmoothStepFunction(coefficients);
    }

    /**
     * Smoothstep function, it is used to do a smooth transition between the "left edge" and the "right edge" with left
     * edge assumed to be smaller than the right edge.
     * <p>
     * By definition, for order n>1 and input x, a smoothstep function respects the following properties :
     * <ul>
     *     <li>f(x <= leftEdge) = leftEdge and f(x >= rightEdge) = rightEdge</li>
     *     <li>f'(leftEdge) = f'(rightEdge) = 0</li>
     * </ul>
     * If x is normalized between edges, we have :
     * <ul>
     *     <li>f(x <= 0) = 0 and f(x >= 1) = 1</li>
     *     <li>f'(0) = f'(1) = 0</li>
     * </ul>
     */
    public static class SmoothStepFunction extends PolynomialFunction {

        /**
         * Construct a smoothstep with the given coefficients.  The first element of the coefficients array is the
         * constant term.  Higher degree coefficients follow in sequence.  The degree of the resulting polynomial is the
         * index of the last non-null element of the array, or 0 if all elements are null.
         * <p>
         * The constructor makes a copy of the input array and assigns the copy to the coefficients property.</p>
         *
         * @param c Smoothstep polynomial coefficients.
         * @throws NullArgumentException if {@code c} is {@code null}.
         * @throws MathIllegalArgumentException if {@code c} is empty.
         */
        private SmoothStepFunction(final double[] c) throws MathIllegalArgumentException, NullArgumentException {
            super(c);
        }

        /**
         * Compute the value of the smoothstep for the given argument normalized between edges.
         *
         * @param xNormalized Normalized argument for which the function value should be computed.
         * @return the value of the polynomial at the given point.
         * @see org.hipparchus.analysis.UnivariateFunction#value(double)
         */
        @Override
        public double value(final double xNormalized) {

            if (xNormalized <= 0) {
                return 0;
            }
            if (xNormalized >= 1) {
                return 1;
            }

            return super.value(xNormalized);
        }

        /**
         * Compute the value of the smoothstep function for the given edges and argument.
         * <p>
         * This implementation is based on "Advanced Real-Time Shader Techniques". AND. p. 94. Retrieved 2022-04-16. by
         * Natalya Tatarchuk (2003).
         * </p>
         *
         * @param x Argument for which the function value should be computed.
         * @return the value of the polynomial at the given point.
         * @see org.hipparchus.analysis.UnivariateFunction#value(double)
         */
        public double value(final double leftEdge, final double rightEdge, final double x)
                throws MathIllegalArgumentException {

            if (leftEdge > rightEdge) {
                throw new MathIllegalArgumentException(LocalizedCoreFormats.RIGHT_EDGE_GREATER_THAN_LEFT_EDGE,
                                                       leftEdge, rightEdge);
            }
            if (x <= leftEdge) {
                return 0;
            }
            if (x >= rightEdge) {
                return 1;
            }

            final double xNormalized = (x - leftEdge) / (rightEdge - leftEdge);

            return super.value(xNormalized);
        }
    }

}
