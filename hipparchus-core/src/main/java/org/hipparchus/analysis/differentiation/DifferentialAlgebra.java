package org.hipparchus.analysis.differentiation;

/** Interface representing an object holding partial derivatives.
 * @see Derivative
 * @see TaylorMap
 * @see FieldDerivative
 * @see FieldTaylorMap
 * @since 3.1
 */
public interface DifferentialAlgebra {

    /** Get the number of free parameters.
     * @return number of free parameters
     */
    int getFreeParameters();

    /** Get the maximum derivation order.
     * @return maximum derivation order
     */
    int getOrder();

}
