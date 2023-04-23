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

package org.hipparchus.ode.nonstiff;

import org.hipparchus.Field;
import org.hipparchus.CalculusFieldElement;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.linear.Array2DRowFieldMatrix;
import org.hipparchus.linear.FieldMatrix;
import org.hipparchus.ode.FieldEquationsMapper;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.util.FastMath;


/**
 * This class implements explicit Adams-Bashforth integrators for Ordinary
 * Differential Equations.
 *
 * <p>Adams-Bashforth methods (in fact due to Adams alone) are explicit
 * multistep ODE solvers. This implementation is a variation of the classical
 * one: it uses adaptive stepsize to implement error control, whereas
 * classical implementations are fixed step size. The value of state vector
 * at step n+1 is a simple combination of the value at step n and of the
 * derivatives at steps n, n-1, n-2 ... Depending on the number k of previous
 * steps one wants to use for computing the next value, different formulas
 * are available:</p>
 * <ul>
 *   <li>k = 1: y<sub>n+1</sub> = y<sub>n</sub> + h y'<sub>n</sub></li>
 *   <li>k = 2: y<sub>n+1</sub> = y<sub>n</sub> + h (3y'<sub>n</sub>-y'<sub>n-1</sub>)/2</li>
 *   <li>k = 3: y<sub>n+1</sub> = y<sub>n</sub> + h (23y'<sub>n</sub>-16y'<sub>n-1</sub>+5y'<sub>n-2</sub>)/12</li>
 *   <li>k = 4: y<sub>n+1</sub> = y<sub>n</sub> + h (55y'<sub>n</sub>-59y'<sub>n-1</sub>+37y'<sub>n-2</sub>-9y'<sub>n-3</sub>)/24</li>
 *   <li>...</li>
 * </ul>
 *
 * <p>A k-steps Adams-Bashforth method is of order k.</p>
 *
 * <p> There must be sufficient time for the {@link #setStarterIntegrator(org.hipparchus.ode.FieldODEIntegrator)
 * starter integrator} to take several steps between the the last reset event, and the end
 * of integration, otherwise an exception may be thrown during integration. The user can
 * adjust the end date of integration, or the step size of the starter integrator to
 * ensure a sufficient number of steps can be completed before the end of integration.
 * </p>
 *
 * <p><strong>Implementation details</strong></p>
 *
 * <p>We define scaled derivatives s<sub>i</sub>(n) at step n as:
 * \[
 *   \left\{\begin{align}
 *   s_1(n) &amp;= h y'_n \text{ for first derivative}\\
 *   s_2(n) &amp;= \frac{h^2}{2} y_n'' \text{ for second derivative}\\
 *   s_3(n) &amp;= \frac{h^3}{6} y_n''' \text{ for third derivative}\\
 *   &amp;\cdots\\
 *   s_k(n) &amp;= \frac{h^k}{k!} y_n^{(k)} \text{ for } k^\mathrm{th} \text{ derivative}
 *   \end{align}\right.
 * \]</p>
 *
 * <p>The definitions above use the classical representation with several previous first
 * derivatives. Lets define
 * \[
 *   q_n = [ s_1(n-1) s_1(n-2) \ldots s_1(n-(k-1)) ]^T
 * \]
 * (we omit the k index in the notation for clarity). With these definitions,
 * Adams-Bashforth methods can be written:
 * \[
 *   \left\{\begin{align}
 *   k = 1: &amp; y_{n+1} = y_n +               s_1(n) \\
 *   k = 2: &amp; y_{n+1} = y_n + \frac{3}{2}   s_1(n) + [ \frac{-1}{2} ] q_n \\
 *   k = 3: &amp; y_{n+1} = y_n + \frac{23}{12} s_1(n) + [ \frac{-16}{12} \frac{5}{12} ] q_n \\
 *   k = 4: &amp; y_{n+1} = y_n + \frac{55}{24} s_1(n) + [ \frac{-59}{24} \frac{37}{24} \frac{-9}{24} ] q_n \\
 *          &amp; \cdots
 *   \end{align}\right.
 * \]
 *
 * <p>Instead of using the classical representation with first derivatives only (y<sub>n</sub>,
 * s<sub>1</sub>(n) and q<sub>n</sub>), our implementation uses the Nordsieck vector with
 * higher degrees scaled derivatives all taken at the same step (y<sub>n</sub>, s<sub>1</sub>(n)
 * and r<sub>n</sub>) where r<sub>n</sub> is defined as:
 * \[
 * r_n = [ s_2(n), s_3(n) \ldots s_k(n) ]^T
 * \]
 * (here again we omit the k index in the notation for clarity)
 * </p>
 *
 * <p>Taylor series formulas show that for any index offset i, s<sub>1</sub>(n-i) can be
 * computed from s<sub>1</sub>(n), s<sub>2</sub>(n) ... s<sub>k</sub>(n), the formula being exact
 * for degree k polynomials.
 * \[
 * s_1(n-i) = s_1(n) + \sum_{j\gt 0} (j+1) (-i)^j s_{j+1}(n)
 * \]
 * The previous formula can be used with several values for i to compute the transform between
 * classical representation and Nordsieck vector. The transform between r<sub>n</sub>
 * and q<sub>n</sub> resulting from the Taylor series formulas above is:
 * \[
 * q_n = s_1(n) u + P r_n
 * \]
 * where u is the [ 1 1 ... 1 ]<sup>T</sup> vector and P is the (k-1)&times;(k-1) matrix built
 * with the \((j+1) (-i)^j\) terms with i being the row number starting from 1 and j being
 * the column number starting from 1:
 * \[
 *   P=\begin{bmatrix}
 *   -2  &amp;  3 &amp;   -4 &amp;    5 &amp; \ldots \\
 *   -4  &amp; 12 &amp;  -32 &amp;   80 &amp; \ldots \\
 *   -6  &amp; 27 &amp; -108 &amp;  405 &amp; \ldots \\
 *   -8  &amp; 48 &amp; -256 &amp; 1280 &amp; \ldots \\
 *       &amp;    &amp;  \ldots\\
 *    \end{bmatrix}
 * \]
 *
 * <p>Using the Nordsieck vector has several advantages:</p>
 * <ul>
 *   <li>it greatly simplifies step interpolation as the interpolator mainly applies
 *   Taylor series formulas,</li>
 *   <li>it simplifies step changes that occur when discrete events that truncate
 *   the step are triggered,</li>
 *   <li>it allows to extend the methods in order to support adaptive stepsize.</li>
 * </ul>
 *
 * <p>The Nordsieck vector at step n+1 is computed from the Nordsieck vector at step n as follows:
 * <ul>
 *   <li>y<sub>n+1</sub> = y<sub>n</sub> + s<sub>1</sub>(n) + u<sup>T</sup> r<sub>n</sub></li>
 *   <li>s<sub>1</sub>(n+1) = h f(t<sub>n+1</sub>, y<sub>n+1</sub>)</li>
 *   <li>r<sub>n+1</sub> = (s<sub>1</sub>(n) - s<sub>1</sub>(n+1)) P<sup>-1</sup> u + P<sup>-1</sup> A P r<sub>n</sub></li>
 * </ul>
 * <p>where A is a rows shifting matrix (the lower left part is an identity matrix):</p>
 * <pre>
 *        [ 0 0   ...  0 0 | 0 ]
 *        [ ---------------+---]
 *        [ 1 0   ...  0 0 | 0 ]
 *    A = [ 0 1   ...  0 0 | 0 ]
 *        [       ...      | 0 ]
 *        [ 0 0   ...  1 0 | 0 ]
 *        [ 0 0   ...  0 1 | 0 ]
 * </pre>
 *
 * <p>The P<sup>-1</sup>u vector and the P<sup>-1</sup> A P matrix do not depend on the state,
 * they only depend on k and therefore are precomputed once for all.</p>
 *
 * @param <T> the type of the field elements
 */
public class AdamsBashforthFieldIntegrator<T extends CalculusFieldElement<T>> extends AdamsFieldIntegrator<T> {

    /** Integrator method name. */
    private static final String METHOD_NAME = "Adams-Bashforth";

    /**
     * Build an Adams-Bashforth integrator with the given order and step control parameters.
     * @param field field to which the time and state vector elements belong
     * @param nSteps number of steps of the method excluding the one being computed
     * @param minStep minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param maxStep maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param scalAbsoluteTolerance allowed absolute error
     * @param scalRelativeTolerance allowed relative error
     * @exception MathIllegalArgumentException if order is 1 or less
     */
    public AdamsBashforthFieldIntegrator(final Field<T> field, final int nSteps,
                                         final double minStep, final double maxStep,
                                         final double scalAbsoluteTolerance,
                                         final double scalRelativeTolerance)
        throws MathIllegalArgumentException {
        super(field, METHOD_NAME, nSteps, nSteps, minStep, maxStep,
              scalAbsoluteTolerance, scalRelativeTolerance);
    }

    /**
     * Build an Adams-Bashforth integrator with the given order and step control parameters.
     * @param field field to which the time and state vector elements belong
     * @param nSteps number of steps of the method excluding the one being computed
     * @param minStep minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param maxStep maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param vecAbsoluteTolerance allowed absolute error
     * @param vecRelativeTolerance allowed relative error
     * @exception IllegalArgumentException if order is 1 or less
     */
    public AdamsBashforthFieldIntegrator(final Field<T> field, final int nSteps,
                                         final double minStep, final double maxStep,
                                         final double[] vecAbsoluteTolerance,
                                         final double[] vecRelativeTolerance)
        throws IllegalArgumentException {
        super(field, METHOD_NAME, nSteps, nSteps, minStep, maxStep,
              vecAbsoluteTolerance, vecRelativeTolerance);
    }

    /** {@inheritDoc} */
    @Override
    protected double errorEstimation(final T[] previousState, final T predictedTime,
                                     final T[] predictedState, final T[] predictedScaled,
                                     final FieldMatrix<T> predictedNordsieck) {

        final StepsizeHelper helper = getStepSizeHelper();
        double error = 0;
        for (int i = 0; i < helper.getMainSetDimension(); ++i) {
            final double tol = helper.getTolerance(i, FastMath.abs(predictedState[i].getReal()));

            // apply Taylor formula from high order to low order,
            // for the sake of numerical accuracy
            double variation = 0;
            int sign = predictedNordsieck.getRowDimension() % 2 == 0 ? -1 : 1;
            for (int k = predictedNordsieck.getRowDimension() - 1; k >= 0; --k) {
                variation += sign * predictedNordsieck.getEntry(k, i).getReal();
                sign       = -sign;
            }
            variation -= predictedScaled[i].getReal();

            final double ratio  = (predictedState[i].getReal() - previousState[i].getReal() + variation) / tol;
            error              += ratio * ratio;

        }

        return FastMath.sqrt(error / helper.getMainSetDimension());

    }

    /** {@inheritDoc} */
    @Override
    protected AdamsFieldStateInterpolator<T> finalizeStep(final T stepSize, final T[] predictedY,
                                                          final T[] predictedScaled, final Array2DRowFieldMatrix<T> predictedNordsieck,
                                                          final boolean isForward,
                                                          final FieldODEStateAndDerivative<T> globalPreviousState,
                                                          final FieldODEStateAndDerivative<T> globalCurrentState,
                                                          final FieldEquationsMapper<T> equationsMapper) {
        return new AdamsFieldStateInterpolator<>(getStepSize(), globalCurrentState,
                                                 predictedScaled, predictedNordsieck, isForward,
                                                 getStepStart(), globalCurrentState, equationsMapper);
    }

}
