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

import java.util.Arrays;

import org.hipparchus.Field;
import org.hipparchus.CalculusFieldElement;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.linear.Array2DRowFieldMatrix;
import org.hipparchus.linear.FieldMatrix;
import org.hipparchus.linear.FieldMatrixPreservingVisitor;
import org.hipparchus.ode.FieldEquationsMapper;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.ode.LocalizedODEFormats;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.MathUtils;


/**
 * This class implements implicit Adams-Moulton integrators for Ordinary
 * Differential Equations.
 *
 * <p>Adams-Moulton methods (in fact due to Adams alone) are implicit
 * multistep ODE solvers. This implementation is a variation of the classical
 * one: it uses adaptive stepsize to implement error control, whereas
 * classical implementations are fixed step size. The value of state vector
 * at step n+1 is a simple combination of the value at step n and of the
 * derivatives at steps n+1, n, n-1 ... Since y'<sub>n+1</sub> is needed to
 * compute y<sub>n+1</sub>, another method must be used to compute a first
 * estimate of y<sub>n+1</sub>, then compute y'<sub>n+1</sub>, then compute
 * a final estimate of y<sub>n+1</sub> using the following formulas. Depending
 * on the number k of previous steps one wants to use for computing the next
 * value, different formulas are available for the final estimate:</p>
 * <ul>
 *   <li>k = 1: y<sub>n+1</sub> = y<sub>n</sub> + h y'<sub>n+1</sub></li>
 *   <li>k = 2: y<sub>n+1</sub> = y<sub>n</sub> + h (y'<sub>n+1</sub>+y'<sub>n</sub>)/2</li>
 *   <li>k = 3: y<sub>n+1</sub> = y<sub>n</sub> + h (5y'<sub>n+1</sub>+8y'<sub>n</sub>-y'<sub>n-1</sub>)/12</li>
 *   <li>k = 4: y<sub>n+1</sub> = y<sub>n</sub> + h (9y'<sub>n+1</sub>+19y'<sub>n</sub>-5y'<sub>n-1</sub>+y'<sub>n-2</sub>)/24</li>
 *   <li>...</li>
 * </ul>
 *
 * <p>A k-steps Adams-Moulton method is of order k+1.</p>
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
 * Adams-Moulton methods can be written:</p>
 * <ul>
 *   <li>k = 1: y<sub>n+1</sub> = y<sub>n</sub> + s<sub>1</sub>(n+1)</li>
 *   <li>k = 2: y<sub>n+1</sub> = y<sub>n</sub> + 1/2 s<sub>1</sub>(n+1) + [ 1/2 ] q<sub>n+1</sub></li>
 *   <li>k = 3: y<sub>n+1</sub> = y<sub>n</sub> + 5/12 s<sub>1</sub>(n+1) + [ 8/12 -1/12 ] q<sub>n+1</sub></li>
 *   <li>k = 4: y<sub>n+1</sub> = y<sub>n</sub> + 9/24 s<sub>1</sub>(n+1) + [ 19/24 -5/24 1/24 ] q<sub>n+1</sub></li>
 *   <li>...</li>
 * </ul>
 *
 * <p>Instead of using the classical representation with first derivatives only (y<sub>n</sub>,
 * s<sub>1</sub>(n+1) and q<sub>n+1</sub>), our implementation uses the Nordsieck vector with
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
 * <p>The predicted Nordsieck vector at step n+1 is computed from the Nordsieck vector at step
 * n as follows:
 * <ul>
 *   <li>Y<sub>n+1</sub> = y<sub>n</sub> + s<sub>1</sub>(n) + u<sup>T</sup> r<sub>n</sub></li>
 *   <li>S<sub>1</sub>(n+1) = h f(t<sub>n+1</sub>, Y<sub>n+1</sub>)</li>
 *   <li>R<sub>n+1</sub> = (s<sub>1</sub>(n) - S<sub>1</sub>(n+1)) P<sup>-1</sup> u + P<sup>-1</sup> A P r<sub>n</sub></li>
 * </ul>
 * where A is a rows shifting matrix (the lower left part is an identity matrix):
 * <pre>
 *        [ 0 0   ...  0 0 | 0 ]
 *        [ ---------------+---]
 *        [ 1 0   ...  0 0 | 0 ]
 *    A = [ 0 1   ...  0 0 | 0 ]
 *        [       ...      | 0 ]
 *        [ 0 0   ...  1 0 | 0 ]
 *        [ 0 0   ...  0 1 | 0 ]
 * </pre>
 * From this predicted vector, the corrected vector is computed as follows:
 * <ul>
 *   <li>y<sub>n+1</sub> = y<sub>n</sub> + S<sub>1</sub>(n+1) + [ -1 +1 -1 +1 ... &plusmn;1 ] r<sub>n+1</sub></li>
 *   <li>s<sub>1</sub>(n+1) = h f(t<sub>n+1</sub>, y<sub>n+1</sub>)</li>
 *   <li>r<sub>n+1</sub> = R<sub>n+1</sub> + (s<sub>1</sub>(n+1) - S<sub>1</sub>(n+1)) P<sup>-1</sup> u</li>
 * </ul>
 * <p>where the upper case Y<sub>n+1</sub>, S<sub>1</sub>(n+1) and R<sub>n+1</sub> represent the
 * predicted states whereas the lower case y<sub>n+1</sub>, s<sub>n+1</sub> and r<sub>n+1</sub>
 * represent the corrected states.</p>
 *
 * <p>The P<sup>-1</sup>u vector and the P<sup>-1</sup> A P matrix do not depend on the state,
 * they only depend on k and therefore are precomputed once for all.</p>
 *
 * @param <T> the type of the field elements
 */
public class AdamsMoultonFieldIntegrator<T extends CalculusFieldElement<T>> extends AdamsFieldIntegrator<T> {

    /** Integrator method name. */
    private static final String METHOD_NAME = "Adams-Moulton";

    /**
     * Build an Adams-Moulton integrator with the given order and error control parameters.
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
    public AdamsMoultonFieldIntegrator(final Field<T> field, final int nSteps,
                                       final double minStep, final double maxStep,
                                       final double scalAbsoluteTolerance,
                                       final double scalRelativeTolerance)
        throws MathIllegalArgumentException {
        super(field, METHOD_NAME, nSteps, nSteps + 1, minStep, maxStep,
              scalAbsoluteTolerance, scalRelativeTolerance);
    }

    /**
     * Build an Adams-Moulton integrator with the given order and error control parameters.
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
    public AdamsMoultonFieldIntegrator(final Field<T> field, final int nSteps,
                                       final double minStep, final double maxStep,
                                       final double[] vecAbsoluteTolerance,
                                       final double[] vecRelativeTolerance)
        throws IllegalArgumentException {
        super(field, METHOD_NAME, nSteps, nSteps + 1, minStep, maxStep,
              vecAbsoluteTolerance, vecRelativeTolerance);
    }

    /** {@inheritDoc} */
    @Override
    protected double errorEstimation(final T[] previousState, final T predictedTime,
                                     final T[] predictedState, final T[] predictedScaled,
                                     final FieldMatrix<T> predictedNordsieck) {
        final double error = predictedNordsieck.walkInOptimizedOrder(new Corrector(previousState, predictedScaled, predictedState)).getReal();
        if (Double.isNaN(error)) {
            throw new MathIllegalStateException(LocalizedODEFormats.NAN_APPEARING_DURING_INTEGRATION,
                                                predictedTime.getReal());
        }
        return error;
    }

    /** {@inheritDoc} */
    @Override
    protected AdamsFieldStateInterpolator<T> finalizeStep(final T stepSize, final T[] predictedY,
                                                          final T[] predictedScaled, final Array2DRowFieldMatrix<T> predictedNordsieck,
                                                          final boolean isForward,
                                                          final FieldODEStateAndDerivative<T> globalPreviousState,
                                                          final FieldODEStateAndDerivative<T> globalCurrentState,
                                                          final FieldEquationsMapper<T> equationsMapper) {

        final T[] correctedYDot = computeDerivatives(globalCurrentState.getTime(), predictedY);

        // update Nordsieck vector
        final T[] correctedScaled = MathArrays.buildArray(getField(), predictedY.length);
        for (int j = 0; j < correctedScaled.length; ++j) {
            correctedScaled[j] = getStepSize().multiply(correctedYDot[j]);
        }
        updateHighOrderDerivativesPhase2(predictedScaled, correctedScaled, predictedNordsieck);

        final FieldODEStateAndDerivative<T> updatedStepEnd =
                        equationsMapper.mapStateAndDerivative(globalCurrentState.getTime(), predictedY, correctedYDot);
        return new AdamsFieldStateInterpolator<>(getStepSize(), updatedStepEnd,
                                                          correctedScaled, predictedNordsieck, isForward,
                                                          getStepStart(), updatedStepEnd,
                                                          equationsMapper);

    }

    /** Corrector for current state in Adams-Moulton method.
     * <p>
     * This visitor implements the Taylor series formula:
     * <pre>
     * Y<sub>n+1</sub> = y<sub>n</sub> + s<sub>1</sub>(n+1) + [ -1 +1 -1 +1 ... &plusmn;1 ] r<sub>n+1</sub>
     * </pre>
     * </p>
     */
    private class Corrector implements FieldMatrixPreservingVisitor<T> {

        /** Previous state. */
        private final T[] previous;

        /** Current scaled first derivative. */
        private final T[] scaled;

        /** Current state before correction. */
        private final T[] before;

        /** Current state after correction. */
        private final T[] after;

        /** Simple constructor.
         * <p>
         * All arrays will be stored by reference to caller arrays.
         * </p>
         * @param previous previous state
         * @param scaled current scaled first derivative
         * @param state state to correct (will be overwritten after visit)
         */
        Corrector(final T[] previous, final T[] scaled, final T[] state) {
            this.previous = previous; // NOPMD - array reference storage is intentional and documented here
            this.scaled   = scaled;   // NOPMD - array reference storage is intentional and documented here
            this.after    = state;    // NOPMD - array reference storage is intentional and documented here
            this.before   = state.clone();
        }

        /** {@inheritDoc} */
        @Override
        public void start(int rows, int columns,
                          int startRow, int endRow, int startColumn, int endColumn) {
            Arrays.fill(after, getField().getZero());
        }

        /** {@inheritDoc} */
        @Override
        public void visit(int row, int column, T value) {
            if ((row & 0x1) == 0) {
                after[column] = after[column].subtract(value);
            } else {
                after[column] = after[column].add(value);
            }
        }

        /**
         * End visiting the Nordsieck vector.
         * <p>The correction is used to control stepsize. So its amplitude is
         * considered to be an error, which must be normalized according to
         * error control settings. If the normalized value is greater than 1,
         * the correction was too large and the step must be rejected.</p>
         * @return the normalized correction, if greater than 1, the step
         * must be rejected
         */
        @Override
        public T end() {

            final StepsizeHelper helper = getStepSizeHelper();
            T error = getField().getZero();
            for (int i = 0; i < after.length; ++i) {
                after[i] = after[i].add(previous[i].add(scaled[i]));
                if (i < helper.getMainSetDimension()) {
                    final T tol   = helper.getTolerance(i, MathUtils.max(previous[i].abs(), after[i].abs()));
                    final T ratio = after[i].subtract(before[i]).divide(tol); // (corrected-predicted)/tol
                    error = error.add(ratio.multiply(ratio));
                }
            }

            return error.divide(helper.getMainSetDimension()).sqrt();

        }
    }

}
