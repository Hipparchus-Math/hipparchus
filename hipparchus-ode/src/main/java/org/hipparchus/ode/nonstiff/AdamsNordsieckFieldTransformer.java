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
import java.util.HashMap;
import java.util.Map;

import org.hipparchus.Field;
import org.hipparchus.CalculusFieldElement;
import org.hipparchus.linear.Array2DRowFieldMatrix;
import org.hipparchus.linear.ArrayFieldVector;
import org.hipparchus.linear.FieldDecompositionSolver;
import org.hipparchus.linear.FieldLUDecomposition;
import org.hipparchus.linear.FieldMatrix;
import org.hipparchus.util.MathArrays;

/** Transformer to Nordsieck vectors for Adams integrators.
 * <p>This class is used by {@link AdamsBashforthIntegrator Adams-Bashforth} and
 * {@link AdamsMoultonIntegrator Adams-Moulton} integrators to convert between
 * classical representation with several previous first derivatives and Nordsieck
 * representation with higher order scaled derivatives.</p>
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
 * <p>With the previous definition, the classical representation of multistep methods
 * uses first derivatives only, i.e. it handles y<sub>n</sub>, s<sub>1</sub>(n) and
 * q<sub>n</sub> where q<sub>n</sub> is defined as:
 * \[
 *   q_n = [ s_1(n-1) s_1(n-2) \ldots s_1(n-(k-1)) ]^T
 * \]
 * (we omit the k index in the notation for clarity).</p>
 *
 * <p>Another possible representation uses the Nordsieck vector with
 * higher degrees scaled derivatives all taken at the same step, i.e it handles y<sub>n</sub>,
 * s<sub>1</sub>(n) and r<sub>n</sub>) where r<sub>n</sub> is defined as:
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
 * classical representation and Nordsieck vector at step end. The transform between r<sub>n</sub>
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
 * <p>Changing -i into +i in the formula above can be used to compute a similar transform between
 * classical representation and Nordsieck vector at step start. The resulting matrix is simply
 * the absolute value of matrix P.</p>
 *
 * <p>For {@link AdamsBashforthIntegrator Adams-Bashforth} method, the Nordsieck vector
 * at step n+1 is computed from the Nordsieck vector at step n as follows:
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
 * <p>For {@link AdamsMoultonIntegrator Adams-Moulton} method, the predicted Nordsieck vector
 * at step n+1 is computed from the Nordsieck vector at step n as follows:
 * <ul>
 *   <li>Y<sub>n+1</sub> = y<sub>n</sub> + s<sub>1</sub>(n) + u<sup>T</sup> r<sub>n</sub></li>
 *   <li>S<sub>1</sub>(n+1) = h f(t<sub>n+1</sub>, Y<sub>n+1</sub>)</li>
 *   <li>R<sub>n+1</sub> = (s<sub>1</sub>(n) - s<sub>1</sub>(n+1)) P<sup>-1</sup> u + P<sup>-1</sup> A P r<sub>n</sub></li>
 * </ul>
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
 * <p>We observe that both methods use similar update formulas. In both cases a P<sup>-1</sup>u
 * vector and a P<sup>-1</sup> A P matrix are used that do not depend on the state,
 * they only depend on k. This class handles these transformations.</p>
 *
 * @param <T> the type of the field elements
 */
public class AdamsNordsieckFieldTransformer<T extends CalculusFieldElement<T>> {

    /** Cache for already computed coefficients. */
    private static final Map<Integer,
                         Map<Field<? extends CalculusFieldElement<?>>,
                                   AdamsNordsieckFieldTransformer<? extends CalculusFieldElement<?>>>> CACHE = new HashMap<>();

    /** Field to which the time and state vector elements belong. */
    private final Field<T> field;

    /** Update matrix for the higher order derivatives h<sup>2</sup>/2 y'', h<sup>3</sup>/6 y''' ... */
    private final Array2DRowFieldMatrix<T> update;

    /** Update coefficients of the higher order derivatives wrt y'. */
    private final T[] c1;

    /** Simple constructor.
     * @param field field to which the time and state vector elements belong
     * @param n number of steps of the multistep method
     * (excluding the one being computed)
     */
    private AdamsNordsieckFieldTransformer(final Field<T> field, final int n) {

        this.field = field;
        final int rows = n - 1;

        // compute coefficients
        FieldMatrix<T> bigP = buildP(rows);
        FieldDecompositionSolver<T> pSolver =
                new FieldLUDecomposition<>(bigP).getSolver();

        T[] u = MathArrays.buildArray(field, rows);
        Arrays.fill(u, field.getOne());
        c1 = pSolver.solve(new ArrayFieldVector<>(u, false)).toArray();

        // update coefficients are computed by combining transform from
        // Nordsieck to multistep, then shifting rows to represent step advance
        // then applying inverse transform
        T[][] shiftedP = bigP.getData();
        for (int i = shiftedP.length - 1; i > 0; --i) {
            // shift rows
            shiftedP[i] = shiftedP[i - 1];
        }
        shiftedP[0] = MathArrays.buildArray(field, rows);
        Arrays.fill(shiftedP[0], field.getZero());
        update = new Array2DRowFieldMatrix<>(pSolver.solve(new Array2DRowFieldMatrix<>(shiftedP, false)).getData());

    }

    /** Get the Nordsieck transformer for a given field and number of steps.
     * @param field field to which the time and state vector elements belong
     * @param nSteps number of steps of the multistep method
     * (excluding the one being computed)
     * @return Nordsieck transformer for the specified field and number of steps
     * @param <T> the type of the field elements
     */
    public static <T extends CalculusFieldElement<T>> AdamsNordsieckFieldTransformer<T>
    getInstance(final Field<T> field, final int nSteps) { // NOPMD - PMD false positive
        synchronized(CACHE) {
            Map<Field<? extends CalculusFieldElement<?>>,
                    AdamsNordsieckFieldTransformer<? extends CalculusFieldElement<?>>> map = CACHE.computeIfAbsent(nSteps, k -> new HashMap<>());
            @SuppressWarnings("unchecked")
            AdamsNordsieckFieldTransformer<T> t = (AdamsNordsieckFieldTransformer<T>) map.get(field);
            if (t == null) {
                t = new AdamsNordsieckFieldTransformer<>(field, nSteps);
                map.put(field, t);
            }
            return t;

        }
    }

    /** Build the P matrix.
     * <p>The P matrix general terms are shifted \((j+1) (-i)^j\) terms
     * with i being the row number starting from 1 and j being the column
     * number starting from 1:
     * <pre>
     *        [  -2   3   -4    5  ... ]
     *        [  -4  12  -32   80  ... ]
     *   P =  [  -6  27 -108  405  ... ]
     *        [  -8  48 -256 1280  ... ]
     *        [          ...           ]
     * </pre></p>
     * @param rows number of rows of the matrix
     * @return P matrix
     */
    private FieldMatrix<T> buildP(final int rows) {

        final T[][] pData = MathArrays.buildArray(field, rows, rows);

        for (int i = 1; i <= pData.length; ++i) {
            // build the P matrix elements from Taylor series formulas
            final T[] pI = pData[i - 1];
            final int factor = -i;
            T aj = field.getZero().add(factor);
            for (int j = 1; j <= pI.length; ++j) {
                pI[j - 1] = aj.multiply(j + 1);
                aj = aj.multiply(factor);
            }
        }

        return new Array2DRowFieldMatrix<>(pData, false);

    }

    /** Initialize the high order scaled derivatives at step start.
     * @param h step size to use for scaling
     * @param t first steps times
     * @param y first steps states
     * @param yDot first steps derivatives
     * @return Nordieck vector at start of first step (h<sup>2</sup>/2 y''<sub>n</sub>,
     * h<sup>3</sup>/6 y'''<sub>n</sub> ... h<sup>k</sup>/k! y<sup>(k)</sup><sub>n</sub>)
     */

    public Array2DRowFieldMatrix<T> initializeHighOrderDerivatives(final T h, final T[] t,
                                                                   final T[][] y,
                                                                   final T[][] yDot) {

        // using Taylor series with di = ti - t0, we get:
        //  y(ti)  - y(t0)  - di y'(t0) =   di^2 / h^2 s2 + ... +   di^k     / h^k sk + O(h^k)
        //  y'(ti) - y'(t0)             = 2 di   / h^2 s2 + ... + k di^(k-1) / h^k sk + O(h^(k-1))
        // we write these relations for i = 1 to i= 1+n/2 as a set of n + 2 linear
        // equations depending on the Nordsieck vector [s2 ... sk rk], so s2 to sk correspond
        // to the appropriately truncated Taylor expansion, and rk is the Taylor remainder.
        // The goal is to have s2 to sk as accurate as possible considering the fact the sum is
        // truncated and we don't want the error terms to be included in s2 ... sk, so we need
        // to solve also for the remainder
        final T[][] a     = MathArrays.buildArray(field, c1.length + 1, c1.length + 1);
        final T[][] b     = MathArrays.buildArray(field, c1.length + 1, y[0].length);
        final T[]   y0    = y[0];
        final T[]   yDot0 = yDot[0];
        for (int i = 1; i < y.length; ++i) {

            final T di    = t[i].subtract(t[0]);
            final T ratio = di.divide(h);
            T dikM1Ohk    = h.reciprocal();

            // linear coefficients of equations
            // y(ti) - y(t0) - di y'(t0) and y'(ti) - y'(t0)
            final T[] aI    = a[2 * i - 2];
            final T[] aDotI = (2 * i - 1) < a.length ? a[2 * i - 1] : null;
            for (int j = 0; j < aI.length; ++j) {
                dikM1Ohk = dikM1Ohk.multiply(ratio);
                aI[j]    = di.multiply(dikM1Ohk);
                if (aDotI != null) {
                    aDotI[j]  = dikM1Ohk.multiply(j + 2);
                }
            }

            // expected value of the previous equations
            final T[] yI    = y[i];
            final T[] yDotI = yDot[i];
            final T[] bI    = b[2 * i - 2];
            final T[] bDotI = (2 * i - 1) < b.length ? b[2 * i - 1] : null;
            for (int j = 0; j < yI.length; ++j) {
                bI[j]    = yI[j].subtract(y0[j]).subtract(di.multiply(yDot0[j]));
                if (bDotI != null) {
                    bDotI[j] = yDotI[j].subtract(yDot0[j]);
                }
            }

        }

        // solve the linear system to get the best estimate of the Nordsieck vector [s2 ... sk],
        // with the additional terms s(k+1) and c grabbing the parts after the truncated Taylor expansion
        final FieldLUDecomposition<T> decomposition = new FieldLUDecomposition<>(new Array2DRowFieldMatrix<>(a, false));
        final FieldMatrix<T> x = decomposition.getSolver().solve(new Array2DRowFieldMatrix<>(b, false));

        // extract just the Nordsieck vector [s2 ... sk]
        final Array2DRowFieldMatrix<T> truncatedX =
                        new Array2DRowFieldMatrix<>(field, x.getRowDimension() - 1, x.getColumnDimension());
        for (int i = 0; i < truncatedX.getRowDimension(); ++i) {
            for (int j = 0; j < truncatedX.getColumnDimension(); ++j) {
                truncatedX.setEntry(i, j, x.getEntry(i, j));
            }
        }
        return truncatedX;

    }

    /** Update the high order scaled derivatives for Adams integrators (phase 1).
     * <p>The complete update of high order derivatives has a form similar to:
     * \[
     * r_{n+1} = (s_1(n) - s_1(n+1)) P^{-1} u + P^{-1} A P r_n
     * \]
     * this method computes the P<sup>-1</sup> A P r<sub>n</sub> part.</p>
     * @param highOrder high order scaled derivatives
     * (h<sup>2</sup>/2 y'', ... h<sup>k</sup>/k! y(k))
     * @return updated high order derivatives
     * @see #updateHighOrderDerivativesPhase2(CalculusFieldElement[], CalculusFieldElement[], Array2DRowFieldMatrix)
     */
    public Array2DRowFieldMatrix<T> updateHighOrderDerivativesPhase1(final Array2DRowFieldMatrix<T> highOrder) {
        return update.multiply(highOrder);
    }

    /** Update the high order scaled derivatives Adams integrators (phase 2).
     * <p>The complete update of high order derivatives has a form similar to:
     * \[
     * r_{n+1} = (s_1(n) - s_1(n+1)) P^{-1} u + P^{-1} A P r_n
     * \]
     * this method computes the (s<sub>1</sub>(n) - s<sub>1</sub>(n+1)) P<sup>-1</sup> u part.</p>
     * <p>Phase 1 of the update must already have been performed.</p>
     * @param start first order scaled derivatives at step start
     * @param end first order scaled derivatives at step end
     * @param highOrder high order scaled derivatives, will be modified
     * (h<sup>2</sup>/2 y'', ... h<sup>k</sup>/k! y(k))
     * @see #updateHighOrderDerivativesPhase1(Array2DRowFieldMatrix)
     */
    public void updateHighOrderDerivativesPhase2(final T[] start,
                                                 final T[] end,
                                                 final Array2DRowFieldMatrix<T> highOrder) {
        final T[][] data = highOrder.getDataRef();
        for (int i = 0; i < data.length; ++i) {
            final T[] dataI = data[i];
            final T c1I = c1[i];
            for (int j = 0; j < dataI.length; ++j) {
                dataI[j] = dataI[j].add(c1I.multiply(start[j].subtract(end[j])));
            }
        }
    }

}
