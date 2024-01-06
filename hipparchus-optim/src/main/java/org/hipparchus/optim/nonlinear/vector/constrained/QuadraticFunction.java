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
package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.util.MathUtils;

/** Given h, c, d, implements \(\frac{1}{2}x^T P X + Q^T x + d\).
 * The gradient is P x + Q^T, and the Hessian is P
 * @since 3.1
 */
public class QuadraticFunction extends TwiceDifferentiableFunction {

    /** Square matrix of weights for quadratic terms. */
    private final RealMatrix p;

    /** Vector of weights for linear terms. */
    private final RealVector q;

    /** Constant term. */
    private final double d;

    /** Dimension of the vector. */
    private final int n;

    /** Construct quadratic function \(\frac{1}{2}x^T P X + Q^T x + d\).
     * @param p square matrix of weights for quadratic terms.
     * Typically expected to be positive definite or positive semi-definite.
     * @param q vector of weights for linear terms.
     * @param d constant term
     */
    public QuadraticFunction(RealMatrix p, RealVector q, double d) {
        this.n = q.getDimension();
        if (n < 1) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.INSUFFICIENT_DIMENSION, d, 1);
        }
        MathUtils.checkDimension(p.getRowDimension(), n);
        this.p = p.copy();
        this.q = q.copy();
        this.d = d;
   }

    /** Construct quadratic function \(\frac{1}{2}x^T P X + Q^T x + d\).
     * @param h square matrix of weights for quadratic terms.
     * Typically expected to be positive definite or positive semi-definite.
     * @param c vector of weights for linear terms.
     * @param d constant term
     */
    public QuadraticFunction(double[][] h, double[] c, double d) {
        this(new Array2DRowRealMatrix(h), new ArrayRealVector(c), d);
    }

    /** Get square matrix of weights for quadratic terms.
     * @return square matrix of weights for quadratic terms
     */
    public RealMatrix getP() {
        return this.p;
    }

    /** Get vector of weights for linear terms.
     * @return vector of weights for linear terms
     */
    public RealVector getQ() {
        return this.q;
    }

    /** Get constant term.
     * @return constant term
     */
    public double getD() {
        return d;
    }

    /** {@inheritDoc} */
    @Override
    public int dim() {
        return n;
    }

    /** {@inheritDoc} */
    @Override
    public double value(final RealVector x) {
        return 0.5 * p.operate(x).dotProduct(x) + q.dotProduct(x) + d;
    }

    /** {@inheritDoc} */
    @Override
    public RealVector gradient(final RealVector x) {
        return p.operate(x).add(q);
    }

    /** {@inheritDoc} */
    @Override
    public RealMatrix hessian(final RealVector x) {
        return p.copy();
    }

}
