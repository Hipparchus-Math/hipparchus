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

/** Given h, c, d, implements 0.5*(x^T)h(x) + c.x + d.
 * The gradient is h(x) + c, and the Hessian is h
 * @since 3.1
 */
public class QuadraticFunction extends TwiceDifferentiableFunction {

    /** Square matrix of weights for quadratic terms. */
    private final RealMatrix h;

    /** Vector of weights for linear terms. */
    private final RealVector c;

    /** Constant term. */
    private final double d;

    /** Dimension of the vector. */
    private final int n;

    /** Construct quadratic function 0.5*(x^T)h(x) + c.x + d.
     * @param h square matrix of weights for quadratic terms.
     * Typically expected to be positive definite or positive semi-definite.
     * @param c vector of weights for linear terms.
     * @param d constant term
     */
    public QuadraticFunction(RealMatrix h, RealVector c, double d) {
        this.n = c.getDimension();
        if (n < 1) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.INSUFFICIENT_DIMENSION, d, 1);
        }
        MathUtils.checkDimension(h.getRowDimension(), n);
        this.h = h.copy();
        this.c = c.copy();
        this.d = d;
   }

    /** Construct quadratic function 0.5*(x^T)h(x) + c.x + d.
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
    public RealMatrix getH() {
        return this.h;
    }

    /** Get vector of weights for linear terms.
     * @return vector of weights for linear terms
     */
    public RealVector getC() {
        return this.c;
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
        return 0.5 * h.operate(x).dotProduct(x) + c.dotProduct(x) + d;
    }

    /** {@inheritDoc} */
    @Override
    public RealVector gradient(final RealVector x) {
        return h.operate(x).add(c);
    }

    /** {@inheritDoc} */
    @Override
    public RealMatrix hessian(final RealVector x) {
        return h.copy();
    }

}
