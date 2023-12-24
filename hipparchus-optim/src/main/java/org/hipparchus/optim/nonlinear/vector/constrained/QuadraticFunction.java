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

/**
 * Given A, b, c, implements 0.5*(x^T)A(x) + b.x + c;
 * The gradient is A(x) + b, and the Hessian is A
 * @since 3.1
 */
public class QuadraticFunction extends TwiceDifferentiableFunction {
    private final RealMatrix A;
    private final RealVector b;
    private final double c;
    private final int n;

    /**
     * Construct quadratic function 0.5*(x^T)A(x) + b.x + c
     *
     * @param A square matrix of weights for quadratic terms.
     * Typically expected to be positive definite or positive semi-definite.
     * @param b vector of weights for linear terms.
     * @param c a constant
     */
    public QuadraticFunction(RealMatrix A, RealVector b, double c) {
        int d = b.getDimension();
        if (d < 1) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.INSUFFICIENT_DIMENSION, d, 1);
        }
        MathUtils.checkDimension(A.getRowDimension(), d);
       // MatrixUtils.checkSymmetric(A, 1e-6);
        this.A = A.copy();
        this.b = b.copy();
        this.c = c;
        this.n = d;
    }

    /**
     * Construct quadratic function 0.5*(x^T)A(x) + b.x + c
     *
     * @param A square matrix of weights for quadratic terms.
     * Typically expected to be positive definite or positive semi-definite.
     * @param b vector of weights for linear terms.
     * @param c a constant
     */
    public QuadraticFunction(double[][] A, double[] b, double c) {
        this(new Array2DRowRealMatrix(A), new ArrayRealVector(b), c);
    }

    public RealMatrix getH()
    {
        return this.A;
    }
     public RealVector getC()
    {
        return this.b;
    }

      public double getD()
    {
        return this.c;
    }
    @Override
    public int dim() { return n; }

    @Override
    public double value(final RealVector x) {
        double v = 0.5 * (A.operate(x)).dotProduct(x);
        v += b.dotProduct(x);
        v += c;
        return v;
    }

    @Override
    public RealVector gradient(final RealVector x) {
        return (A.operate(x)).add(b);
    }

    @Override
    public RealMatrix hessian(final RealVector x) {
        return A.copy();
    }


}
