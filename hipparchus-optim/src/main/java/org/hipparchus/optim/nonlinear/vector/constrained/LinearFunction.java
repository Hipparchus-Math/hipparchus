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

import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.linear.OpenMapRealMatrix;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.linear.ArrayRealVector;
/** Linear function b.x + c; The gradient is b, and the Hessian is |0| (all zeros).
 * @since 3.1
 */
public class LinearFunction extends TwiceDifferentiableFunction {
    private final RealVector b;
    private final double c;
    private final int n;

    /**
     * Construct a linear function b.x + c
     *
     * @param b a weight vector
     * @param c a constant
     */
    public LinearFunction(RealVector b, double c) {
        int d = b.getDimension();
        if (d < 1) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.INSUFFICIENT_DIMENSION, d, 1);
        }
        this.b = b;
        this.c = c;
        this.n = d;
    }

    /**
     * Construct a linear function b.x + c
     *
     * @param b a weight vector
     * @param c a constant
     */
    public LinearFunction(double[] b, double c) {
        this(new ArrayRealVector(b), c);
    }

    @Override
    public int dim() { return n; }

    @Override
    public double value(final RealVector x) {
        return c + b.dotProduct(x);
    }

    @Override
    public RealVector gradient(final RealVector x) {
        return b.copy();
    }

    @Override
    public RealMatrix hessian(final RealVector x) {
        // the Hessian is just zero for a linear function
        return new OpenMapRealMatrix(n, n);
    }

    @Override
    public String toString() {
        return String.format("LinearFunction(%g, %s)", c, b.toString());
    }
}
