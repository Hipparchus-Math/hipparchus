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

import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealVector;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;

/** Constraint with lower and upper bounds: \(l \le f(x) \le u\).
 * @since 3.1
 */
public abstract class BoundedConstraint implements Constraint {

    /** Lower bound. */
    private final RealVector lower;

    /** Upper bound. */
    private final RealVector upper;


    /** Simple constructor.
     * <p>
     * At least one of the bounds must be non-null.
     * </p>
     * @param lower lower bound (null if no lower bound)
     * @param upper upper bound (null if no upper bound)
     */
     public BoundedConstraint(final RealVector lower, final RealVector upper) {

         // ensure lower is always properly set, even when there are no lower bounds
         if (lower == null) {
             MathUtils.checkNotNull(upper);
             this.lower = MatrixUtils.createRealVector(upper.getDimension());
             this.lower.set(Double.NEGATIVE_INFINITY);
         } else {
             this.lower = lower;
         }

         // ensure upper is always properly set, even when there are no upper bounds
         if (upper == null) {
             this.upper = MatrixUtils.createRealVector(lower.getDimension());
             this.upper.set(Double.POSITIVE_INFINITY);
         } else {
             this.upper = upper;
         }

         // safety check on dimensions
         MathUtils.checkDimension(this.lower.getDimension(), this.upper.getDimension());

     }

     /** {@inheritDoc} */
     @Override
     public int dimY() {
         return lower.getDimension();
     }

     /** {@inheritDoc} */
     @Override
     public RealVector getLowerBound() {
         return lower;
     }

     /** {@inheritDoc} */
     @Override
     public RealVector getUpperBound() {
         return upper;
     }

     /** {@inheritDoc} */
     @Override
     public double overshoot(final RealVector y) {

         double overshoot = 0;
         for (int i = 0; i < y.getDimension(); ++i) {
             overshoot += FastMath.max(0, lower.getEntry(i) - y.getEntry(i));
             overshoot += FastMath.max(0, y.getEntry(i) - upper.getEntry(i));
         }

         return overshoot;

     }

}
