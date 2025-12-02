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

import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;

public class RosenbrookConstraint extends InequalityConstraint {

    public RosenbrookConstraint(final double[] b) {
        super(MatrixUtils.createRealVector(b));
    }

     @Override
    public RealVector value(final RealVector x) {
        final double x0 = x.getEntry(0);
        final double x1 = x.getEntry(1);
        final RealVector a = new ArrayRealVector(5);
        a.setEntry(0, -x0 * x0 - x1 * x1);
        a.setEntry(1,  x0);
        a.setEntry(2, -x0);
        a.setEntry(3,  x1);
        a.setEntry(4, -x1);
        return a;
    }

     @Override
    public RealMatrix jacobian(final RealVector x) {
      final RealMatrix a = new Array2DRowRealMatrix(5, 2);
      a.setEntry(0, 0, -2 * x.getEntry(0));
      a.setEntry(0, 1, -2 * x.getEntry(1));
      a.setEntry(1, 0,  1);
      a.setEntry(2, 0, -1);
      a.setEntry(3, 1,  1);
      a.setEntry(4, 1, -1);
      return a;
    }

    @Override
    public int dim(){
         return 2;
     }

    @Override
    public int dimY(){
         return 5;
     }

}
