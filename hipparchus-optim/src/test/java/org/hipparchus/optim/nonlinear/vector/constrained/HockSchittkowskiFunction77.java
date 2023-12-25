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
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;

public class HockSchittkowskiFunction77 extends TwiceDifferentiableFunction {

    @Override
    public int dim() {
        return 5;
    }

    @Override
    public double value(RealVector x) {
        double x1 = x.getEntry(0);
        double x2 = x.getEntry(1);
        double x3 = x.getEntry(2);
        double x4 = x.getEntry(3);
        double x5 = x.getEntry(4);
        // obj = (x[1]-1)^2 + (x[1] - x[2])^2 + (x[3]-1)^2 + (x[4]-1)^4 + (x[5]-1)^6
        return (x1-1)*(x1-1) + (x1 - x2)*(x1 - x2) + (x3-1)*(x3-1) + (x4-1)*(x4-1)*(x4-1)*(x4-1) + (x5-1)*(x5-1)*(x5-1)*(x5-1)*(x5-1)*(x5-1);
    }

    @Override
    public RealVector gradient(RealVector x) {
        RealVector grad=new ArrayRealVector(x.getDimension());
        double x1 = x.getEntry(0);
        double x2 = x.getEntry(1);
        double x3 = x.getEntry(2);
        double x4 = x.getEntry(3);
        double x5 = x.getEntry(4);
        grad.setEntry(0,2.0*(x1-1.0)+2.0*(x1-x2));
        grad.setEntry(1,-2.0*(x1-x2));
        grad.setEntry(2,2.0*(x3-1));
         grad.setEntry(3,4.0*(x4-1)*(x4-1)*(x4-1));
          grad.setEntry(4,6.0*(x5-1)*(x5-1)*(x5-1)*(x5-1)*(x5-1));
        return grad;
    }

    @Override
    public RealMatrix hessian(RealVector x) {
        RealMatrix h=new Array2DRowRealMatrix(x.getDimension(),x.getDimension());
        double x4 = x.getEntry(3);
        double x5 = x.getEntry(4);
        h.setEntry(0,0, 4.0);
        h.setEntry(0,1, -2.0);
        h.setEntry(0,2, 0);
        h.setEntry(0,3, 0);
        h.setEntry(0,4, 0);

         h.setEntry(1,0, -2.0);
        h.setEntry(1,1, 2.0);
        h.setEntry(1,2, 0);
        h.setEntry(1,3, 0);
        h.setEntry(1,4, 0);

        h.setEntry(2,0, 0);
        h.setEntry(2,1, 0);
        h.setEntry(2,2, 2.0);
        h.setEntry(2,3, 0);
        h.setEntry(2,4, 0);

        h.setEntry(3,0, 0);
        h.setEntry(3,1, 0);
        h.setEntry(3,2, 0.0);
        h.setEntry(3,3, 12.0*(x4-1)*(x4-1));
        h.setEntry(3,4, 0);

        h.setEntry(4,0, 0);
        h.setEntry(4,1, 0);
        h.setEntry(4,2, 0.0);
        h.setEntry(4,3, 0.0);
        h.setEntry(4,4, 30.0*(x5-1)*(x5-1)*(x5-1)*(x5-1));

        return h;
    }

}
