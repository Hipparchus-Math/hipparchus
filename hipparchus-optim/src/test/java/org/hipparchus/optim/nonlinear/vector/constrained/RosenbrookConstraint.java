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

public class RosenbrookConstraint extends LinearInequalityConstraint {

    public RosenbrookConstraint(RealMatrix m, RealVector b) {
        super(m, b);
    }

     @Override
    public RealVector value(RealVector x) {
        RealVector a=new ArrayRealVector(5);
       a.setEntry(0,-x.getEntry(0)*x.getEntry(0)-x.getEntry(1)*x.getEntry(1));
       a.setEntry(1,x.getEntry(0));
       a.setEntry(2,-x.getEntry(0));
        a.setEntry(3,x.getEntry(1));
       a.setEntry(4,-x.getEntry(1));

         return a;
    }

     @Override
    public RealMatrix jacobian(RealVector x) {
      RealMatrix a= new Array2DRowRealMatrix(5,2);
      a.setEntry(0, 0,-2*x.getEntry(0));
      a.setEntry(0, 1,-2*x.getEntry(1));
      a.setEntry(1, 0,1);
       a.setEntry(2, 0,-1);
        a.setEntry(3, 1,1);
          a.setEntry(4, 1,-1);
        return a;
    }
//
//    @Override
//      public RealVector getLowerBound()
//    {   RealVector lb=new ArrayRealVector(9);
//        lb.setEntry(0,25.0);
//        lb.setEntry(1,1.0);
//        lb.setEntry(2,1.0);
//        lb.setEntry(3,1.0);
//        lb.setEntry(4,1.0);
//        lb.setEntry(5,-5.0);
//        lb.setEntry(6,-5.0);
//        lb.setEntry(7,-5.0);
//        lb.setEntry(8,-5.0);
//        return lb;
//    }
//
//    /**
// * Return Upper Bound .
// * @return Upper Bound Vector
// */
//      @Override
//     public RealVector getUpperBound()
//    {
//        return new ArrayRealVector(9,Double.POSITIVE_INFINITY);
//    }
//      @Override
     public int dimY(){
         return 5;
     }
}
