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

public class HockSchittkowskiConstraintInequality72 extends InequalityConstraint {
    double a11 = 4;
    double a12 = 2.25;
    double a13 = 1.0;
    double a14 = 0.25;
    double a21 = 0.16;
    double a22 = 0.36;
    double a23 = 0.64;
    double a24 = 0.64;
    double b1 = 0.0401;
    double b2=0.010085;

    public HockSchittkowskiConstraintInequality72() {

    }

     @Override
    public RealVector value(RealVector x) {
        double x1 = x.getEntry(0);
        double x2 = x.getEntry(1);
        double x3 = x.getEntry(2);
        double x4 = x.getEntry(3);

       // a[1][1]/x[1] + a[1][2]/x[2] + a[1][3]/x[3] + a[1][4]/x[4] <= b[1]
    //a[2][1]/x[1] + a[2][2]/x[2] + a[2][3]/x[3] + a[2][4]/x[4] <= b[2]
        RealVector a=new ArrayRealVector(2);
        a.setEntry(0,-a11/x1-a12/x2-a13/x3-a14/x4);
         a.setEntry(1,-a21/x1-a22/x2-a23/x3-a24/x4);


         return a;
    }

     @Override
    public RealMatrix jacobian(RealVector x) {
      double x1 = x.getEntry(0);
      double x2 = x.getEntry(1);
      double x3 = x.getEntry(2);
      double x4 = x.getEntry(3);

      RealMatrix a= new Array2DRowRealMatrix(2,4);

      a.setEntry(0,0,a11/(x1*x1));
      a.setEntry(0,1,a12/(x2*x2));
      a.setEntry(0,2,a13/(x3*x3));
      a.setEntry(0,3,a14/(x4*x4));

      a.setEntry(1,0,a21/(x1*x1));
      a.setEntry(1,1,a22/(x2*x2));
      a.setEntry(1,2,a23/(x3*x3));
      a.setEntry(1,3,a24/(x4*x4));


        return a;
    }
    @Override
      public RealVector getLowerBound()
    {   RealVector lb=new ArrayRealVector(2);
        lb.setEntry(0,-b1);
        lb.setEntry(1,-b2);

        return lb;
    }

    /**
 * Return Upper Bound .
 * @return Upper Bound Vector
 */
      @Override
     public RealVector getUpperBound()
    {   RealVector lb=new ArrayRealVector(2);
        lb.setEntry(0,Double.POSITIVE_INFINITY);
        lb.setEntry(1,Double.POSITIVE_INFINITY);

        return lb;
    }

      @Override
     public int dimY(){
         return 2;
     }
      @Override
     public int dim(){
         return 4;
     }
}
