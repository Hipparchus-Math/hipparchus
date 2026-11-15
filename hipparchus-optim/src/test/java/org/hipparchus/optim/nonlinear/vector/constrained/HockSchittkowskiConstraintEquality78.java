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

public class HockSchittkowskiConstraintEquality78 extends EqualityConstraint {

    public HockSchittkowskiConstraintEquality78() {
        super(MatrixUtils.createRealVector(new double[] {
            10, 0, -1
        }));
    }

    @Override
    public RealVector value(RealVector x) {
        double x1 = x.getEntry(0);
        double x2 = x.getEntry(1);
        double x3 = x.getEntry(2);
        double x4 = x.getEntry(3);
        double x5 = x.getEntry(4);
        // x[1]^2+x[2]^2+x[3]^2+x[4]^2+x[5]^2 = 10
        //x[2]*x[3] - 5*x[4]*x[5] = 0
        //x[1]^3 + x[2]^3 = -1
        RealVector a=new ArrayRealVector(3);
        a.setEntry(0,(x1*x1)+(x2*x2)+(x3*x3)+(x4*x4)+x5*x5);
        a.setEntry(1,x2*x3-5.0*x4*x5);
        a.setEntry(2,x1*x1*x1+x2*x2*x2);



        return a;
    }

    @Override
    public RealMatrix jacobian(RealVector x) {
        double x1 = x.getEntry(0);
        double x2 = x.getEntry(1);
        double x3 = x.getEntry(2);
        double x4 = x.getEntry(3);
        double x5 = x.getEntry(4);
        RealMatrix a= new Array2DRowRealMatrix(3,5);

        a.setEntry(0,0,2.0*x1);
        a.setEntry(0,1,2.0*x2);
        a.setEntry(0,2,2.0*x3);
        a.setEntry(0,3,2.0*x4);
        a.setEntry(0,4,2.0*x5);

        a.setEntry(1,0,0.0);
        a.setEntry(1,1,x3);
        a.setEntry(1,2,x2);
        a.setEntry(1,3,-5.0*x5);
        a.setEntry(1,4,-5.0*x4);

        a.setEntry(2,0,3.0*x1*x1);
        a.setEntry(2,1,3.0*x2*x2);


        return a;
    }

    @Override
    public int dim(){
        return 5;
    }

}
