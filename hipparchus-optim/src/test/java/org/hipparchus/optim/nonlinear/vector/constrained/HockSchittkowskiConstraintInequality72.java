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

public class HockSchittkowskiConstraintInequality72 extends InequalityConstraint {

    private static final double A11 = 4;
    private static final double A12 = 2.25;
    private static final double A13 = 1.0;
    private static final double A14 = 0.25;
    private static final double A21 = 0.16;
    private static final double A22 = 0.36;
    private static final double A23 = 0.64;
    private static final double A24 = 0.64;
    private static final double B1  = 0.0401;
    private static final double B2  = 0.010085;

    public HockSchittkowskiConstraintInequality72() {
        super(MatrixUtils.createRealVector(new double[] { -B1, -B2 }));
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
        a.setEntry(0,-A11/x1-A12/x2-A13/x3-A14/x4);
        a.setEntry(1,-A21/x1-A22/x2-A23/x3-A24/x4);


        return a;
    }

    @Override
    public RealMatrix jacobian(RealVector x) {
        double x1 = x.getEntry(0);
        double x2 = x.getEntry(1);
        double x3 = x.getEntry(2);
        double x4 = x.getEntry(3);

        RealMatrix a= new Array2DRowRealMatrix(2,4);

        a.setEntry(0,0,A11/(x1*x1));
        a.setEntry(0,1,A12/(x2*x2));
        a.setEntry(0,2,A13/(x3*x3));
        a.setEntry(0,3,A14/(x4*x4));

        a.setEntry(1,0,A21/(x1*x1));
        a.setEntry(1,1,A22/(x2*x2));
        a.setEntry(1,2,A23/(x3*x3));
        a.setEntry(1,3,A24/(x4*x4));


        return a;
    }

    @Override
    public int dim() {
        return 4;
    }

}
