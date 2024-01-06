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

public class RosenbrockFunction extends TwiceDifferentiableFunction {

    @Override
    public int dim() {
        return 2;
    }

    @Override
    public double value(RealVector x) {
        return (1.0-x.getEntry(0))*(1.0-x.getEntry(0))+100.0*(x.getEntry(1)-x.getEntry(0)*x.getEntry(0))*(x.getEntry(1)-x.getEntry(0)*x.getEntry(0));
    }

    @Override
    public RealVector gradient(RealVector x) {
        RealVector grad=new ArrayRealVector(x.getDimension());
        grad.setEntry(0,-2.0*(1.0-x.getEntry(0))-400.0*(x.getEntry(1)-x.getEntry(0)*x.getEntry(0))*x.getEntry(0));
        grad.setEntry(1,200.0*(x.getEntry(1)-x.getEntry(0)*x.getEntry(0)));
        return grad;
    }

    @Override
    public RealMatrix hessian(RealVector x) {
        RealMatrix h=new Array2DRowRealMatrix(2,2);
        double a00=800.0*x.getEntry(0)*x.getEntry(0)-400.0*(x.getEntry(1)-x.getEntry(0)*x.getEntry(0))+2;
        double a01=-400.0*x.getEntry(0);
        double a10=-400*x.getEntry(0);
        double a11=200;
        h.setEntry(0,0, a00);
        h.setEntry(0,1, a01);
        h.setEntry(1,0, a10);
        h.setEntry(1,1, a11);
        return h;
    }

}
