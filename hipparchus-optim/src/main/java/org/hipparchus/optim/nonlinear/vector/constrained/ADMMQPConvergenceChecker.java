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
import org.hipparchus.optim.ConvergenceChecker;
import org.hipparchus.optim.OptimizationData;
import org.hipparchus.util.FastMath;

/** Convergence Checker for ADMM QP Optimizer.
 * @since 3.1
 */
public class ADMMQPConvergenceChecker implements ConvergenceChecker<LagrangeSolution>, OptimizationData  {

    private  RealMatrix H;
    private  RealMatrix A;
    private  RealVector q;
    private  double epsAbs;
    private  double epsRel;
    private  boolean coverged = false;
     ADMMQPConvergenceChecker(RealMatrix H,RealMatrix A,RealVector q,double epsAbs,double epsRel)
     {
       this.H = H;
       this.A = A;
       this.q = q;
       this.epsAbs = epsAbs;
       this.epsRel = epsRel;
     }


    @Override
    public boolean converged(int i, LagrangeSolution previous, LagrangeSolution current) {
        return this.coverged;
    }

    //TO AVOID DOUBLE CALCULATION
    public boolean converged(int iteration,double rp,double rd,double maxPrimal,double maxDual)
    {   boolean result = false;


        //if (rp<=epsPrimalDual(maxPrimal))
       if (rp<=epsPrimalDual(maxPrimal) && rd<=epsPrimalDual(maxDual))
       {     result = true;
         this.coverged = true;
       }
         return result;
    }

     public double residualPrime(RealVector x, RealVector z) {
        return ((A.operate(x)).subtract(z)).getLInfNorm();
    }

     public double residualDual(RealVector x, RealVector y) {

       return (q.add(A.transpose().operate(y)).add(H.operate(x))).getLInfNorm();
    }

    public double maxPrimal(RealVector x, RealVector z) {
       return FastMath.max((A.operate(x)).getLInfNorm(),z.getLInfNorm());
    }

    public double maxDual(RealVector x, RealVector y) {

        return FastMath.max(FastMath.max((H.operate(x)).getLInfNorm(),(A.transpose().operate(y)).getLInfNorm()),q.getLInfNorm());
    }

    private double epsPrimalDual(double maxPrimalDual) {
        return epsAbs + epsRel * maxPrimalDual;
    }


}
