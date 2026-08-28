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
package org.hipparchus.optim.nonlinear.vector.constrained.sqp.hock_shittkowski;

import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.vector.constrained.GradientMode;
import org.hipparchus.optim.nonlinear.vector.constrained.LagrangeSolution;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOptimizerS2;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPOption;
import org.hipparchus.optim.nonlinear.vector.constrained.SQPProblem;
import org.junit.jupiter.api.Test;

public class HS001TestProblem {
    
    private static class HS001Problem implements SQPProblem{

        @Override
        public int getdim() {
            return 2; 
        }

       

        @Override
        public double[] getBoundsLB() {
            return null;
        }

        @Override
        public double[] getBoundsUB() {
            return null;
        }

        @Override
        public double getObjectiveEvaluation(RealVector x) {
            
            final double x0     = x.getEntry(0);
            final double x1     = x.getEntry(1);
            final double x1Mx02 = x1 - x0  * x0;
            final double oMx0   = 1 - x0;
            return 100 * x1Mx02 * x1Mx02 + oMx0 * oMx0;
        }

        @Override
        public RealVector getObjectiveGradient(RealVector x) {
             final double x0     = x.getEntry(0);
            final double x1     = x.getEntry(1);
            final double x1Mx02 = x1 - x0  * x0;
            final double a      = 200 * x1Mx02;
            return MatrixUtils.createRealVector(new double[] {
                    2 * (x0 * (1 - a) - 1),
                    a
            });
        }

        @Override
        public RealVector getEqCostraintEvaluation(RealVector rv) {
            return null;
        }

        @Override
        public RealMatrix getEqCostraintJacobian(RealVector rv) {
             return null;
        }

        @Override
        public RealVector getEqCostraintLB() {
             return null;
        }

        @Override
        public RealVector getIneqConstraintEvaluation(RealVector rv) {
             return null;
        }

        @Override
        public RealMatrix getIneqCostraintJacobian(RealVector rv) {
             return null;
        }

        @Override
        public RealVector getIneqCostraintLB() {
             return null;
        }

       

        @Override
        public boolean hasBounds() {
            return false;
        }

        @Override
        public boolean hasEquality() {
            return false;
        }

        @Override
        public boolean hasInequality() {
           return false;
        }
        
    }

   
    
    @Test
    public void testHS001ExternalGradient() {
        doTestHS001(GradientMode.EXTERNAL);
    }

    @Test
    public void testHS001ForwardGradient() {
        doTestHS001(GradientMode.FORWARD);
    }

    @Test
    public void testHS001CentralGradient() {
        doTestHS001(GradientMode.CENTRAL);
    }

    private void doTestHS001(final GradientMode gradientMode) {
        SQPOption sqpOption=new SQPOption();
        sqpOption.setGradientMode(gradientMode);
        SQPOptimizerS2 optimizer = HSProblemTestUtils.newOptimizer();
        double val = 0.0;
         InitialGuess guess = new InitialGuess(new double[]{-2, 1});
        LagrangeSolution sol = optimizer.optimize(sqpOption,guess,new HS001Problem());
        
        HSProblemTestUtils.assertExpectedObjective(val, sol);
    }
}
