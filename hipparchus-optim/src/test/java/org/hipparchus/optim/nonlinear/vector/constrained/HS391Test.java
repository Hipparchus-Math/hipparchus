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

import org.hipparchus.linear.RealVector;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.optim.InitialGuess;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;


public class HS391Test {

    
    static final class HS391Obj extends TwiceDifferentiableFunction {

        @Override
        public int dim() { return 30; }

        @Override
        public double value(RealVector x) {
            double fx = 0.0;

            for (int i = 0; i < 30; i++) {
                int I = i + 1; // Fortran indexing

                double sum = 0.0;
                for (int j = 0; j < 30; j++) {
                    if (j == i) continue;

                    int J = j + 1;
                    double wurz = Math.sqrt(x.getEntry(j) * x.getEntry(j) + (double) I / (double) J);
                    double t = Math.log(wurz);
                    double part = Math.sin(t);
                    double part2 = Math.cos(t);
                    sum += wurz * (Math.pow(part, 5) + Math.pow(part2, 5));
                }

                double term = 420.0 * x.getEntry(i)
                             + Math.pow((I - 15), 3)
                             + sum;

                fx += term * term;
            }

            return fx;
        }

        
        @Override
        public RealVector gradient(RealVector x) {
            throw new UnsupportedOperationException("Gradient not provided for TP391");
        }

        @Override
        public RealMatrix hessian(RealVector x) {
            throw new UnsupportedOperationException("Hessian not provided");
        }
    }


   
    private double[] initialPoint() {
        double[] x = new double[30];

        for (int i = 1; i <= 30; i++) {
            double sum = 0.0;

            for (int j = 1; j <= 30; j++) {
                if (j == i) continue;

                double wurz = Math.sqrt((double) i / (double) j);
                double t = Math.log(wurz);
                sum += wurz * (Math.pow(Math.sin(t), 5) + Math.pow(Math.cos(t), 5));
            }

            x[i - 1] = -2.8742711 * ((double) ((i - 15) * (i - 15) * (i - 15)) + sum);
        }

        return x;
    }


    private LagrangeSolution solve() {
        SQPOptimizerS2 opt = new SQPOptimizerS2();
         if (Boolean.getBoolean("hipparchus.debug.sqp")) {
          opt.setDebugPrinter(System.out::println);
          }

        return opt.optimize(
            new InitialGuess(initialPoint()),
            new ObjectiveFunction(new HS391Obj())
        );
    }


    @Test
    public void testHS391() {
        

        LagrangeSolution sol = solve();
        double f = sol.getValue();
        final double expected = 0.0;
        assertEquals(expected, f, 1e-6);
    }
}
