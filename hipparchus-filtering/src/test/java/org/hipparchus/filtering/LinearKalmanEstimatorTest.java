/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.hipparchus.filtering;

import java.util.stream.IntStream;
import java.util.stream.Stream;

import org.hipparchus.linear.CholeskyDecomposer;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well1024a;
import org.junit.Assert;
import org.junit.Test;

public class LinearKalmanEstimatorTest {

    @Test
    public void testWelshBishopExactR() {
        doTestWelshBishop(0xd30a8f811e2f7c61l, -0.37727, 0.1,
                          0.0, 1.0, 1.0e-5, 0.1 * 0.1,
                          50, -0.389117, 1.0e-6);
    }

    @Test
    public void testWelshBishopBigR() {
        doTestWelshBishop(0xd30a8f811e2f7c61l, -0.37727, 0.1,
                          0.0, 1.0, 1.0e-5, 1.0 * 1.0,
                          50, -0.385613, 1.0e-6);
    }

    @Test
    public void testWelshBishopSmallR() {
        doTestWelshBishop(0xd30a8f811e2f7c61l, -0.37727, 0.1,
                          0.0, 1.0, 1.0e-5, 0.01 * 0.01,
                          50, -0.403015, 1.0e-6);
    }

    private void doTestWelshBishop(final long seed,
                                   final double trueConstant, final double trueStdv,
                                   final double initialEstimate, final double initialCovariance,
                                   final double q, final double r,
                                   final int nbMeasurements,
                                   final double expected, final double tolerance) {

        // this is the constant voltage example from paper
        // An Introduction to the Kalman Filter, Greg Welch and Gary Bishop
        // available from http://www.cs.unc.edu/~welch/media/pdf/kalman_intro.pdf
        final LinearProcess process = new ConstantMatricesProcess(MatrixUtils.createRealIdentityMatrix(1),
                                                            null, null,
                                                            MatrixUtils.createRealDiagonalMatrix(new double[] {
                                                                q
                                                            }));
        final ProcessEstimate initial = new ProcessEstimate(0,
                                                      MatrixUtils.createRealVector(new double[] { initialEstimate }),
                                                      MatrixUtils.createRealDiagonalMatrix(new double[] { initialCovariance }));
        final RandomGenerator generator = new Well1024a(seed);
        final Stream<Measurement> measurements =
                        IntStream.
                        range(0, nbMeasurements).
                        mapToObj(i -> new Measurement(i,
                                                      MatrixUtils.createRealVector(new double[] {
                                                          trueConstant + generator.nextGaussian() * trueStdv,
                                                      }),
                                                      MatrixUtils.createRealMatrix(new double[][] { { 1.0 } }),
                                                      MatrixUtils.createRealDiagonalMatrix(new double[] { r })));

        // set up Kalan estimator
        final LinearKalmanEstimator estimator =
                        new LinearKalmanEstimator(new CholeskyDecomposer(1.0e-15, 1.0e-15), process, initial);

        // sequentially process all measurements and get only the lat one
        final Stream<ProcessEstimate> estimates = estimator.estimate(measurements);
        final ProcessEstimate finalEstimate = estimates.reduce((first, second) -> second).get();

        Assert.assertEquals(expected, finalEstimate.getState().getEntry(0), tolerance);

    }

    private static class ConstantMatricesProcess implements LinearProcess {

        private final RealMatrix a;
        private final RealMatrix b;
        private final RealVector u;
        private final RealMatrix q;

        
        public ConstantMatricesProcess(final RealMatrix a,
                                       final RealMatrix b,
                                       final RealVector u,
                                       final RealMatrix q) {
            this.a = a;
            this.b = b;
            this.u = u;
            this.q = q;
        }

        @Override
        public RealMatrix getStateTransitionMatrix(double time) {
            return a;
        }

        @Override
        public RealMatrix getControlMatrix(double time) {
            return b;
        }

        @Override
        public RealVector getCommand(double time) {
            return u;
        }

        @Override
        public RealMatrix getProcessNoiseMatrix(double time) {
            return q;
        }

    };
}
