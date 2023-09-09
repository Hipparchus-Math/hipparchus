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

package org.hipparchus.filtering.kalman.extended;

import java.util.List;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import org.hipparchus.filtering.kalman.ProcessEstimate;
import org.hipparchus.filtering.kalman.Reference;
import org.hipparchus.filtering.kalman.SimpleMeasurement;
import org.hipparchus.linear.CholeskyDecomposer;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well1024a;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

public class ExtendedKalmanFilterTest {

    @Test
    public void testConstant() {

        ConstantProcess process = new ConstantProcess();

        // initial estimate is perfect, and process noise is perfectly known
        final ProcessEstimate initial = new ProcessEstimate(0,
                                                            MatrixUtils.createRealVector(new double[] { 10.0 }),
                                                            process.q);
        Assert.assertNull(initial.getInnovationCovariance());

        // reference values from Apache Commons Math 3.6.1 unit test
        final List<Reference> referenceData = Reference.loadReferenceData(1, 1, "constant-value.txt");
        final Stream<SimpleMeasurement> measurements =
                        referenceData.stream().
                        map(r -> new SimpleMeasurement(r.getTime(),
                                                       r.getZ(),
                                                       MatrixUtils.createRealDiagonalMatrix(new double[] { 0.1 })));

        // set up Kalman filter
        final ExtendedKalmanFilter<SimpleMeasurement> filter =
                        new ExtendedKalmanFilter<>(new CholeskyDecomposer(1.0e-15, 1.0e-15),
                                                   process, initial);

        // sequentially process all measurements and check against the reference estimated state and covariance
        measurements.
        map(measurement -> filter.estimationStep(measurement)).
        forEach(estimate -> {
            for (Reference r : referenceData) {
                if (r.sameTime(estimate.getTime())) {
                    r.checkState(estimate.getState(), 1.0e-15);
                    r.checkCovariance(estimate.getCovariance(), 3.0e-19);
                    return;
                }
            }
        });

    }

    private static class ConstantProcess implements NonLinearProcess<SimpleMeasurement> {

        private RealMatrix q = MatrixUtils.createRealDiagonalMatrix(new double[] {
            1.0e-5
        });

        @Override
        public NonLinearEvolution getEvolution(double previousTime, RealVector previousState, SimpleMeasurement measurement) {
            return new NonLinearEvolution(measurement.getTime(),
                                          previousState,
                                          MatrixUtils.createRealIdentityMatrix(1),
                                          q,
                                          MatrixUtils.createRealMatrix(new double[][] { { 1.0 } }));
        }

        @Override
        public RealVector getInnovation(SimpleMeasurement measurement, NonLinearEvolution evolution, RealMatrix innovationCovarianceMatrix) {
            return measurement.getValue().subtract(evolution.getCurrentState());
        }

    }

    @Test
    public void testConstantAcceleration() {
        doTestConstantAcceleration("constant-acceleration.txt");
    }

    @Test
    public void testConstantAccelerationWithIntermediateData() {
        doTestConstantAcceleration("constant-acceleration-with-intermediate-data.txt");
    }

    @Test
    public void testConstantAccelerationWithOutlier() {
        doTestConstantAcceleration("constant-acceleration-with-outlier.txt");
    }

    private void doTestConstantAcceleration(String name) {

        final double acc    = 0.1;
        final double aNoise = 0.2;
        final double mNoise = 10.0;
        final NonLinearProcess<SimpleMeasurement> process = new ConstantAccelerationProcess(acc, aNoise);

        // initial state is estimated to be at rest on origin
        final ProcessEstimate initial = new ProcessEstimate(0,
                                                            MatrixUtils.createRealVector(new double[] { 0.0, 0.0 }),
                                                            MatrixUtils.createRealMatrix(new double[][] {
                                                                { 1.0, 1.0 },
                                                                { 1.0, 1.0 }
                                                            }));

        // reference values from Apache Commons Math 3.6.1 unit test
        final List<Reference> referenceData = Reference.loadReferenceData(2, 1, name);
        final Stream<SimpleMeasurement> measurements =
                        referenceData.stream().
                        map(r -> new SimpleMeasurement(r.getTime(),
                                                       r.getZ(),
                                                       MatrixUtils.createRealDiagonalMatrix(new double[] { mNoise * mNoise })));

        // set up Kalman filter
        final ExtendedKalmanFilter<SimpleMeasurement> filter =
        new ExtendedKalmanFilter<>(new CholeskyDecomposer(1.0e-15, 1.0e-15), process, initial);

        // sequentially process all measurements and check against the reference estimate
        measurements.
        map(measurement -> filter.estimationStep(measurement)).
        forEach(estimate -> {
            for (Reference r : referenceData) {
                if (r.sameTime(estimate.getTime())) {
                    r.checkState(estimate.getState(), 6.0e-15);
                    r.checkCovariance(estimate.getCovariance(), 5.0e-15);
                    if (r.hasIntermediateData()) {
                        r.checkStateTransitionMatrix(estimate.getStateTransitionMatrix(), 1.0e-14);
                        r.checkMeasurementJacobian(estimate.getMeasurementJacobian(),     1.0e-15);
                        r.checkInnovationCovariance(estimate.getInnovationCovariance(),   1.0e-12);
                        r.checkKalmanGain(estimate.getKalmanGain(),                       1.0e-12);
                        r.checkKalmanGain(estimate.getKalmanGain(),                       1.0e-15);
                    }
                    return;
                }
            }
        });

    }

    private static class ConstantAccelerationProcess implements NonLinearProcess<SimpleMeasurement> {
        private final double acc;
        private final double aNoise2;
        
        public ConstantAccelerationProcess(final double acc, final double aNoise) {
            this.acc     = acc;
            this.aNoise2 = aNoise * aNoise;
        }

        @Override
        public NonLinearEvolution getEvolution(double previousTime, RealVector previousState, SimpleMeasurement measurement) {
            final double     dt    = measurement.getTime() - previousTime;
            final double     dt2   = dt  * dt;
            final double     dt3   = dt2 * dt;
            final double     dt4   = dt2 * dt2;
            final RealVector state = MatrixUtils.createRealVector(new double[] {
                previousState.getEntry(0) + previousState.getEntry(1) * dt + 0.5 * acc * dt * dt,
                previousState.getEntry(1) + acc * dt
            });
            final RealMatrix stm = MatrixUtils.createRealMatrix(new double[][] {
                { 1.0,  dt },
                { 0.0, 1.0 }
            });
            final RealMatrix processNoiseMatrix = MatrixUtils.createRealMatrix(new double[][] {
                { 0.25 * dt4 * aNoise2, 0.5 * dt3 * aNoise2 },
                { 0.5  * dt3 * aNoise2, dt2 * aNoise2 }
            });
            RealMatrix h = (measurement.getValue().getEntry(0) > 1.0e6) ?
                           null : MatrixUtils.createRealMatrix(new double[][] { { 1.0, 0.0 } });
            return new NonLinearEvolution(measurement.getTime(), state, stm, processNoiseMatrix, h);
        }

        @Override
        public RealVector getInnovation(SimpleMeasurement measurement,
                                        NonLinearEvolution evolution,
                                        RealMatrix innovationCovarianceMatrix) {
            return measurement.getValue().subtract(evolution.getCurrentState().getSubVector(0, 1));
        }

    }

    @Test
    public void testCannonballZeroProcessNoise() {
        doTestCannonball(new double[][] {
                            { 0.00, 0.00, 0.00, 0.00 },
                            { 0.00, 0.00, 0.00, 0.00 },
                            { 0.00, 0.00, 0.00, 0.00 },
                            { 0.00, 0.00, 0.00, 0.00 },
                         }, "cannonball-zero-process-noise.txt",
                         5.0e-13, 6.0e-14);
    }

    @Test
    public void testCannonballNonZeroProcessNoise() {
        doTestCannonball(new double[][] {
                            { 0.01, 0.00, 0.00, 0.00 },
                            { 0.00, 0.10, 0.00, 0.00 },
                            { 0.00, 0.00, 0.01, 0.00 },
                            { 0.00, 0.00, 0.00, 0.10 },
                         }, "cannonball-non-zero-process-noise.txt",
                         4.0e-13, 2.0e-13);
    }

    private void doTestCannonball(final double[][] q, final String name,
                                  final double tolState, final double tolCovariance) {

        final double mNoise   = 30.0;
        final double vIni     = 100.0;
        final double alphaIni = FastMath.PI / 4;
        final NonLinearProcess<SimpleMeasurement> process = new CannonballProcess(9.81, q);

        // initial state is estimated to be a shot from origin with known angle and velocity
        final ProcessEstimate initial = new ProcessEstimate(0.0,
                                                            MatrixUtils.createRealVector(new double[] {
                                                                0.0, vIni * FastMath.cos(alphaIni),
                                                                0.0, vIni * FastMath.sin(alphaIni)
                                                            }),
                                                            MatrixUtils.createRealDiagonalMatrix(new double[] {
                                                                mNoise * mNoise, 1.0e-3, mNoise * mNoise, 1.0e-3
                                                            }));

        // reference values from Apache Commons Math 3.6.1 unit test
        // we have changed the test slightly, setting up a non-zero process noise
        final List<Reference> referenceData = Reference.loadReferenceData(4, 2, name);
        final Stream<SimpleMeasurement> measurements =
                        referenceData.stream().
                        map(r -> new SimpleMeasurement(r.getTime(),
                                                       r.getZ(),
                                                       MatrixUtils.createRealDiagonalMatrix(new double[] {
                                                           mNoise * mNoise, mNoise * mNoise
                                                       })));

        // set up Kalman filter
        final ExtendedKalmanFilter<SimpleMeasurement> filter =
        new ExtendedKalmanFilter<>(new CholeskyDecomposer(1.0e-15, 1.0e-15), process, initial);

        // sequentially process all measurements and check against the reference estimate
        measurements.
        map(measurement -> filter.estimationStep(measurement)).
        map(estimate -> {
            final ProcessEstimate p = filter.getPredicted();
            final ProcessEstimate c = filter.getCorrected();
            Assert.assertEquals(p.getTime(), c.getTime(), 1.0e-15);
            Assert.assertTrue(p.getState().getDistance(c.getState()) > 0.005);
            return estimate;
        }).
        forEach(estimate -> {
            for (Reference r : referenceData) {
                if (r.sameTime(estimate.getTime())) {
                    r.checkState(estimate.getState(), tolState);
                    r.checkCovariance(estimate.getCovariance(), tolCovariance);
                    return;
                }
            }
        });

    }

    private static class CannonballProcess implements NonLinearProcess<SimpleMeasurement> {
        private final double g;
        private final RealMatrix q;
        
        public CannonballProcess(final double g, final double[][] qData) {
            this.g = g;
            this.q = MatrixUtils.createRealMatrix(qData);
        }

        @Override
        public NonLinearEvolution getEvolution(double previousTime, RealVector previousState, SimpleMeasurement measurement) {
            final double dt = measurement.getTime() - previousTime;
            final RealVector state = MatrixUtils.createRealVector(new double[] {
                previousState.getEntry(0) + previousState.getEntry(1) * dt,
                previousState.getEntry(1),
                previousState.getEntry(2) + previousState.getEntry(3) * dt - 0.5 * g * dt * dt,
                previousState.getEntry(3) - g * dt
            });
            final RealMatrix stm = MatrixUtils.createRealMatrix(new double[][] {
                { 1.0,  dt, 0.0, 0.0 },
                { 0.0, 1.0, 0.0, 0.0 },
                { 0.0, 0.0, 1.0,  dt },
                { 0.0, 0.0, 0.0, 1.0 },
            });
            return new NonLinearEvolution(measurement.getTime(), state, stm, q,
                                          MatrixUtils.createRealMatrix(new double[][] {
                                              { 1.0, 0.0, 0.0, 0.0 },
                                              { 0.0, 0.0, 1.0, 0.0 }
                                          }));
        }

        @Override
        public RealVector getInnovation(SimpleMeasurement measurement, NonLinearEvolution evolution,
                                        RealMatrix innovationCovarianceMatrix) {
            return measurement.getValue().
                            subtract(MatrixUtils.createRealVector(new double[] {
                                evolution.getCurrentState().getEntry(0),
                                evolution.getCurrentState().getEntry(2)
                            }));
        }

    }

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

        WelshBishopProcess process = new WelshBishopProcess(q);

        // this is the constant voltage example from paper
        // An Introduction to the Kalman Filter, Greg Welch and Gary Bishop
        // available from http://www.cs.unc.edu/~welch/media/pdf/kalman_intro.pdf
        final ProcessEstimate initial = new ProcessEstimate(0,
                                                            MatrixUtils.createRealVector(new double[] { initialEstimate }),
                                                            MatrixUtils.createRealDiagonalMatrix(new double[] { initialCovariance }));
        final RandomGenerator generator = new Well1024a(seed);
        final Stream<SimpleMeasurement> measurements =
                        IntStream.
                        range(0, nbMeasurements).
                        mapToObj(i -> new SimpleMeasurement(i,
                                                            MatrixUtils.createRealVector(new double[] {
                                                                trueConstant + generator.nextGaussian() * trueStdv,
                                                            }),
                                                            MatrixUtils.createRealDiagonalMatrix(new double[] { r })));

        // set up Kalman filter
        final ExtendedKalmanFilter<SimpleMeasurement> filter =
                        new ExtendedKalmanFilter<>(new CholeskyDecomposer(1.0e-15, 1.0e-15),
                                                   process, initial);

        // sequentially process all measurements and get only the last one
        ProcessEstimate finalEstimate = measurements.
                        map(measurement -> filter.estimationStep(measurement)).
                        reduce((first, second) -> second).get();

        Assert.assertEquals(expected, finalEstimate.getState().getEntry(0), tolerance);

    }

    private final class WelshBishopProcess implements NonLinearProcess<SimpleMeasurement> {

        private RealMatrix q;

        WelshBishopProcess(double qValue) {
            q = MatrixUtils.createRealDiagonalMatrix(new double[] {
                qValue
            });
        }
        @Override
        public NonLinearEvolution getEvolution(double previousTime,
                                               RealVector previousState,
                                               SimpleMeasurement measurement) {
            return new NonLinearEvolution(measurement.getTime(),
                                          previousState,
                                          MatrixUtils.createRealIdentityMatrix(1),
                                          q,
                                          MatrixUtils.createRealMatrix(new double[][] { { 1.0 } }));
        }

        @Override
        public RealVector getInnovation(SimpleMeasurement measurement,
                                        NonLinearEvolution evolution,
                                        RealMatrix innovationCovarianceMatrix) {
            return measurement.getValue().subtract(evolution.getCurrentState());
        }

    }

}
