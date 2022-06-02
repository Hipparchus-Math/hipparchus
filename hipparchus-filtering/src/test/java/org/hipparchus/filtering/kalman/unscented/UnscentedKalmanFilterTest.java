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

package org.hipparchus.filtering.kalman.unscented;

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
import org.hipparchus.util.MerweUnscentedTransform;
import org.junit.Assert;
import org.junit.Test;


public class UnscentedKalmanFilterTest {

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
        final UnscentedKalmanFilter<SimpleMeasurement> filter = new UnscentedKalmanFilter<>(new CholeskyDecomposer(1.0e-15, 1.0e-15),
                process, initial);

        // sequentially process all measurements and check against the reference estimated state and covariance
        measurements.
        map(measurement -> filter.estimationStep(measurement)).
        forEach(estimate -> {
            for (Reference r : referenceData) {
                if (r.sameTime(estimate.getTime())) {
                    r.checkState(estimate.getState(), 1.0e-4);
                    r.checkCovariance(estimate.getCovariance(), 3.0e-4);
                    return;
                }
            }
        });

    }

    private static class ConstantProcess implements UnscentedProcess<SimpleMeasurement> {

        private RealMatrix q = MatrixUtils.createRealDiagonalMatrix(new double[] {
            1.0e-5
        });
        @Override
        public UnscentedEvolution getEvolution(double previousTime, final RealVector[] previousStateSamples,  SimpleMeasurement measurement) {
            return new UnscentedEvolution(measurement.getTime(),
                                          previousStateSamples,
                                          previousStateSamples,
                                          q);
        }
        @Override
        public RealVector getInnovation(SimpleMeasurement measurement, RealVector predictedMeasurement, RealVector predictedState,
                RealMatrix innovationCovarianceMatrix) {
            return measurement.getValue().subtract(predictedMeasurement);
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
                         5.0e-4, 6.0e-4);
    }

    @Test
    public void testCannonballNonZeroProcessNoise() {
        doTestCannonball(new double[][] {
                            { 0.01, 0.00, 0.00, 0.00 },
                            { 0.00, 0.10, 0.00, 0.00 },
                            { 0.00, 0.00, 0.01, 0.00 },
                            { 0.00, 0.00, 0.00, 0.10 },
                         }, "cannonball-non-zero-process-noise.txt",
                         4.0e-4, 2.0e-2);
    }

    private void doTestCannonball(final double[][] q, final String name,
                                  final double tolState, final double tolCovariance) {
        final double mNoise   = 30.0;
        final double vIni     = 100.0;
        final double alphaIni = FastMath.PI / 4;
        final UnscentedProcess<SimpleMeasurement> process = new CannonballProcess(9.81, q);

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
        final UnscentedKalmanFilter<SimpleMeasurement> filter =
        new UnscentedKalmanFilter<>(new CholeskyDecomposer(1.0e-15, 1.0e-15), process, initial);

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

    private static class CannonballProcess implements UnscentedProcess<SimpleMeasurement> {
        private final double g;
        private final RealMatrix q;
        
        public CannonballProcess(final double g, final double[][] qData) {
            this.g = g;
            this.q = MatrixUtils.createRealMatrix(qData);
        }

        @Override
        public UnscentedEvolution getEvolution(double previousTime, RealVector[] previousStateSamples, SimpleMeasurement measurement) {
            final double dt = measurement.getTime() - previousTime;
            final RealVector[] states = new RealVector[9];
            final RealVector[] measurementSamples = new RealVector[9];
            for (int i = 0 ; i < 9 ; i++) {
                states[i]= MatrixUtils.createRealVector(new double[] {previousStateSamples[i].getEntry(0) + previousStateSamples[i].getEntry(1) * dt,
                        previousStateSamples[i].getEntry(1),
                        previousStateSamples[i].getEntry(2) + previousStateSamples[i].getEntry(3) * dt - 0.5 * g * dt * dt,
                        previousStateSamples[i].getEntry(3) - g * dt});
                measurementSamples[i]= MatrixUtils.createRealVector(new double[] { states[i].getEntry(0), states[i].getEntry(2) });

            }
            

            
            return new UnscentedEvolution(measurement.getTime(), states, measurementSamples, q);
        }

        @Override
        public RealVector getInnovation(SimpleMeasurement measurement, RealVector predictedMeasurement, RealVector predictedState,
                RealMatrix innovationCovarianceMatrix) {
            return measurement.getValue().subtract(predictedMeasurement);
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
        final UnscentedKalmanFilter<SimpleMeasurement> filter =
                        new UnscentedKalmanFilter<>(new CholeskyDecomposer(1.0e-15, 1.0e-15),
                                                   process, initial);

        // sequentially process all measurements and get only the last one
        ProcessEstimate finalEstimate = measurements.
                        map(measurement -> filter.estimationStep(measurement)).
                        reduce((first, second) -> second).get();

        Assert.assertEquals(expected, finalEstimate.getState().getEntry(0), tolerance);

    }

    private final class WelshBishopProcess implements UnscentedProcess<SimpleMeasurement> {

        private RealMatrix q;

        WelshBishopProcess(double qValue) {
            q = MatrixUtils.createRealDiagonalMatrix(new double[] {
                qValue
            });
        }
        @Override
        public UnscentedEvolution getEvolution(double previousTime,
                                               RealVector[] previousStates,
                                               SimpleMeasurement measurement) {
            return new UnscentedEvolution(measurement.getTime(),
                                          previousStates, previousStates
                                          ,q);
        }
        @Override
        public RealVector getInnovation(SimpleMeasurement measurement, RealVector predictedMeasurement, RealVector predictedState,
                RealMatrix innovationCovarianceMatrix) {
            return measurement.getValue().subtract(predictedMeasurement);
        }


    }
    
    @Test
    public void testRadar() {
        doTestRadar();
    }
    
    private void doTestRadar() {
        
        
        final ProcessEstimate initial = new ProcessEstimate(0.0, MatrixUtils.createRealVector(new double[] {0., 90., 1100.}),
                                                                 MatrixUtils.createRealMatrix(new double[][] {
                                                                 { 100., 0.00, 0.00},
                                                                 { 0.00, 100., 0.00},
                                                                 { 0.00, 0.00, 100.}}));
        final UnscentedProcess<SimpleMeasurement> process = new RadarProcess();
        MerweUnscentedTransform utProvider = new MerweUnscentedTransform(initial.getState().getDimension(), 0.001, 2, 0 );
        final UnscentedKalmanFilter<SimpleMeasurement> filter =
                new UnscentedKalmanFilter<>(new CholeskyDecomposer(1.0e-15, 1.0e-15), process, initial, utProvider);
        
        // Reference values generated from test_radar
        // available from https://github.com/rlabbe/filterpy/blob/master/filterpy/kalman/tests/test_ukf.py
        final List<Reference> referenceData = Reference.loadReferenceData(3, 1, "radar.txt");
        final Stream<SimpleMeasurement> measurements =
                referenceData.stream().
                map(r -> new SimpleMeasurement(r.getTime(),
                                               r.getZ(),
                                               MatrixUtils.createRealDiagonalMatrix(new double[] { 10.0 })));
        
        // sequentially process all measurements and check against the reference estimate
        measurements.
        map(measurement -> filter.estimationStep(measurement)).
        forEach(estimate -> {
            for (Reference r : referenceData) {
                if (r.sameTime(estimate.getTime())) {
                    r.checkState(estimate.getState(), 6.0e-6);
                    r.checkCovariance(estimate.getCovariance(), 5.0e-6);
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
    
    /**
     * Cross validation for the {@link UscentedKalmanFilter unscented kalman filter} with Roger Labbe's filterpy python library.
     * Class implementing test_radar from test_ukf. Data in "radar.txt" were generated using test_radar. 
     * @see "Roger Labbe's tests for Unscented Kalman Filter: https://github.com/rlabbe/filterpy/blob/master/filterpy/kalman/tests/test_ukf.py"
     */
    private static class RadarProcess implements UnscentedProcess<SimpleMeasurement> {

        public RadarProcess() {
        }
        

        @Override
        public UnscentedEvolution getEvolution(double previousTime, RealVector[] sigmaPoints, SimpleMeasurement measurement) {
            final double     dt    = measurement.getTime() - previousTime;

            final RealVector[] states = new RealVector[7];
            final RealVector[] measurementSamples = new RealVector[7];
            for (int i = 0; i < 7; i++) {
                
                states[i]= MatrixUtils.createRealVector(new double[] {
                        sigmaPoints[i].getEntry(0) + sigmaPoints[i].getEntry(1) * dt,
                        sigmaPoints[i].getEntry(1),
                        sigmaPoints[i].getEntry(2)
                    });
                measurementSamples[i]= MatrixUtils.createRealVector(new double[] { FastMath.sqrt(FastMath.pow(states[i].getEntry(0), 2) + FastMath.pow(states[i].getEntry(2), 2)) });
            }


            final RealMatrix processNoiseMatrix = MatrixUtils.createRealMatrix(new double[][] {
                { 0.01, 0.00, 0.00},
                { 0.00, 0.01, 0.00},
                { 0.00, 0.00, 0.01}
            });

            return new UnscentedEvolution(measurement.getTime(), states, measurementSamples,  processNoiseMatrix);
        }


        @Override
        public RealVector getInnovation(SimpleMeasurement measurement, RealVector predictedMeasurement, RealVector predictedState,
                RealMatrix innovationCovarianceMatrix) {
            return measurement.getValue().subtract(predictedMeasurement);
        }
        }




}
