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

package org.hipparchus.filtering.kalman;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.util.FastMath;
import org.junit.Assert;

public class Reference {

    private final double     time;
    private final RealVector z;
    private final RealVector s;
    private final RealMatrix c;
    private final RealMatrix stm;
    private final RealMatrix h;
    private final RealMatrix icm;
    private final RealMatrix k;

    private Reference(final int stateDimension, final int measurementDimension, final String line) {
        final String[] fields = line.split("\\s+");
        int index = 0;
        time = Double.parseDouble(fields[index++]);
        z = MatrixUtils.createRealVector(new double[measurementDimension]);
        for (int i = 0; i < measurementDimension; ++i) {
            z.setEntry(i, Double.parseDouble(fields[index++]));
        }
        s = MatrixUtils.createRealVector(new double[stateDimension]);
        for (int i = 0; i < stateDimension; ++i) {
            s.setEntry(i, Double.parseDouble(fields[index++]));
        }
        c = MatrixUtils.createRealMatrix(stateDimension, stateDimension);
        for (int i = 0; i < stateDimension; ++i) {
            for (int j = i; j < stateDimension; ++j) {
                c.setEntry(i, j, Double.parseDouble(fields[index++]));
                c.setEntry(j, i, c.getEntry(i, j));
            }
        }
        if (fields.length > index) {
            // there are additional data

            stm = MatrixUtils.createRealMatrix(stateDimension, stateDimension);
            for (int i = 0; i < stateDimension; ++i) {
                for (int j = 0; j < stateDimension; ++j) {
                    stm.setEntry(i, j, Double.parseDouble(fields[index++]));
                }
            }

            h = MatrixUtils.createRealMatrix(measurementDimension, stateDimension);
            for (int i = 0; i < measurementDimension; ++i) {
                for (int j = 0; j < stateDimension; ++j) {
                    h.setEntry(i, j, Double.parseDouble(fields[index++]));
                }
            }

            icm = MatrixUtils.createRealMatrix(measurementDimension, measurementDimension);
            for (int i = 0; i < measurementDimension; ++i) {
                for (int j = i; j < measurementDimension; ++j) {
                    icm.setEntry(i, j, Double.parseDouble(fields[index++]));
                    icm.setEntry(j, i, icm.getEntry(i, j));
                }
            }

            k = MatrixUtils.createRealMatrix(stateDimension, measurementDimension);
            for (int i = 0; i < stateDimension; ++i) {
                for (int j = 0; j < measurementDimension; ++j) {
                    k.setEntry(i, j, Double.parseDouble(fields[index++]));
                }
            }

        } else {
            stm = null;
            h   = null;
            icm = null;
            k   = null;
        }
    }

    public static List<Reference> loadReferenceData(final int stateDimension, final int measurementDimension,
                                                    final String name) {
        List<Reference> loaded = new ArrayList<>();
        try (InputStream is = KalmanFilter.class.getResourceAsStream(name);
             InputStreamReader isr = new InputStreamReader(is, "UTF-8");
             BufferedReader br = new BufferedReader(isr)) {
            for (String line = br.readLine(); line != null; line = br.readLine()) {
                line = line.trim();
                if (line.length() > 0 && !line.startsWith("#")) {
                    loaded.add(new Reference(stateDimension, measurementDimension, line));
                }
            }
        } catch (IOException ioe) {
            Assert.fail(ioe.getLocalizedMessage());
        }
        return loaded;
    }

    public boolean sameTime(final double otherTime) {
        return FastMath.abs(time - otherTime) < 1.0e-6;
    }

    public void checkState(final RealVector otherState, final double tolerance) {
        checkVector(s, otherState, tolerance);
    }

    public void checkCovariance(final RealMatrix otherCovariance, final double tolerance) {
        checkMatrix(c, otherCovariance, tolerance);
    }

    public boolean hasIntermediateData() {
        return stm != null;
    }

    public void checkStateTransitionMatrix(final RealMatrix otherSTM, final double tolerance) {
        checkMatrix(stm, otherSTM, tolerance);
    }

    public void checkMeasurementJacobian(final RealMatrix otherMeasurementJacobian, final double tolerance) {
        checkMatrix(h, otherMeasurementJacobian, tolerance);
    }

    public void checkInnovationCovariance(final RealMatrix otherInnovationCovariance, final double tolerance) {
        checkMatrix(icm, otherInnovationCovariance, tolerance);
    }

    public void checkKalmanGain(final RealMatrix otherKalmanGain, final double tolerance) {
        checkMatrix(k, otherKalmanGain, tolerance);
    }

    public double getTime() {
        return time;
    }

    public RealVector getZ() {
        return z;
    }

    private void checkVector(final RealVector referenceVector, final RealVector otherVector, final double tolerance) {
        Assert.assertEquals(referenceVector.getDimension(), otherVector.getDimension());
        for (int i = 0; i < referenceVector.getDimension(); ++i) {
            Assert.assertEquals(time + ": ", referenceVector.getEntry(i), otherVector.getEntry(i), tolerance);
        }
    }

    private void checkMatrix(final RealMatrix referenceMatrix, final RealMatrix otherMatrix, final double tolerance) {
        Assert.assertEquals(referenceMatrix.getRowDimension(), otherMatrix.getRowDimension());
        Assert.assertEquals(referenceMatrix.getColumnDimension(), otherMatrix.getColumnDimension());
        for (int i = 0; i < referenceMatrix.getRowDimension(); ++i) {
            for (int j = i; j < referenceMatrix.getColumnDimension(); ++j) {
                Assert.assertEquals(time + ": ", referenceMatrix.getEntry(i, j), otherMatrix.getEntry(i, j), tolerance);
            }
        }
    }

}
