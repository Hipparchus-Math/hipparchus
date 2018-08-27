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
        Assert.assertEquals(s.getDimension(), otherState.getDimension());
        for (int i = 0; i < s.getDimension(); ++i) {
            Assert.assertEquals(time + ": ", s.getEntry(i), otherState.getEntry(i), tolerance);
        }
    }

    public void checkcovariance(final RealMatrix otherCovariance, final double tolerance) {
        Assert.assertEquals(c.getRowDimension(), otherCovariance.getRowDimension());
        Assert.assertEquals(c.getColumnDimension(), otherCovariance.getColumnDimension());
        for (int i = 0; i < c.getRowDimension(); ++i) {
            for (int j = i; j < c.getColumnDimension(); ++j) {
                Assert.assertEquals(time + ": ", c.getEntry(i, j), otherCovariance.getEntry(i, j), tolerance);
            }
        }
    }

    public double getTime() {
        return time;
    }

    public RealVector getZ() {
        return z;
    }

}
