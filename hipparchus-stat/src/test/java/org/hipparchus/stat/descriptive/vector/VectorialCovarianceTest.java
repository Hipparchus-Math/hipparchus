//Licensed to the Apache Software Foundation (ASF) under one
//or more contributor license agreements.  See the NOTICE file
//distributed with this work for additional information
//regarding copyright ownership.  The ASF licenses this file
//to you under the Apache License, Version 2.0 (the
//"License"); you may not use this file except in compliance
//with the License.  You may obtain a copy of the License at

//https://www.apache.org/licenses/LICENSE-2.0

//Unless required by applicable law or agreed to in writing,
//software distributed under the License is distributed on an
//"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
//KIND, either express or implied.  See the License for the
//specific language governing permissions and limitations
//under the License.

package org.hipparchus.stat.descriptive.vector;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.linear.RealMatrix;
import org.junit.Test;

public class VectorialCovarianceTest {
    private double[][] points;

    public VectorialCovarianceTest() {
        points = new double[][] {
            { 1.2, 2.3,  4.5},
            {-0.7, 2.3,  5.0},
            { 3.1, 0.0, -3.1},
            { 6.0, 1.2,  4.2},
            {-0.7, 2.3,  5.0}
        };
    }

    @Test
    public void testMismatch() {
        try {
            new VectorialCovariance(8, true).increment(new double[5]);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException dme) {
            assertEquals(5, ((Integer) dme.getParts()[0]).intValue());
            assertEquals(8, ((Integer) dme.getParts()[1]).intValue());
        }
    }

    @Test
    public void testSimplistic() {
        VectorialCovariance stat = new VectorialCovariance(2, true);
        stat.increment(new double[] {-1.0,  1.0});
        stat.increment(new double[] { 1.0, -1.0});
        RealMatrix c = stat.getResult();
        assertEquals( 2.0, c.getEntry(0, 0), 1.0e-12);
        assertEquals(-2.0, c.getEntry(1, 0), 1.0e-12);
        assertEquals( 2.0, c.getEntry(1, 1), 1.0e-12);
    }

    @Test
    public void testBasicStats() {

        VectorialCovariance stat = new VectorialCovariance(points[0].length, true);
        for (int i = 0; i < points.length; ++i) {
            stat.increment(points[i]);
        }

        assertEquals(points.length, stat.getN());

        RealMatrix c = stat.getResult();
        double[][] refC    = new double[][] {
                { 8.0470, -1.9195, -3.4445},
                {-1.9195,  1.0470,  3.2795},
                {-3.4445,  3.2795, 12.2070}
        };

        for (int i = 0; i < c.getRowDimension(); ++i) {
            for (int j = 0; j <= i; ++j) {
                assertEquals(refC[i][j], c.getEntry(i, j), 1.0e-12);
            }
        }

    }

    @Test
    public void testSerial(){
        VectorialCovariance stat = new VectorialCovariance(points[0].length, true);
        assertEquals(stat, UnitTestUtils.serializeAndRecover(stat));
    }
}
