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
import org.hipparchus.stat.descriptive.moment.Mean;
import org.junit.Test;

public class VectorialStorelessStatisticTest {
    private double[][] points;

    public VectorialStorelessStatisticTest() {
        points = new double[][] {
            { 1.2, 2.3,  4.5},
            {-0.7, 2.3,  5.0},
            { 3.1, 0.0, -3.1},
            { 6.0, 1.2,  4.2},
            {-0.7, 2.3,  5.0}
        };
    }

    protected VectorialStorelessStatistic createStatistic(int dimension) {
        return new VectorialStorelessStatistic(dimension, new Mean());
    }

    @Test
    public void testMismatch() {
        try {
            createStatistic(8).increment(new double[5]);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException dme) {
            assertEquals(5, ((Integer) dme.getParts()[0]).intValue());
            assertEquals(8, ((Integer) dme.getParts()[1]).intValue());
        }
    }

    @Test
    public void testSimplistic() {
        VectorialStorelessStatistic stat = createStatistic(2);
        stat.increment(new double[] {-1.0,  1.0});
        stat.increment(new double[] { 1.0, -1.0});
        double[] mean = stat.getResult();
        assertEquals(0.0, mean[0], 1.0e-12);
        assertEquals(0.0, mean[1], 1.0e-12);
    }

    @Test
    public void testBasicStats() {

        VectorialStorelessStatistic stat = createStatistic(points[0].length);
        for (int i = 0; i < points.length; ++i) {
            stat.increment(points[i]);
        }

        assertEquals(points.length, stat.getN());

        double[] mean = stat.getResult();
        double[]   refMean = new double[] { 1.78, 1.62,  3.12};

        for (int i = 0; i < mean.length; ++i) {
            assertEquals(refMean[i], mean[i], 1.0e-12);
        }
    }

    @Test
    public void testSerial() {
        VectorialStorelessStatistic stat = createStatistic(points[0].length);
        for (int i = 0; i < points.length; ++i) {
            stat.increment(points[i]);
        }
        assertEquals(stat, UnitTestUtils.serializeAndRecover(stat));
    }
}
