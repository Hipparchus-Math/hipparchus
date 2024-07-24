/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */
package org.hipparchus.stat;

import org.hipparchus.stat.descriptive.DescriptiveStatistics;
import org.hipparchus.stat.descriptive.StreamingStatistics;
import org.junit.jupiter.api.Test;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;

/**
 * Certified data test cases.
 */
class CertifiedDataTest {

    protected double mean = Double.NaN;

    protected double std = Double.NaN;

    /**
     * Test SummaryStatistics - implementations that do not store the data
     * and use single pass algorithms to compute statistics
    */
    @Test
    void testSummaryStatistics() throws Exception {
        StreamingStatistics u = new StreamingStatistics();
        loadStats("data/PiDigits.txt", u);
        assertEquals(std, u.getStandardDeviation(), 1E-13, "PiDigits: std");
        assertEquals(mean, u.getMean(), 1E-13, "PiDigits: mean");

        loadStats("data/Mavro.txt", u);
        assertEquals(std, u.getStandardDeviation(), 1E-14, "Mavro: std");
        assertEquals(mean, u.getMean(), 1E-14, "Mavro: mean");

        loadStats("data/Michelso.txt", u);
        assertEquals(std, u.getStandardDeviation(), 1E-13, "Michelso: std");
        assertEquals(mean, u.getMean(), 1E-13, "Michelso: mean");

        loadStats("data/NumAcc1.txt", u);
        assertEquals(std, u.getStandardDeviation(), 1E-14, "NumAcc1: std");
        assertEquals(mean, u.getMean(), 1E-14, "NumAcc1: mean");

        loadStats("data/NumAcc2.txt", u);
        assertEquals(std, u.getStandardDeviation(), 1E-14, "NumAcc2: std");
        assertEquals(mean, u.getMean(), 1E-14, "NumAcc2: mean");
    }

    /**
     * Test DescriptiveStatistics - implementations that store full array of
     * values and execute multi-pass algorithms
     */
    @Test
    void testDescriptiveStatistics() throws Exception {

        DescriptiveStatistics u = new DescriptiveStatistics();

        loadStats("data/PiDigits.txt", u);
        assertEquals(std, u.getStandardDeviation(), 1E-14, "PiDigits: std");
        assertEquals(mean, u.getMean(), 1E-14, "PiDigits: mean");

        loadStats("data/Mavro.txt", u);
        assertEquals(std, u.getStandardDeviation(), 1E-14, "Mavro: std");
        assertEquals(mean, u.getMean(), 1E-14, "Mavro: mean");

        loadStats("data/Michelso.txt", u);
        assertEquals(std, u.getStandardDeviation(), 1E-14, "Michelso: std");
        assertEquals(mean, u.getMean(), 1E-14, "Michelso: mean");

        loadStats("data/NumAcc1.txt", u);
        assertEquals(std, u.getStandardDeviation(), 1E-14, "NumAcc1: std");
        assertEquals(mean, u.getMean(), 1E-14, "NumAcc1: mean");

        loadStats("data/NumAcc2.txt", u);
        assertEquals(std, u.getStandardDeviation(), 1E-14, "NumAcc2: std");
        assertEquals(mean, u.getMean(), 1E-14, "NumAcc2: mean");
    }

    /**
     * loads a DescriptiveStatistics off of a test file
     */
    private void loadStats(String resource, Object u) throws Exception {

        DescriptiveStatistics d = null;
        StreamingStatistics s = null;
        if (u instanceof DescriptiveStatistics) {
            d = (DescriptiveStatistics) u;
            d.clear();
        } else {
            s = (StreamingStatistics) u;
            s.clear();
        }

        mean = Double.NaN;
        std = Double.NaN;

        InputStream resourceAsStream = CertifiedDataTest.class.getResourceAsStream(resource);
        assertNotNull(resourceAsStream,"Could not find resource "+resource);
        BufferedReader in =
            new BufferedReader(
                    new InputStreamReader(
                            resourceAsStream));

        String line = null;

        for (int j = 0; j < 60; j++) {
            line = in.readLine();
            if (j == 40) {
                mean =
                    Double.parseDouble(
                            line.substring(line.lastIndexOf(":") + 1).trim());
            }
            if (j == 41) {
                std =
                    Double.parseDouble(
                            line.substring(line.lastIndexOf(":") + 1).trim());
            }
        }

        line = in.readLine();

        while (line != null) {
            if (d != null) {
                d.addValue(Double.parseDouble(line.trim()));
            }  else {
                s.addValue(Double.parseDouble(line.trim()));
            }
            line = in.readLine();
        }

        resourceAsStream.close();
        in.close();
    }
}
