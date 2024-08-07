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
package org.hipparchus.stat.inference;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.stat.descriptive.StreamingStatistics;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;


/**
 * Test cases for the OneWayAnovaImpl class.
 *
 */

class OneWayAnovaTest {

    protected OneWayAnova testStatistic = new OneWayAnova();

    private double[] emptyArray = {};

    private double[] classA =
            {93.0, 103.0, 95.0, 101.0, 91.0, 105.0, 96.0, 94.0, 101.0 };
    private double[] classB =
            {99.0, 92.0, 102.0, 100.0, 102.0, 89.0 };
    private double[] classC =
            {110.0, 115.0, 111.0, 117.0, 128.0, 117.0 };

    @Test
    void testAnovaFValue() {
        // Target comparison values computed using R version 2.6.0 (Linux version)
        List<double[]> threeClasses = new ArrayList<double[]>();
        threeClasses.add(classA);
        threeClasses.add(classB);
        threeClasses.add(classC);

        assertEquals(24.67361709460624,
                 testStatistic.anovaFValue(threeClasses), 1E-12, "ANOVA F-value");

        List<double[]> twoClasses = new ArrayList<double[]>();
        twoClasses.add(classA);
        twoClasses.add(classB);

        assertEquals(0.0150579150579,
                 testStatistic.anovaFValue(twoClasses), 1E-12, "ANOVA F-value");

        List<double[]> emptyContents = new ArrayList<double[]>();
        emptyContents.add(emptyArray);
        emptyContents.add(classC);
        try {
            testStatistic.anovaFValue(emptyContents);
            fail("empty array for key classX, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        List<double[]> tooFew = new ArrayList<double[]>();
        tooFew.add(classA);
        try {
            testStatistic.anovaFValue(tooFew);
            fail("less than two classes, MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }


    @Test
    void testAnovaPValue() {
        // Target comparison values computed using R version 2.6.0 (Linux version)
        List<double[]> threeClasses = new ArrayList<double[]>();
        threeClasses.add(classA);
        threeClasses.add(classB);
        threeClasses.add(classC);

        assertEquals(6.959446E-06,
                 testStatistic.anovaPValue(threeClasses), 1E-12, "ANOVA P-value");

        List<double[]> twoClasses = new ArrayList<double[]>();
        twoClasses.add(classA);
        twoClasses.add(classB);

        assertEquals(0.904212960464,
                 testStatistic.anovaPValue(twoClasses), 1E-12, "ANOVA P-value");

    }

    @Test
    void testAnovaPValueSummaryStatistics() {
        // Target comparison values computed using R version 2.6.0 (Linux version)
        List<StreamingStatistics> threeClasses = new ArrayList<StreamingStatistics>();
        StreamingStatistics statsA = new StreamingStatistics();
        for (double a : classA) {
            statsA.addValue(a);
        }
        threeClasses.add(statsA);
        StreamingStatistics statsB = new StreamingStatistics();
        for (double b : classB) {
            statsB.addValue(b);
        }
        threeClasses.add(statsB);
        StreamingStatistics statsC = new StreamingStatistics();
        for (double c : classC) {
            statsC.addValue(c);
        }
        threeClasses.add(statsC);

        assertEquals(6.959446E-06,
                 testStatistic.anovaPValue(threeClasses, true), 1E-12, "ANOVA P-value");

        List<StreamingStatistics> twoClasses = new ArrayList<StreamingStatistics>();
        twoClasses.add(statsA);
        twoClasses.add(statsB);

        assertEquals(0.904212960464,
                 testStatistic.anovaPValue(twoClasses, false), 1E-12, "ANOVA P-value");

    }

    @Test
    void testAnovaTest() {
        // Target comparison values computed using R version 2.3.1 (Linux version)
        List<double[]> threeClasses = new ArrayList<double[]>();
        threeClasses.add(classA);
        threeClasses.add(classB);
        threeClasses.add(classC);

        assertTrue(testStatistic.anovaTest(threeClasses, 0.01), "ANOVA Test P<0.01");

        List<double[]> twoClasses = new ArrayList<double[]>();
        twoClasses.add(classA);
        twoClasses.add(classB);

        assertFalse(testStatistic.anovaTest(twoClasses, 0.01), "ANOVA Test P>0.01");
    }

}
