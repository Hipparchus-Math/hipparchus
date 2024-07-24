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
package org.hipparchus.stat.descriptive.moment;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;


/**
 * Test class for use of external moments.
 */
class InteractionTest {

    protected double mean = 12.40454545454550;
    protected double var = 10.00235930735930;
    protected double skew = 1.437423729196190;
    protected double kurt = 2.377191264804700;

    protected double tolerance = 10E-12;

    protected double[] testArray = {
        12.5, 12, 11.8, 14.2, 14.9, 14.5, 21, 8.2, 10.3, 11.3, 14.1,
        9.9, 12.2, 12, 12.1, 11, 19.8, 11, 10, 8.8, 9, 12.3
    };

    @Test
    void testInteraction() {

        FourthMoment m4 = new FourthMoment();
        Mean m = new Mean(m4);
        Variance v = new Variance(m4);
        Skewness s= new Skewness(m4);
        Kurtosis k = new Kurtosis(m4);

        for (int i = 0; i < testArray.length; i++){
            m4.increment(testArray[i]);
        }

        assertEquals(mean, m.getResult(), tolerance);
        assertEquals(var, v.getResult(), tolerance);
        assertEquals(skew, s.getResult(), tolerance);
        assertEquals(kurt, k.getResult(), tolerance);
    }

}
