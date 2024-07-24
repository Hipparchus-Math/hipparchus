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

package org.hipparchus.optim;

import org.hipparchus.UnitTestUtils;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

class PointVectorValuePairTest {
    @Test
    void testSerial() {
        PointVectorValuePair pv1 = new PointVectorValuePair(new double[] { 1.0, 2.0, 3.0 },
                                                            new double[] { 4.0, 5.0 });
        PointVectorValuePair pv2 = (PointVectorValuePair) UnitTestUtils.serializeAndRecover(pv1);
        assertEquals(pv1.getKey().length, pv2.getKey().length);
        for (int i = 0; i < pv1.getKey().length; ++i) {
            assertEquals(pv1.getKey()[i], pv2.getKey()[i], 1.0e-15);
        }
        assertEquals(pv1.getValue().length, pv2.getValue().length);
        for (int i = 0; i < pv1.getValue().length; ++i) {
            assertEquals(pv1.getValue()[i], pv2.getValue()[i], 1.0e-15);
        }
    }
}
