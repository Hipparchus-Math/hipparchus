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
package org.hipparchus.analysis.integration;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class LebedevQuadratureTest {

    @Test
    void verifyAllOrders() {
        for (int order : LebedevQuadrature.AVAILABLE_ORDERS) {
            LebedevQuadrature.Rule rule = LebedevQuadrature.getRule(order);
            assertEquals(order, rule.getOrder(), "order mismatch for " + order);
            double wsum = 0;
            double maxNormErr = 0;
            for (int i = 0; i < order; i++) {
                double x = rule.getX()[i];
                double y = rule.getY()[i];
                double z = rule.getZ()[i];
                double norm = Math.sqrt(x * x + y * y + z * z);
                maxNormErr = Math.max(maxNormErr, Math.abs(norm - 1.0));
                wsum += rule.getW()[i];
            }
            System.out.println("order=" + order + " weightSum=" + wsum + " maxNormErr=" + maxNormErr);
            assertEquals(1.0, wsum, 1e-9, "weight sum mismatch for order " + order);
            assertTrue(maxNormErr < 1e-9, "points not on unit sphere for order " + order);
        }
    }
}
