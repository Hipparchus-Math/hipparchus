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
package org.hipparchus.linear;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertThrows;

public final class OpenMapRealMatrixTest {

    @Test
    public void testMath679() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            new OpenMapRealMatrix(3, Integer.MAX_VALUE);
        });
    }

    @Test
    public void testMath870() {
        // Caveat: This implementation assumes that, for any {@code x},
        // the equality {@code x * 0d == 0d} holds. But it is is not true for
        // {@code NaN}. Moreover, zero entries will lose their sign.
        // Some operations (that involve {@code NaN} and/or infinities) may
        // thus give incorrect results.
        OpenMapRealMatrix a = new OpenMapRealMatrix(3, 3);
        OpenMapRealMatrix x = new OpenMapRealMatrix(3, 1);
        x.setEntry(0, 0, Double.NaN);
        x.setEntry(2, 0, Double.NEGATIVE_INFINITY);
        OpenMapRealMatrix b = a.multiply(x);
        for (int i = 0; i < b.getRowDimension(); ++i) {
            for (int j = 0; j < b.getColumnDimension(); ++j) {
                // NaNs and infinities have disappeared, this is a limitation of our implementation
                Assertions.assertEquals(0.0, b.getEntry(i, j), 1.0e-20);
            }
        }
    }
}
