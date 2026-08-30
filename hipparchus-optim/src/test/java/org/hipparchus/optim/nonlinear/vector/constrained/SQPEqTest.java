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
package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.linear.MatrixUtils;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

public class SQPEqTest {

    @Test
    public void testNoConstraints() {
        try {
            new SQPEq(new DummyProblem(false, false, false));
            Assertions.fail("an exception should have been thrown");
        } catch (NullArgumentException e) {
            // expected
        }
    }

    @Test
    public void testWithAllConstraints() {
        final SQPEq s = new SQPEq(new DummyProblem(true, true, true));
        Assertions.assertEquals(3, s.dim());
        Assertions.assertEquals(1, s.dimY());
        Assertions.assertArrayEquals(new double[] { 2 },
                                     s.value(MatrixUtils.createRealVector(new double[] { 1.0, 1.0, 1.0 })).toArray(),
                                     1.0e-12);
        Assertions.assertArrayEquals(new double[] { 2 },
                                     s.value(new double[] { 1.0, 1.0, 1.0 }),
                                     1.0e-12);
    }

}
