/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.linear;

import org.junit.Assert;
import org.junit.Test;

public class RealLinearOperatorTest {

    @Test
    public void testDefaultIsTransposable() {
        Assert.assertFalse(new DefaultOperator(MatrixUtils.createRealIdentityMatrix(2)).isTransposable());
    }

    @Test
    public void testDefaultTransposeMultiply() {
        try {
            new DefaultOperator(MatrixUtils.createRealIdentityMatrix(2)).operateTranspose(MatrixUtils.createRealVector(2));
            Assert.fail("an exception should have been thrown");
        } catch (UnsupportedOperationException uoe) {
            // expected
        }
    }

    // local class that does NOT override isTransposable nor operateTranspose
    // so the default methods are called
    private class DefaultOperator implements RealLinearOperator {

        RealMatrix m;
        public DefaultOperator(RealMatrix m) {
            this.m = m;
        }

        @Override
        public int getRowDimension() {
            return m.getRowDimension();
        }

        @Override
        public int getColumnDimension() {
            return m.getColumnDimension();
        }

        @Override
        public RealVector operate(RealVector x) {
            return m.operate(x);
        }

    }

}
