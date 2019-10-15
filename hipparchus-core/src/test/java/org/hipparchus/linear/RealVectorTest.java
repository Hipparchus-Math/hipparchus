/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */

package org.hipparchus.linear;

import java.util.Iterator;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.linear.RealVector.Entry;
import org.junit.Assert;
import org.junit.Test;

/**
 * Tests for {@link RealVector}.
 */
public class RealVectorTest extends RealVectorAbstractTest {

    @Override
    public RealVector create(final double[] data) {
        return new RealVectorTestImpl(data);
    }

    @Test
    @Override
    public void testAppendVector() {
        checkUnsupported(() -> super.testAppendVector());
    }

    @Test
    @Override
    public void testAppendScalar() {
        checkUnsupported(() -> super.testAppendScalar());
    }

    @Test
    @Override
    public void testGetSubVector() {
        checkUnsupported(() -> super.testGetSubVector());
    }

    @Test
    @Override
    public void testGetSubVectorInvalidIndex1() {
        checkUnsupported(() -> super.testGetSubVectorInvalidIndex1());
    }

    @Test
    @Override
    public void testGetSubVectorInvalidIndex2() {
        checkUnsupported(() -> super.testGetSubVectorInvalidIndex2());
    }

    @Test
    @Override
    public void testGetSubVectorInvalidIndex3() {
        checkUnsupported(() -> super.testGetSubVectorInvalidIndex3());
    }

    @Test
    @Override
    public void testGetSubVectorInvalidIndex4() {
        checkUnsupported(() -> super.testGetSubVectorInvalidIndex4());
    }

    @Test
    @Override
    public void testSetSubVectorSameType() {
        checkUnsupported(() -> super.testSetSubVectorSameType());
    }

    @Test
    @Override
    public void testSetSubVectorMixedType() {
        checkUnsupported(() -> super.testSetSubVectorMixedType());
    }

    @Test
    @Override
    public void testSetSubVectorInvalidIndex1() {
        checkUnsupported(() -> super.testSetSubVectorInvalidIndex1());
    }

    @Test
    @Override
    public void testSetSubVectorInvalidIndex2() {
        checkUnsupported(() -> super.testSetSubVectorInvalidIndex2());
    }

    @Test
    @Override
    public void testSetSubVectorInvalidIndex3() {
        checkUnsupported(() -> super.testSetSubVectorInvalidIndex3());
    }

    @Test
    @Override
    public void testIsNaN() {
        checkUnsupported(() -> super.testIsNaN());
    }

    @Test
    @Override
    public void testIsInfinite() {
        checkUnsupported(() -> super.testIsInfinite());
    }

    @Test
    @Override
    public void testEbeMultiplySameType() {
        checkUnsupported(() -> super.testEbeMultiplySameType());
    }

    @Test
    @Override
    public void testEbeMultiplyMixedTypes() {
        checkUnsupported(() -> super.testEbeMultiplyMixedTypes());
    }

    @Test
    @Override
    public void testEbeMultiplyDimensionMismatch() {
        checkUnsupported(() -> super.testEbeMultiplyDimensionMismatch());
    }

    @Test
    @Override
    public void testEbeDivideSameType() {
        checkUnsupported(() -> super.testEbeDivideSameType());
    }

    @Test
    @Override
    public void testEbeDivideMixedTypes() {
        checkUnsupported(() -> super.testEbeDivideMixedTypes());
    }

    @Test
    @Override
    public void testEbeDivideDimensionMismatch() {
        checkUnsupported(() -> super.testEbeDivideDimensionMismatch());
    }

    @Test
    public void testSparseIterator() {
        /*
         * For non-default values, use x + 1, x + 2, etc... to make sure that
         * these values are really different from x.
         */
        final double x = getPreferredEntryValue();
        final double[] data = {
            x, x + 1d, x, x, x + 2d, x + 3d, x + 4d, x, x, x, x + 5d, x + 6d, x
        };

        RealVector v = create(data);
        Entry e;
        int i = 0;
        final double[] nonDefault = {
            x + 1d, x + 2d, x + 3d, x + 4d, x + 5d, x + 6d
        };
        for (Iterator<Entry> it = v.sparseIterator(); it.hasNext(); i++) {
            e = it.next();
            Assert.assertEquals(nonDefault[i], e.getValue(), 0);
        }
        double [] onlyOne = {x, x + 1d, x};
        v = create(onlyOne);
        for(Iterator<Entry> it = v.sparseIterator(); it.hasNext(); ) {
            e = it.next();
            Assert.assertEquals(onlyOne[1], e.getValue(), 0);
        }
    }

    @Test
    @Override
    public void testSerial() {
        checkUnsupported(() -> super.testSerial());
    }

    @Test
    @Override
    public void testEquals() {
        checkUnsupported(() -> super.testEquals());
    }

    interface Thunk {
        void call();
    }

    private void checkUnsupported(final Thunk t) {
        try {
            t.call();
            Assert.fail("an exception should have been thrown");
        } catch (MathRuntimeException mre) {
            Assert.assertEquals(LocalizedCoreFormats.UNSUPPORTED_OPERATION, mre.getSpecifier());
        } catch (UnsupportedOperationException uoe) {
            // expected
        }
    }

}
