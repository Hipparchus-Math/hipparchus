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
package org.hipparchus.analysis.interpolation;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.List;

public class AxisCheckerTest {

    @Test
    public void testMonoDimensional() {

        final List<D1> data = Arrays.asList(new D1(2.0), new D1(1.0), new D1(4.0), new D1(3.0));

        final RegularIndexer indexer = new AxisChecker<>(D1::getX, 1.0e-15).checkGridData(data);

        Assertions.assertEquals(1.0, indexer.getMin(),  1.0e-15);
        Assertions.assertEquals(4.0, indexer.getMax(),  1.0e-15);
        Assertions.assertEquals(  4, indexer.getN());
        Assertions.assertEquals(1.0, indexer.getStep(), 1.0e-15);

        Assertions.assertEquals(0, indexer.index(1.0));
        Assertions.assertEquals(1, indexer.index(2.0));
        Assertions.assertEquals(2, indexer.index(3.0));
        Assertions.assertEquals(3, indexer.index(4.0));

        Assertions.assertEquals(1.0, indexer.coordinate(0), 1.0e-15);
        Assertions.assertEquals(2.0, indexer.coordinate(1), 1.0e-15);
        Assertions.assertEquals(3.0, indexer.coordinate(2), 1.0e-15);
        Assertions.assertEquals(4.0, indexer.coordinate(3), 1.0e-15);

    }

    @Test
    public void testBiDimensional() {

        final List<D2> data = Arrays.asList(new D2(0.0, 2.0), new D2(4.0, 1.0), new D2(2.0, 4.0), new D2(2.0, 2.0),
                                            new D2(2.0, 3.0), new D2(2.0, 1.0), new D2(0.0, 4.0), new D2(0.0, 3.0),
                                            new D2(4.0, 2.0), new D2(0.0, 1.0), new D2(4.0, 4.0), new D2(4.0, 3.0));

        final RegularIndexer xIndexer = new AxisChecker<>(D2::getX, 1.0e-15).checkGridData(data);
        final RegularIndexer yIndexer = new AxisChecker<>(D2::getY, 1.0e-15).checkGridData(data);

        Assertions.assertEquals(0.0, xIndexer.getMin(),  1.0e-15);
        Assertions.assertEquals(4.0, xIndexer.getMax(),  1.0e-15);
        Assertions.assertEquals(  3, xIndexer.getN());
        Assertions.assertEquals(2.0, xIndexer.getStep(), 1.0e-15);

        Assertions.assertEquals(0, xIndexer.index(0.0));
        Assertions.assertEquals(1, xIndexer.index(2.0));
        Assertions.assertEquals(2, xIndexer.index(4.0));

        Assertions.assertEquals(0.0, xIndexer.coordinate(0), 1.0e-15);
        Assertions.assertEquals(2.0, xIndexer.coordinate(1), 1.0e-15);
        Assertions.assertEquals(4.0, xIndexer.coordinate(2), 1.0e-15);

        Assertions.assertEquals(1.0, yIndexer.getMin(),  1.0e-15);
        Assertions.assertEquals(4.0, yIndexer.getMax(),  1.0e-15);
        Assertions.assertEquals(  4, yIndexer.getN());
        Assertions.assertEquals(1.0, yIndexer.getStep(), 1.0e-15);

        Assertions.assertEquals(0, yIndexer.index(1.0));
        Assertions.assertEquals(1, yIndexer.index(2.0));
        Assertions.assertEquals(2, yIndexer.index(3.0));
        Assertions.assertEquals(3, yIndexer.index(4.0));

        Assertions.assertEquals(1.0, yIndexer.coordinate(0), 1.0e-15);
        Assertions.assertEquals(2.0, yIndexer.coordinate(1), 1.0e-15);
        Assertions.assertEquals(3.0, yIndexer.coordinate(2), 1.0e-15);
        Assertions.assertEquals(4.0, yIndexer.coordinate(3), 1.0e-15);

    }

    @Test
    public void testNaN() {
        final List<D1> data = Arrays.asList(new D1(2.0), new D1(1.0), new D1(Double.NaN), new D1(3.0));
        try {
            new AxisChecker<>(D1::getX, 1.0e-15).checkGridData(data);
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assertions.assertEquals(LocalizedCoreFormats.NAN_ELEMENT_AT_INDEX, miae.getSpecifier());
            Assertions.assertEquals(2, miae.getParts()[0]);
        }
    }

    @Test
    public void testDuplicatedNode() {

        final List<D2> data = Arrays.asList(new D2(0.0, 2.0), new D2(4.0, 1.0), new D2(2.0, 4.0), new D2(2.0, 2.0),
                                            new D2(2.0, 3.0), new D2(2.0, 1.0), new D2(0.0, 3.0), new D2(0.0, 3.0),
                                            new D2(4.0, 2.0), new D2(0.0, 1.0), new D2(4.0, 4.0), new D2(4.0, 3.0));

        final RegularIndexer xIndexer = new AxisChecker<>(D2::getX, 1.0e-15).checkGridData(data);
        Assertions.assertEquals(0.0, xIndexer.getMin(),  1.0e-15);
        Assertions.assertEquals(4.0, xIndexer.getMax(),  1.0e-15);
        Assertions.assertEquals(  3, xIndexer.getN());
        Assertions.assertEquals(2.0, xIndexer.getStep(), 1.0e-15);

        Assertions.assertEquals(0, xIndexer.index(0.0));
        Assertions.assertEquals(1, xIndexer.index(2.0));
        Assertions.assertEquals(2, xIndexer.index(4.0));

        try {
            new AxisChecker<>(D2::getY, 1.0e-15).checkGridData(data);
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assertions.assertEquals(LocalizedCoreFormats.IRREGULAR_GRID, miae.getSpecifier());
            Assertions.assertEquals(3, miae.getParts()[0]);
            Assertions.assertEquals(2, miae.getParts()[1]);
            Assertions.assertEquals(4, miae.getParts()[2]);
        }
    }

    @Test
    public void testOffsetPoint() {

        final List<D2> data = Arrays.asList(new D2(0.0,   2.0), new D2(4.0, 1.0), new D2(2.0, 4.0), new D2(2.0, 2.0),
                                            new D2(2.001, 3.0), new D2(2.0, 1.0), new D2(0.0, 3.0), new D2(0.0, 4.0),
                                            new D2(4.0,   2.0), new D2(0.0, 1.0), new D2(4.0, 4.0), new D2(4.0, 3.0));
        try {
            new AxisChecker<>(D2::getX, 1.0e-15).checkGridData(data);
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assertions.assertEquals(LocalizedCoreFormats.MISALIGNED_GRID_POINT, miae.getSpecifier());
            Assertions.assertEquals(1,      miae.getParts()[0]);
            Assertions.assertEquals(2.0,    (Double) miae.getParts()[1], 1.0e-15);
            Assertions.assertEquals(1.0e-3, (Double) miae.getParts()[2], 1.0e-15);
        }
    }

    private static class D1 {

        private final double x;

        public D1(final double x) {
            this.x = x;
        }

        public double getX() {
            return x;
        }

    }

    private static class D2 extends D1 {

        private final double y;

        public D2(final double x, final double y) {
            super(x);
            this.y = y;
        }

        public double getY() {
            return y;
        }

    }

}
