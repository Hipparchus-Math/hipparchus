/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.geometry.euclidean.twod;

import org.hipparchus.Field;
import org.hipparchus.RealFieldElement;
import org.hipparchus.util.Decimal64Field;
import org.junit.Assert;
import org.junit.Test;

public class FieldVector2DTest {

    @Test
    public void testCrossProduct() {
        doTestCrossProduct(Decimal64Field.getInstance());
    }

    @Test
    public void testOrientation() {
        doTestOrientation(Decimal64Field.getInstance());
    }

    private <T extends RealFieldElement<T>> void doTestCrossProduct(final Field<T> field) {
        final double epsilon = 1e-10;

        FieldVector2D<T> p1 = new FieldVector2D<>(field, new Vector2D(1, 1));
        FieldVector2D<T> p2 = new FieldVector2D<>(field, new Vector2D(2, 2));

        FieldVector2D<T> p3 = new FieldVector2D<>(field, new Vector2D(3, 3));
        Assert.assertEquals(0.0, p3.crossProduct(p1, p2).getReal(), epsilon);

        FieldVector2D<T> p4 = new FieldVector2D<>(field, new Vector2D(1, 2));
        Assert.assertEquals(1.0, p4.crossProduct(p1, p2).getReal(), epsilon);

        FieldVector2D<T> p5 = new FieldVector2D<>(field, new Vector2D(2, 1));
        Assert.assertEquals(-1.0, p5.crossProduct(p1, p2).getReal(), epsilon);
    }

    private <T extends RealFieldElement<T>> void doTestOrientation(final Field<T> field) {
        Assert.assertTrue(FieldVector2D.orientation(new FieldVector2D<>(field, new Vector2D(0, 0)),
                                                    new FieldVector2D<>(field, new Vector2D(1, 0)),
                                                    new FieldVector2D<>(field, new Vector2D(1, 1))).getReal() > 0);
        Assert.assertTrue(FieldVector2D.orientation(new FieldVector2D<>(field, new Vector2D(1, 0)),
                                                    new FieldVector2D<>(field, new Vector2D(0, 0)),
                                                    new FieldVector2D<>(field, new Vector2D(1, 1))).getReal() < 0);
        Assert.assertEquals(0.0,
                            FieldVector2D.orientation(new FieldVector2D<>(field, new Vector2D(0, 0)),
                                                      new FieldVector2D<>(field, new Vector2D(1, 0)),
                                                      new FieldVector2D<>(field, new Vector2D(1, 0))).getReal(),
                            1.0e-15);
        Assert.assertEquals(0.0,
                            FieldVector2D.orientation(new FieldVector2D<>(field, new Vector2D(0, 0)),
                                                      new FieldVector2D<>(field, new Vector2D(1, 0)),
                                                      new FieldVector2D<>(field, new Vector2D(2, 0))).getReal(),
                            1.0e-15);
    }

}
