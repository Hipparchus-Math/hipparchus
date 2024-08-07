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
package org.hipparchus.geometry.euclidean.twod.hull;

import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.geometry.euclidean.twod.Vector2D;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Collection;

import static org.junit.jupiter.api.Assertions.assertThrows;

/**
 * Test class for MonotoneChain.
 */
class MonotoneChainTest extends ConvexHullGenerator2DAbstractTest {

    @Override
    protected ConvexHullGenerator2D createConvexHullGenerator(boolean includeCollinearPoints) {
        return new MonotoneChain(includeCollinearPoints);
    }

    // ------------------------------------------------------------------------------

    @Test
    void testConvergenceException() {
        assertThrows(MathIllegalStateException.class, () -> {
            final Collection<Vector2D> points = new ArrayList<Vector2D>();

            points.add(new Vector2D(1, 1));
            points.add(new Vector2D(1, 5));
            points.add(new Vector2D(0, 7));
            points.add(new Vector2D(1, 10));
            points.add(new Vector2D(1, 20));
            points.add(new Vector2D(20, 20));
            points.add(new Vector2D(20, 40));
            points.add(new Vector2D(40, 1));

            @SuppressWarnings("unused")
            final ConvexHull2D hull = new MonotoneChain(true, 2).generate(points);
        });
    }

}
