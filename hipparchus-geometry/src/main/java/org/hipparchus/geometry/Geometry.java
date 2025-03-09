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
package org.hipparchus.geometry;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;

import java.util.List;

/** Utilities for geometry.
 * @since 4.0
 */
public class Geometry {

    /**
     * Private constructor for a utility class.
     */
    private Geometry() {
        // nothing to do
    }

    /**
     * Compute the barycenter of n points.
     * @param <S> Type of the space.
     * @param <P> Type of the points in space.
     * @param points points generating the barycenter
     * @return barycenter of the points
     */
    public static <S extends Space, P extends Point<S, P>> P barycenter(final List<P> points) {

        // safety check
        if (points.isEmpty()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_OF_ELEMENTS_SHOULD_BE_POSITIVE, 0);
        }

        // compute barycenter by moving from point to point
        P current = points.get(0);
        for (int i = 1; i < points.size(); i++) {
            current = points.get(i).moveTowards(current, ((double) i) / (i + 1));
        }

        return current;

    }

}
