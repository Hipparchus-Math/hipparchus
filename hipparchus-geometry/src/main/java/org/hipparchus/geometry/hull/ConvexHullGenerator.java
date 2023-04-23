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
package org.hipparchus.geometry.hull;

import java.util.Collection;

import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.geometry.Point;
import org.hipparchus.geometry.Space;

/**
 * Interface for convex hull generators.
 *
 * @param <S> Type of the {@link Space}
 * @param <P> Type of the {@link Point}
 *
 * @see <a href="http://en.wikipedia.org/wiki/Convex_hull">Convex Hull (Wikipedia)</a>
 * @see <a href="http://mathworld.wolfram.com/ConvexHull.html">Convex Hull (MathWorld)</a>
 *
 */
public interface ConvexHullGenerator<S extends Space, P extends Point<S>> {

    /**
     * Builds the convex hull from the set of input points.
     *
     * @param points the set of input points
     * @return the convex hull
     * @throws MathIllegalStateException if generator fails to generate a convex hull for
     * the given set of input points
     */
    ConvexHull<S, P> generate(Collection<P> points) throws MathIllegalStateException;
}
