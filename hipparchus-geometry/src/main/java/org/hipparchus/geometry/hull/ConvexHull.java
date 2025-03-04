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
package org.hipparchus.geometry.hull;

import java.io.Serializable;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.geometry.Point;
import org.hipparchus.geometry.Space;
import org.hipparchus.geometry.partitioning.Hyperplane;
import org.hipparchus.geometry.partitioning.Region;
import org.hipparchus.geometry.partitioning.SubHyperplane;

/**
 * This class represents a convex hull.
 *
 * @param <S> Space type.
 * @param <P> Point type.
 * @param <H> Type of the hyperplane.
 * @param <I> Type of the sub-hyperplane.
*/
public interface ConvexHull<S extends Space, P extends Point<S, P>, H extends Hyperplane<S, P, H, I>, I extends SubHyperplane<S, P, H, I>>
        extends Serializable {

    /**
     * Get the vertices of the convex hull.
     * @return vertices of the convex hull
     */
    P[] getVertices();

    /**
     * Returns a new region that is enclosed by the convex hull.
     * @return the region enclosed by the convex hull
     * @throws MathIllegalArgumentException if the number of vertices is not enough to
     * build a region in the respective space
     */
    Region<S, P, H, I> createRegion() throws MathIllegalArgumentException;
}
