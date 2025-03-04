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
package org.hipparchus.geometry;

import java.io.Serializable;

/** This interface represents a generic geometrical point.
 * @param <S> Type of the space.
 * @param <P> Type of the points in space.
 * @see Space
 * @see Vector
 */
public interface Point<S extends Space, P extends Point<S, P>> extends Serializable {

    /** Get the space to which the point belongs.
     * @return containing space
     */
    Space getSpace();

    /**
     * Returns true if any coordinate of this point is NaN; false otherwise
     * @return  true if any coordinate of this point is NaN; false otherwise
     */
    boolean isNaN();

    /** Compute the distance between the instance and another point.
     * @param p second point
     * @return the distance between the instance and p
     */
    double distance(P p);

    /** Move towards another point.
     * <p>
     * Motion is linear (along space curvature) and based on a ratio
     * where 0.0 stands for not moving at all, 0.5 stands for moving halfway
     * towards other point, and 1.0 stands for moving fully to the other point.
     * </p>
     * @param other other point
     * @param ratio motion ratio,
     * @return moved point
     * @since 4.0
     */
    P moveTowards(P other, double ratio);

}
