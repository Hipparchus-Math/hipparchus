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
package org.hipparchus.geometry.partitioning;

import org.hipparchus.geometry.Point;
import org.hipparchus.geometry.Space;

/** This interface defines mappers between a space and one of its sub-spaces.

 * <p>Sub-spaces are the lower dimensions subsets of a n-dimensions
 * space. The (n-1)-dimension sub-spaces are specific sub-spaces known
 * as {@link Hyperplane hyperplanes}. This interface can be used regardless
 * of the dimensions differences. As an example, {@link
 * org.hipparchus.geometry.euclidean.threed.Line Line} in 3D
 * implements Embedding&lt; {@link
 * org.hipparchus.geometry.euclidean.threed.Vector3D Vector3D}, {@link
 * org.hipparchus.geometry.euclidean.oned.Vector1D Vector1D}&gt;, i.e. it
 * maps directly dimensions 3 and 1.</p>

 * <p>In the 3D euclidean space, hyperplanes are 2D planes, and the 1D
 * sub-spaces are lines.</p>

 * <p>
 * Note that this interface is <em>not</em> intended to be implemented
 * by Hipparchus users, it is only intended to be implemented
 * within the library itself. New methods may be added even for minor
 * versions, which breaks compatibility for external implementations.
 * </p>

 * @param <S> Type of the embedding space.
 * @param <P> Type of the points in the embedding space.
 * @param <T> Type of the embedded sub-space.
 * @param <Q> Type of the points in the embedded sub-space.

 * @see Hyperplane
 */
public interface Embedding<S extends Space, P extends Point<S>,
                           T extends Space, Q extends Point<T>> {

    /** Transform a space point into a sub-space point.
     * @param point n-dimension point of the space
     * @return (n-1)-dimension point of the sub-space corresponding to
     * the specified space point
     * @see #toSpace
     */
    Q toSubSpace(P point);

    /** Transform a sub-space point into a space point.
     * @param point (n-1)-dimension point of the sub-space
     * @return n-dimension point of the space corresponding to the
     * specified sub-space point
     * @see #toSubSpace
     */
    P toSpace(Q point);

}
