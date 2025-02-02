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


/** This interface represents an inversible affine transform in a space.
 * <p>Inversible affine transform include for example scalings,
 * translations, rotations.</p>

 * <p>Transforms are dimension-specific. The consistency rules between
 * the three {@code apply} methods are the following ones for a
 * transformed defined for dimension D:</p>
 * <ul>
 *   <li>
 *     the transform can be applied to a point in the
 *     D-dimension space using its {@link #apply(Point)}
 *     method
 *   </li>
 *   <li>
 *     the transform can be applied to a (D-1)-dimension
 *     hyperplane in the D-dimension space using its
 *     {@link #apply(Hyperplane)} method
 *   </li>
 *   <li>
 *     the transform can be applied to a (D-2)-dimension
 *     sub-hyperplane in a (D-1)-dimension hyperplane using
 *     its {@link #apply(SubHyperplane, Hyperplane, Hyperplane)}
 *     method
 *   </li>
 * </ul>

 * @param <S> Type of the origin space.
 * @param <P> Type of the points in the origin space.
 * @param <H> Type of the hyperplane in the origin space.
 * @param <I> Type of the sub-hyperplane in the origin space.
 * @param <T> Type of the destination sub-space.
 * @param <Q> Type of the points in the destination sub-space.
 * @param <F> Type of the hyperplane in the destination sub-space.
 * @param <J> Type of the sub-hyperplane in the destination sub-space.

 */
public interface Transform<S extends Space,
                           P extends Point<S>,
                           H extends Hyperplane<S, P, H, I>,
                           I extends SubHyperplane<S, P, H, I>,
                           T extends Space,
                           Q extends Point<T>,
                           F extends Hyperplane<T, Q, F, J>,
                           J extends SubHyperplane<T, Q, F, J>> {

    /** Transform a point of a space.
     * @param point point to transform
     * @return a new object representing the transformed point
     */
    P apply(P point);

    /** Transform an hyperplane of a space.
     * @param hyperplane hyperplane to transform
     * @return a new object representing the transformed hyperplane
     */
    H apply(H hyperplane);

    /** Transform a sub-hyperplane embedded in an hyperplane.
     * @param sub sub-hyperplane to transform
     * @param original hyperplane in which the sub-hyperplane is
     * defined (this is the original hyperplane, the transform has
     * <em>not</em> been applied to it)
     * @param transformed hyperplane in which the sub-hyperplane is
     * defined (this is the transformed hyperplane, the transform
     * <em>has</em> been applied to it)
     * @return a new object representing the transformed sub-hyperplane
     */
    J apply(J sub, H original, H transformed);

}
