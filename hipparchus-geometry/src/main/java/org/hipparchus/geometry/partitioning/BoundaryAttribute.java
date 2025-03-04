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

/** Class holding boundary attributes.
 * <p>This class is used for the attributes associated with the
 * nodes of region boundary shell trees returned by the {@link
 * Region#getTree(boolean) Region.getTree(includeBoundaryAttributes)}
 * when the boolean {@code includeBoundaryAttributes} parameter is
 * set to {@code true}. It contains the parts of the node cut
 * sub-hyperplane that belong to the boundary.</p>
 * <p>This class is a simple placeholder, it does not provide any
 * processing methods.</p>
 * @param <S> Type of the space.
 * @param <P> Type of the points in space.
 * @param <H> Type of the hyperplane.
 * @param <I> Type of the sub-hyperplane.
 * @see Region#getTree
 */
public class BoundaryAttribute<S extends Space,
                               P extends Point<S, P>,
                               H extends Hyperplane<S, P, H, I>,
                               I extends SubHyperplane<S, P, H, I>> {

    /** Part of the node cut sub-hyperplane that belongs to the
     * boundary and has the outside of the region on the plus side of
     * its underlying hyperplane (may be null).
     */
    private final I plusOutside;

    /** Part of the node cut sub-hyperplane that belongs to the
     * boundary and has the inside of the region on the plus side of
     * its underlying hyperplane (may be null).
     */
    private final I plusInside;

    /** Sub-hyperplanes that were used to split the boundary part. */
    private final NodesSet<S, P, H, I> splitters;

    /** Simple constructor.
     * @param plusOutside part of the node cut sub-hyperplane that
     * belongs to the boundary and has the outside of the region on
     * the plus side of its underlying hyperplane (may be null)
     * @param plusInside part of the node cut sub-hyperplane that
     * belongs to the boundary and has the inside of the region on the
     * plus side of its underlying hyperplane (may be null)
     * @param splitters sub-hyperplanes that were used to
     * split the boundary part (may be null)
     */
    BoundaryAttribute(final I plusOutside,
                      final I plusInside,
                      final NodesSet<S, P, H, I> splitters) {
        this.plusOutside = plusOutside;
        this.plusInside  = plusInside;
        this.splitters   = splitters;
    }

    /** Get the part of the node cut sub-hyperplane that belongs to the
     * boundary and has the outside of the region on the plus side of
     * its underlying hyperplane.
     * @return part of the node cut sub-hyperplane that belongs to the
     * boundary and has the outside of the region on the plus side of
     * its underlying hyperplane
     */
    public I getPlusOutside() {
        return plusOutside;
    }

    /** Get the part of the node cut sub-hyperplane that belongs to the
     * boundary and has the inside of the region on the plus side of
     * its underlying hyperplane.
     * @return part of the node cut sub-hyperplane that belongs to the
     * boundary and has the inside of the region on the plus side of
     * its underlying hyperplane
     */
    public I getPlusInside() {
        return plusInside;
    }

    /** Get the sub-hyperplanes that were used to split the boundary part.
     * @return sub-hyperplanes that were used to split the boundary part
     */
    public NodesSet<S, P, H, I> getSplitters() {
        return splitters;
    }

}
