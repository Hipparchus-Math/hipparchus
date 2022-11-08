/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */
package org.hipparchus.geometry.euclidean.threed;

import org.hipparchus.geometry.Point;
import org.hipparchus.geometry.euclidean.oned.Euclidean1D;
import org.hipparchus.geometry.euclidean.oned.Vector1D;
import org.hipparchus.geometry.euclidean.twod.Euclidean2D;
import org.hipparchus.geometry.euclidean.twod.PolygonsSet;
import org.hipparchus.geometry.euclidean.twod.Vector2D;
import org.hipparchus.geometry.partitioning.AbstractSubHyperplane;
import org.hipparchus.geometry.partitioning.BSPTree;
import org.hipparchus.geometry.partitioning.Hyperplane;
import org.hipparchus.geometry.partitioning.Region;
import org.hipparchus.geometry.partitioning.SubHyperplane;

/** This class represents a sub-hyperplane for {@link Plane}.
 */
public class SubPlane extends AbstractSubHyperplane<Euclidean3D, Euclidean2D> {

    /** Simple constructor.
     * @param hyperplane underlying hyperplane
     * @param remainingRegion remaining region of the hyperplane
     */
    public SubPlane(final Hyperplane<Euclidean3D> hyperplane,
                    final Region<Euclidean2D> remainingRegion) {
        super(hyperplane, remainingRegion);
    }

    /** {@inheritDoc} */
    @Override
    protected AbstractSubHyperplane<Euclidean3D, Euclidean2D> buildNew(final Hyperplane<Euclidean3D> hyperplane,
                                                                       final Region<Euclidean2D> remainingRegion) {
        return new SubPlane(hyperplane, remainingRegion);
    }

    /** Split the instance in two parts by an hyperplane.
     * @param hyperplane splitting hyperplane
     * @return an object containing both the part of the instance
     * on the plus side of the instance and the part of the
     * instance on the minus side of the instance
     */
    @Override
    public SplitSubHyperplane<Euclidean3D> split(Hyperplane<Euclidean3D> hyperplane) {

        final Plane otherPlane = (Plane) hyperplane;
        final Plane thisPlane  = (Plane) getHyperplane();
        final Line  inter      = otherPlane.intersection(thisPlane);
        final double tolerance = thisPlane.getTolerance();

        if (inter == null) {
            // the hyperplanes are parallel
            final double global = otherPlane.getOffset(thisPlane);
            if (global < -tolerance) {
                return new SplitSubHyperplane<>(null, this);
            } else if (global > tolerance) {
                return new SplitSubHyperplane<>(this, null);
            } else {
                return new SplitSubHyperplane<>(null, null);
            }
        }

        // the hyperplanes do intersect
        Vector2D p = thisPlane.toSubSpace((Point<Euclidean3D>) inter.toSpace((Point<Euclidean1D>) Vector1D.ZERO));
        Vector2D q = thisPlane.toSubSpace((Point<Euclidean3D>) inter.toSpace((Point<Euclidean1D>) Vector1D.ONE));
        Vector3D crossP = Vector3D.crossProduct(inter.getDirection(), thisPlane.getNormal());
        if (crossP.dotProduct(otherPlane.getNormal()) < 0) {
            final Vector2D tmp = p;
            p           = q;
            q           = tmp;
        }
        final SubHyperplane<Euclidean2D> l2DMinus =
            new org.hipparchus.geometry.euclidean.twod.Line(p, q, tolerance).wholeHyperplane();
        final SubHyperplane<Euclidean2D> l2DPlus =
            new org.hipparchus.geometry.euclidean.twod.Line(q, p, tolerance).wholeHyperplane();

        final BSPTree<Euclidean2D> splitTree = getRemainingRegion().getTree(false).split(l2DMinus);
        final BSPTree<Euclidean2D> plusTree  = getRemainingRegion().isEmpty(splitTree.getPlus()) ?
                                               new BSPTree<>(Boolean.FALSE) :
                                               new BSPTree<>(l2DPlus, new BSPTree<>(Boolean.FALSE),
                                                             splitTree.getPlus(), null);

        final BSPTree<Euclidean2D> minusTree = getRemainingRegion().isEmpty(splitTree.getMinus()) ?
                                               new BSPTree<>(Boolean.FALSE) :
                                               new BSPTree<>(l2DMinus, new BSPTree<>(Boolean.FALSE),
                                                                            splitTree.getMinus(), null);

        return new SplitSubHyperplane<>(new SubPlane(thisPlane.copySelf(), new PolygonsSet(plusTree, tolerance)),
                                        new SubPlane(thisPlane.copySelf(), new PolygonsSet(minusTree, tolerance)));

    }

}
