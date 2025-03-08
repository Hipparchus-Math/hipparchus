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
package org.hipparchus.geometry.euclidean.threed;

import org.hipparchus.geometry.euclidean.oned.Vector1D;
import org.hipparchus.geometry.euclidean.twod.Euclidean2D;
import org.hipparchus.geometry.euclidean.twod.PolygonsSet;
import org.hipparchus.geometry.euclidean.twod.SubLine;
import org.hipparchus.geometry.euclidean.twod.Vector2D;
import org.hipparchus.geometry.partitioning.AbstractSubHyperplane;
import org.hipparchus.geometry.partitioning.BSPTree;
import org.hipparchus.geometry.partitioning.Region;

/** This class represents a sub-hyperplane for {@link Plane}.
 */
public class SubPlane extends AbstractSubHyperplane<Euclidean3D, Vector3D, Plane, SubPlane,
                                                    Euclidean2D, Vector2D, org.hipparchus.geometry.euclidean.twod.Line, SubLine> {

    /** Simple constructor.
     * @param hyperplane underlying hyperplane
     * @param remainingRegion remaining region of the hyperplane
     */
    public SubPlane(final Plane hyperplane,
                    final Region<Euclidean2D, Vector2D, org.hipparchus.geometry.euclidean.twod.Line, SubLine> remainingRegion) {
        super(hyperplane, remainingRegion);
    }

    /** {@inheritDoc} */
    @Override
    protected SubPlane buildNew(final Plane hyperplane,
                                final Region<Euclidean2D, Vector2D, org.hipparchus.geometry.euclidean.twod.Line, SubLine> remainingRegion) {
        return new SubPlane(hyperplane, remainingRegion);
    }

    /** {@inheritDoc} */
    @Override
    public Vector3D getInteriorPoint() {
        return isEmpty() ? null : getHyperplane().toSpace(getRemainingRegion().getInteriorPoint());
    }

    /** Split the instance in two parts by an hyperplane.
     * @param hyperplane splitting hyperplane
     * @return an object containing both the part of the instance
     * on the plus side of the instance and the part of the
     * instance on the minus side of the instance
     */
    @Override
    public SplitSubHyperplane<Euclidean3D, Vector3D, Plane, SubPlane> split(Plane hyperplane) {

        final Line  inter      = hyperplane.intersection(getHyperplane());
        final double tolerance = getHyperplane().getTolerance();

        if (inter == null) {
            // the hyperplanes are parallel
            final double global = hyperplane.getOffset(getHyperplane());
            if (global < -tolerance) {
                return new SplitSubHyperplane<>(null, this);
            } else if (global > tolerance) {
                return new SplitSubHyperplane<>(this, null);
            } else {
                return new SplitSubHyperplane<>(null, null);
            }
        }

        // the hyperplanes do intersect
        Vector2D p = getHyperplane().toSubSpace(inter.toSpace(Vector1D.ZERO));
        Vector2D q = getHyperplane().toSubSpace(inter.toSpace(Vector1D.ONE));
        Vector3D crossP = Vector3D.crossProduct(inter.getDirection(), getHyperplane().getNormal());
        if (crossP.dotProduct(hyperplane.getNormal()) < 0) {
            final Vector2D tmp = p;
            p           = q;
            q           = tmp;
        }
        final SubLine l2DMinus = new org.hipparchus.geometry.euclidean.twod.Line(p, q, tolerance).wholeHyperplane();
        final SubLine l2DPlus  = new org.hipparchus.geometry.euclidean.twod.Line(q, p, tolerance).wholeHyperplane();

        final BSPTree<Euclidean2D, Vector2D, org.hipparchus.geometry.euclidean.twod.Line, SubLine> splitTree =
                getRemainingRegion().getTree(false).split(l2DMinus);
        final BSPTree<Euclidean2D, Vector2D, org.hipparchus.geometry.euclidean.twod.Line, SubLine> plusTree  =
                getRemainingRegion().isEmpty(splitTree.getPlus()) ?
                                             new BSPTree<>(Boolean.FALSE) :
                                             new BSPTree<>(l2DPlus, new BSPTree<>(Boolean.FALSE), splitTree.getPlus(), null);

        final BSPTree<Euclidean2D, Vector2D, org.hipparchus.geometry.euclidean.twod.Line, SubLine> minusTree =
                getRemainingRegion().isEmpty(splitTree.getMinus()) ?
                                             new BSPTree<>(Boolean.FALSE) :
                                             new BSPTree<>(l2DMinus, new BSPTree<>(Boolean.FALSE), splitTree.getMinus(), null);

        return new SplitSubHyperplane<>(new SubPlane(getHyperplane().copySelf(), new PolygonsSet(plusTree, tolerance)),
                                        new SubPlane(getHyperplane().copySelf(), new PolygonsSet(minusTree, tolerance)));

    }

}
