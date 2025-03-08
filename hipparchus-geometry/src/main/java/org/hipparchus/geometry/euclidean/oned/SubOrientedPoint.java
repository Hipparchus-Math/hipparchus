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
package org.hipparchus.geometry.euclidean.oned;

import org.hipparchus.geometry.partitioning.AbstractSubHyperplane;
import org.hipparchus.geometry.partitioning.Region;

/** This class represents sub-hyperplane for {@link OrientedPoint}.
 * <p>An hyperplane in 1D is a simple point, its orientation being a
 * boolean.</p>
 */
public class SubOrientedPoint
    extends AbstractSubHyperplane<Euclidean1D, Vector1D, OrientedPoint, SubOrientedPoint,
                                  Euclidean1D, Vector1D, OrientedPoint, SubOrientedPoint> {

    /** Simple constructor.
     * @param hyperplane underlying hyperplane
     * @param remainingRegion remaining region of the hyperplane
     */
    public SubOrientedPoint(final OrientedPoint hyperplane,
                            final Region<Euclidean1D, Vector1D, OrientedPoint, SubOrientedPoint> remainingRegion) {
        super(hyperplane, remainingRegion);
    }

    /** {@inheritDoc} */
    @Override
    public double getSize() {
        return 0;
    }

    /** {@inheritDoc} */
    @Override
    public boolean isEmpty() {
        return false;
    }

    /** {@inheritDoc} */
    @Override
    protected SubOrientedPoint buildNew(final OrientedPoint hyperplane,
                                        final Region<Euclidean1D, Vector1D, OrientedPoint, SubOrientedPoint> remainingRegion) {
        return new SubOrientedPoint(hyperplane, remainingRegion);
    }

    /** {@inheritDoc} */
    @Override
    public Vector1D getInteriorPoint() {
        return getHyperplane().getLocation();
    }

    /** {@inheritDoc} */
    @Override
    public SplitSubHyperplane<Euclidean1D, Vector1D, OrientedPoint, SubOrientedPoint> split(final OrientedPoint hyperplane) {
        final double global = hyperplane.getOffset(getHyperplane().getLocation());
        if (global < -hyperplane.getTolerance()) {
            return new SplitSubHyperplane<>(null, this);
        } else if (global > hyperplane.getTolerance()) {
            return new SplitSubHyperplane<>(this, null);
        } else {
            return new SplitSubHyperplane<>(null, null);
        }
    }

}
