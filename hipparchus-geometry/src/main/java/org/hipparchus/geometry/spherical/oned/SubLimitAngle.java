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
package org.hipparchus.geometry.spherical.oned;

import org.hipparchus.geometry.partitioning.AbstractSubHyperplane;
import org.hipparchus.geometry.partitioning.Hyperplane;
import org.hipparchus.geometry.partitioning.Region;

/** This class represents sub-hyperplane for {@link LimitAngle}.
 * <p>Instances of this class are guaranteed to be immutable.</p>
 */
public class SubLimitAngle extends AbstractSubHyperplane<Sphere1D, Sphere1D> {

    /** Simple constructor.
     * @param hyperplane underlying hyperplane
     * @param remainingRegion remaining region of the hyperplane
     */
    public SubLimitAngle(final Hyperplane<Sphere1D> hyperplane,
                         final Region<Sphere1D> remainingRegion) {
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
    protected AbstractSubHyperplane<Sphere1D, Sphere1D> buildNew(final Hyperplane<Sphere1D> hyperplane,
                                                                 final Region<Sphere1D> remainingRegion) {
        return new SubLimitAngle(hyperplane, remainingRegion);
    }

    /** {@inheritDoc} */
    @Override
    public SplitSubHyperplane<Sphere1D> split(final Hyperplane<Sphere1D> hyperplane) {
        final double global = hyperplane.getOffset(((LimitAngle) getHyperplane()).getLocation());
        if (global < -hyperplane.getTolerance())  {
            return new SplitSubHyperplane<Sphere1D>(null, this);
        } else if (global > hyperplane.getTolerance()) {
            return new SplitSubHyperplane<Sphere1D>(this, null);
        } else {
            return new SplitSubHyperplane<Sphere1D>(null, null);
        }
    }

}
