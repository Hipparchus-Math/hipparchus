/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.geometry.partitioning;

import org.hipparchus.geometry.Point;
import org.hipparchus.geometry.Space;
import org.junit.jupiter.api.Assertions;

public class InteriorChecker<S extends Space, P extends Point<S, P>, H extends Hyperplane<S, P, H, I>, I extends SubHyperplane<S, P, H, I>>
    implements BSPTreeVisitor<S, P, H, I> {

    private final P start;
    private final double   tolerance;

    public InteriorChecker(final P start, final double tolerance) {
        this.start     = start;
        this.tolerance = tolerance;
    }

    @Override
    public Order visitOrder(final BSPTree<S, P, H, I> node) {
        return Order.MINUS_PLUS_SUB;
    }

    @Override
    public void visitInternalNode(final BSPTree<S, P, H, I> node) {
        check(node);
    }

    @Override
    public void visitLeafNode(final BSPTree<S, P, H, I> node) {
        check(node);
    }

    private void check(final BSPTree<S, P, H, I> node) {
        final P point = node.pointInsideCell(start);
        final BSPTree<S, P, H, I> convex = node.pruneAroundConvexCell(Boolean.TRUE, Boolean.FALSE, null);
        //System.out.println(start + " " + point);
        Assertions.assertTrue((Boolean) convex.getCell(point, tolerance).getAttribute());
    }

}
