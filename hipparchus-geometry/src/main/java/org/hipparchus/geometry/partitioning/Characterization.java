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

import java.util.ArrayList;
import java.util.List;

import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.geometry.Point;
import org.hipparchus.geometry.Space;

/** Cut sub-hyperplanes characterization with respect to inside/outside cells.
 * @see BoundaryBuilder
 * @param <S> Type of the space.
 * @param <P> Type of the points in space.
 * @param <H> Type of the hyperplane.
 * @param <I> Type of the sub-hyperplane.
 */
class Characterization<S extends Space,
                       P extends Point<S, P>,
                       H extends Hyperplane<S, P, H, I>,
                       I extends SubHyperplane<S, P, H, I>> {

    /** Part of the cut sub-hyperplane that touch outside cells. */
    private I outsideTouching;

    /** Part of the cut sub-hyperplane that touch inside cells. */
    private I insideTouching;

    /** Nodes that were used to split the outside touching part. */
    private final NodesSet<S, P, H, I> outsideSplitters;

    /** Nodes that were used to split the inside touching part. */
    private final NodesSet<S, P, H, I> insideSplitters;

    /** Simple constructor.
     * <p>Characterization consists in splitting the specified
     * sub-hyperplane into several parts lying in inside and outside
     * cells of the tree. The principle is to compute characterization
     * twice for each cut sub-hyperplane in the tree, once on the plus
     * node and once on the minus node. The parts that have the same flag
     * (inside/inside or outside/outside) do not belong to the boundary
     * while parts that have different flags (inside/outside or
     * outside/inside) do belong to the boundary.</p>
     * @param node current BSP tree node
     * @param sub sub-hyperplane to characterize
     */
    Characterization(final BSPTree<S, P, H, I> node, final I sub) {
        outsideTouching  = null;
        insideTouching   = null;
        outsideSplitters = new NodesSet<>();
        insideSplitters  = new NodesSet<>();
        characterize(node, sub, new ArrayList<>());
    }

    /** Filter the parts of an hyperplane belonging to the boundary.
     * <p>The filtering consist in splitting the specified
     * sub-hyperplane into several parts lying in inside and outside
     * cells of the tree. The principle is to call this method twice for
     * each cut sub-hyperplane in the tree, once on the plus node and
     * once on the minus node. The parts that have the same flag
     * (inside/inside or outside/outside) do not belong to the boundary
     * while parts that have different flags (inside/outside or
     * outside/inside) do belong to the boundary.</p>
     * @param node current BSP tree node
     * @param sub sub-hyperplane to characterize
     * @param splitters nodes that did split the current one
     */
    private void characterize(final BSPTree<S, P, H, I> node, final I sub,
                              final List<BSPTree<S, P, H, I>> splitters) {
        if (node.getCut() == null) {
            // we have reached a leaf node
            final boolean inside = (Boolean) node.getAttribute();
            if (inside) {
                addInsideTouching(sub, splitters);
            } else {
                addOutsideTouching(sub, splitters);
            }
        } else {
            final H hyperplane = node.getCut().getHyperplane();
            final SubHyperplane.SplitSubHyperplane<S, P, H, I> split = sub.split(hyperplane);
            switch (split.getSide()) {
            case PLUS:
                characterize(node.getPlus(),  sub, splitters);
                break;
            case MINUS:
                characterize(node.getMinus(), sub, splitters);
                break;
            case BOTH:
                splitters.add(node);
                characterize(node.getPlus(),  split.getPlus(),  splitters);
                characterize(node.getMinus(), split.getMinus(), splitters);
                splitters.remove(splitters.size() - 1);
                break;
            default:
                // this should not happen
                throw MathRuntimeException.createInternalError();
            }
        }
    }

    /** Add a part of the cut sub-hyperplane known to touch an outside cell.
     * @param sub part of the cut sub-hyperplane known to touch an outside cell
     * @param splitters sub-hyperplanes that did split the current one
     */
    private void addOutsideTouching(final I sub, final List<BSPTree<S, P, H, I>> splitters) {
        if (outsideTouching == null) {
            outsideTouching = sub;
        } else {
            outsideTouching = outsideTouching.reunite(sub);
        }
        outsideSplitters.addAll(splitters);
    }

    /** Add a part of the cut sub-hyperplane known to touch an inside cell.
     * @param sub part of the cut sub-hyperplane known to touch an inside cell
     * @param splitters sub-hyperplanes that did split the current one
     */
    private void addInsideTouching(final I sub, final List<BSPTree<S, P, H, I>> splitters) {
        if (insideTouching == null) {
            insideTouching = sub;
        } else {
            insideTouching = insideTouching.reunite(sub);
        }
        insideSplitters.addAll(splitters);
    }

    /** Check if the cut sub-hyperplane touches outside cells.
     * @return true if the cut sub-hyperplane touches outside cells
     */
    public boolean touchOutside() {
        return outsideTouching != null && !outsideTouching.isEmpty();
    }

    /** Get all the parts of the cut sub-hyperplane known to touch outside cells.
     * @return parts of the cut sub-hyperplane known to touch outside cells
     * (may be null or empty)
     */
    public I outsideTouching() {
        return outsideTouching;
    }

    /** Get the nodes that were used to split the outside touching part.
     * <p>
     * Splitting nodes are internal nodes (i.e. they have a non-null
     * cut sub-hyperplane).
     * </p>
     * @return nodes that were used to split the outside touching part
     */
    public NodesSet<S, P, H, I> getOutsideSplitters() {
        return outsideSplitters;
    }

    /** Check if the cut sub-hyperplane touches inside cells.
     * @return true if the cut sub-hyperplane touches inside cells
     */
    public boolean touchInside() {
        return insideTouching != null && !insideTouching.isEmpty();
    }

    /** Get all the parts of the cut sub-hyperplane known to touch inside cells.
     * @return parts of the cut sub-hyperplane known to touch inside cells
     * (may be null or empty)
     */
    public I insideTouching() {
        return insideTouching;
    }

    /** Get the nodes that were used to split the inside touching part.
     * <p>
     * Splitting nodes are internal nodes (i.e. they have a non-null
     * cut sub-hyperplane).
     * </p>
     * @return nodes that were used to split the inside touching part
     */
    public NodesSet<S, P, H, I> getInsideSplitters() {
        return insideSplitters;
    }

}
