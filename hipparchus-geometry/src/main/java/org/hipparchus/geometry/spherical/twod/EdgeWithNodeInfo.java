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

package org.hipparchus.geometry.spherical.twod;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.geometry.partitioning.BSPTree;

/** Specialized version of {@link Edge} with tree node information.
 * <p>
 * The tree nodes information is used to set up connection between
 * edges (i.e. combine the end vertex of an edge and the start vertex
 * of the next edge) using topological information, thus avoiding
 * inaccuracies due to angular distance only checks.
 * </p>
 * @since 1.4
 */
class EdgeWithNodeInfo extends Edge {

    /** Node containing edge. */
    private final BSPTree<Sphere2D> node;

    /** Node whose intersection with current node defines start point. */
    private final BSPTree<Sphere2D> startNode;

    /** Node whose intersection with current node defines end point. */
    private final BSPTree<Sphere2D> endNode;

    /** Indicator for completely processed edge. */
    private boolean processed;

    /** Build an edge.
     * @param start start point
     * @param end end point
     * @param length length of the arc (it can be greater than π)
     * @param circle circle supporting the edge
     * @param node node containing the edge
     * @param startNode node whose intersection with current node defines start point
     * @param endNode node whose intersection with current node defines end point
     */
    EdgeWithNodeInfo(final Vertex start, final Vertex end,
                     final double length, final Circle circle,
                     final BSPTree<Sphere2D> node,
                     final BSPTree<Sphere2D> startNode,
                     final BSPTree<Sphere2D> endNode) {
        super(start, end, length, circle);
        this.node      = node;
        this.startNode = startNode;
        this.endNode   = endNode;
        this.processed = false;
    }

    /** Check if two edges follow each other naturally.
     * @param previous candidate previous edge
     * @param next candidate next edge
     * @return true if {@code edge} is a natural follower for instance
     */
    public static boolean areNaturalFollowers(final EdgeWithNodeInfo previous, final EdgeWithNodeInfo next) {
        return next.getStart().getIncoming() == null &&
               previous.endNode              == next.node &&
               previous.node                 == next.startNode &&
               Vector3D.dotProduct(previous.getEnd().getLocation().getVector(),
                                   next.getStart().getLocation().getVector()) > 0.0;
    }

    /** Check if two edges result from a single edged having split by a circle.
     * @param previous candidate previous edge
     * @param next candidate next edge
     * @return true if {@code edge} is a natural follower for instance
     */
    public static boolean resultFromASplit(final EdgeWithNodeInfo previous, final EdgeWithNodeInfo next) {
        return next.getStart().getIncoming()          == null &&
               previous.node.getCut().getHyperplane() == next.node.getCut().getHyperplane() &&
               previous.endNode                       == next.startNode &&
               Vector3D.dotProduct(previous.getEnd().getLocation().getVector(),
                                   next.getStart().getLocation().getVector()) > 0.0;
    }

    /** Set the processed flag.
     * @param processed processed flag to set
     */
    public void setProcessed(final boolean processed) {
        this.processed = processed;
    }

    /** Check if the edge has been processed.
     * @return true if the edge has been processed
     */
    public boolean isProcessed() {
        return processed;
    }

}
