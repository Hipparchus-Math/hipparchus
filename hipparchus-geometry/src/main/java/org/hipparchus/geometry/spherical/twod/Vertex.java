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
package org.hipparchus.geometry.spherical.twod;

/** Spherical polygons boundary vertex.
 * @see SphericalPolygonsSet#getBoundaryLoops()
 * @see Edge
 */
public class Vertex {

    /** Vertex location. */
    private final S2Point location;

    /** Incoming edge. */
    private Edge incoming;

    /** Outgoing edge. */
    private Edge outgoing;

    /** Build a non-processed vertex not owned by any node yet.
     * @param location vertex location
     */
    Vertex(final S2Point location) {
        this.location = location;
        this.incoming = null;
        this.outgoing = null;
    }

    /** Get Vertex location.
     * @return vertex location
     */
    public S2Point getLocation() {
        return location;
    }

    /** Set incoming edge.
     * <p>
     * The circle supporting the incoming edge is automatically bound
     * with the instance.
     * </p>
     * @param incoming incoming edge
     */
    void setIncoming(final Edge incoming) {
        this.incoming = incoming;
    }

    /** Get incoming edge.
     * @return incoming edge
     */
    public Edge getIncoming() {
        return incoming;
    }

    /** Set outgoing edge.
     * <p>
     * The circle supporting the outgoing edge is automatically bound
     * with the instance.
     * </p>
     * @param outgoing outgoing edge
     */
    void setOutgoing(final Edge outgoing) {
        this.outgoing = outgoing;
    }

    /** Get outgoing edge.
     * @return outgoing edge
     */
    public Edge getOutgoing() {
        return outgoing;
    }

}
