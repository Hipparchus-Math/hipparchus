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

import java.util.Formatter;
import java.util.Locale;

import org.hipparchus.geometry.Point;
import org.hipparchus.geometry.Space;
import org.hipparchus.geometry.euclidean.oned.Euclidean1D;
import org.hipparchus.geometry.euclidean.oned.IntervalsSet;
import org.hipparchus.geometry.euclidean.oned.OrientedPoint;
import org.hipparchus.geometry.euclidean.oned.SubOrientedPoint;
import org.hipparchus.geometry.euclidean.oned.Vector1D;
import org.hipparchus.geometry.euclidean.threed.Euclidean3D;
import org.hipparchus.geometry.euclidean.threed.Plane;
import org.hipparchus.geometry.euclidean.threed.PolyhedronsSet;
import org.hipparchus.geometry.euclidean.threed.SubPlane;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.geometry.euclidean.twod.Euclidean2D;
import org.hipparchus.geometry.euclidean.twod.Line;
import org.hipparchus.geometry.euclidean.twod.PolygonsSet;
import org.hipparchus.geometry.euclidean.twod.SubLine;
import org.hipparchus.geometry.euclidean.twod.Vector2D;
import org.hipparchus.geometry.spherical.oned.ArcsSet;
import org.hipparchus.geometry.spherical.oned.LimitAngle;
import org.hipparchus.geometry.spherical.oned.S1Point;
import org.hipparchus.geometry.spherical.oned.Sphere1D;
import org.hipparchus.geometry.spherical.oned.SubLimitAngle;
import org.hipparchus.geometry.spherical.twod.Circle;
import org.hipparchus.geometry.spherical.twod.S2Point;
import org.hipparchus.geometry.spherical.twod.Sphere2D;
import org.hipparchus.geometry.spherical.twod.SphericalPolygonsSet;
import org.hipparchus.geometry.spherical.twod.SubCircle;

/** Class dumping a string representation of an {@link AbstractRegion}.
 * <p>
 * This class is intended for tests and debug purposes only.
 * </p>
 * @see RegionParser
 */
public class RegionDumper {

    /** Private constructor for a utility class
     */
    private RegionDumper() {
    }

    /** Get a string representation of an {@link ArcsSet}.
     * @param arcsSet region to dump
     * @return string representation of the region
     */
    public static String dump(final ArcsSet arcsSet) {
        final TreeDumper<Sphere1D, S1Point, LimitAngle, SubLimitAngle> visitor =
                new TreeDumper<Sphere1D, S1Point, LimitAngle, SubLimitAngle>("ArcsSet", arcsSet.getTolerance()) {

            /** {@inheritDoc} */
            @Override
            protected void formatHyperplane(final LimitAngle hyperplane) {
                getFormatter().format("%22.15e %b %22.15e",
                                      hyperplane.getLocation().getAlpha(),
                                      hyperplane.isDirect(),
                                      hyperplane.getTolerance());
            }

        };
        arcsSet.getTree(false).visit(visitor);
        return visitor.getDump();
    }

    /** Get a string representation of a {@link SphericalPolygonsSet}.
     * @param sphericalPolygonsSet region to dump
     * @return string representation of the region
     */
    public static String dump(final SphericalPolygonsSet sphericalPolygonsSet) {
        final TreeDumper<Sphere2D, S2Point, Circle, SubCircle> visitor =
                new TreeDumper<Sphere2D, S2Point, Circle, SubCircle>("SphericalPolygonsSet", sphericalPolygonsSet.getTolerance()) {

            /** {@inheritDoc} */
            @Override
            protected void formatHyperplane(final Circle hyperplane) {
                getFormatter().format("%22.15e %22.15e %22.15e %22.15e",
                                      hyperplane.getPole().getX(), hyperplane.getPole().getY(), hyperplane.getPole().getZ(),
                                      hyperplane.getTolerance());
            }

        };
        sphericalPolygonsSet.getTree(false).visit(visitor);
        return visitor.getDump();
    }

    /** Get a string representation of an {@link IntervalsSet}.
     * @param intervalsSet region to dump
     * @return string representation of the region
     */
    public static String dump(final IntervalsSet intervalsSet) {
        final TreeDumper<Euclidean1D, Vector1D, OrientedPoint, SubOrientedPoint> visitor =
                new TreeDumper<Euclidean1D, Vector1D, OrientedPoint, SubOrientedPoint>("IntervalsSet", intervalsSet.getTolerance()) {

            /** {@inheritDoc} */
            @Override
            protected void formatHyperplane(final OrientedPoint hyperplane) {
                getFormatter().format("%22.15e %b %22.15e",
                                      hyperplane.getLocation().getX(),
                                      hyperplane.isDirect(),
                                      hyperplane.getTolerance());
            }

        };
        intervalsSet.getTree(false).visit(visitor);
        return visitor.getDump();
    }

    /** Get a string representation of a {@link PolygonsSet}.
     * @param polygonsSet region to dump
     * @return string representation of the region
     */
    public static String dump(final PolygonsSet polygonsSet) {
        final TreeDumper<Euclidean2D, Vector2D, Line, SubLine> visitor =
                new TreeDumper<Euclidean2D, Vector2D, Line, SubLine>("PolygonsSet", polygonsSet.getTolerance()) {

            /** {@inheritDoc} */
            @Override
            protected void formatHyperplane(final Line hyperplane) {
                final Vector2D p = hyperplane.toSpace(Vector1D.ZERO);
                getFormatter().format("%22.15e %22.15e %22.15e %22.15e",
                                      p.getX(), p.getY(), hyperplane.getAngle(), hyperplane.getTolerance());
            }

        };
        polygonsSet.getTree(false).visit(visitor);
        return visitor.getDump();
    }

    /** Get a string representation of a {@link PolyhedronsSet}.
     * @param polyhedronsSet region to dump
     * @return string representation of the region
     */
    public static String dump(final PolyhedronsSet polyhedronsSet) {
        final TreeDumper<Euclidean3D, Vector3D, Plane, SubPlane> visitor =
                new TreeDumper<Euclidean3D, Vector3D, Plane, SubPlane>("PolyhedronsSet", polyhedronsSet.getTolerance()) {

            /** {@inheritDoc} */
            @Override
            protected void formatHyperplane(final Plane hyperplane) {
                final Vector3D p = hyperplane.toSpace(Vector2D.ZERO);
                getFormatter().format("%22.15e %22.15e %22.15e %22.15e %22.15e %22.15e %22.15e",
                                      p.getX(), p.getY(), p.getZ(),
                                      hyperplane.getNormal().getX(), hyperplane.getNormal().getY(), hyperplane.getNormal().getZ(),
                                      hyperplane.getTolerance());
            }

        };
        polyhedronsSet.getTree(false).visit(visitor);
        return visitor.getDump();
    }

    /** Dumping visitor.
     * @param <S> Type of the space.
     * @param <P> Type of the points in space.
     * @param <H> Type of the hyperplane.
     * @param <I> Type of the sub-hyperplane.
     */
    private abstract static class TreeDumper<S extends Space,
                                             P extends Point<S, P>,
                                             H extends Hyperplane<S, P, H, I>,
                                             I extends SubHyperplane<S, P, H, I>>
        implements BSPTreeVisitor<S, P, H, I> {

        /** Builder for the string representation of the dumped tree. */
        private final StringBuilder dump;

        /** Formatter for strings. */
        private final Formatter formatter;

        /** Current indentation prefix. */
        private String prefix;

        /** Simple constructor.
         * @param type type of the region to dump
         * @param tolerance tolerance of the region
         */
        public TreeDumper(final String type, final double tolerance) {
            this.dump      = new StringBuilder();
            this.formatter = new Formatter(dump, Locale.US);
            this.prefix    = "";
            formatter.format("%s%n", type);
            formatter.format("tolerance %22.15e%n", tolerance);
        }

        /** Get the string representation of the tree.
         * @return string representation of the tree.
         */
        public String getDump() {
            return dump.toString();
        }

        /** Get the formatter to use.
         * @return formatter to use
         */
        protected Formatter getFormatter() {
            return formatter;
        }

        /** Format a string representation of the hyperplane underlying a cut sub-hyperplane.
         * @param hyperplane hyperplane to format
         */
        protected abstract void formatHyperplane(H hyperplane);

        /** {@inheritDoc} */
        @Override
        public Order visitOrder(final BSPTree<S, P, H, I> node) {
            return Order.SUB_MINUS_PLUS;
        }

        /** {@inheritDoc} */
        @Override
        public void visitInternalNode(final BSPTree<S, P, H, I> node) {
            formatter.format("%s %s internal ", prefix, type(node));
            formatHyperplane(node.getCut().getHyperplane());
            formatter.format("%n");
            prefix = prefix + "  ";
        }

        /** {@inheritDoc} */
        @Override
        public void visitLeafNode(final BSPTree<S, P, H, I> node) {
            formatter.format("%s %s leaf %s%n",
                             prefix, type(node), node.getAttribute());
            for (BSPTree<S, P, H, I> n = node;
                 n.getParent() != null && n == n.getParent().getPlus();
                 n = n.getParent()) {
                prefix = prefix.substring(0, prefix.length() - 2);
            }
        }

        /** Get the type of the node.
         * @param node node to check
         * @return "plus " or "minus" depending on the node being the plus or minus
         * child of its parent ("plus " is arbitrarily returned for the root node)
         */
        private String type(final BSPTree<S, P, H, I> node) {
            return (node.getParent() != null && node == node.getParent().getMinus()) ? "minus" : "plus ";
        }

    }

}
