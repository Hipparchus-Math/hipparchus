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

import java.io.IOException;
import java.text.ParseException;
import java.util.StringTokenizer;

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

/** Class parsing a string representation of an {@link AbstractRegion}.
 * <p>
 * This class is intended for tests and debug purposes only.
 * </p>
 * @see RegionDumper
 */
public class RegionParser {

    /** Private constructor for a utility class
     */
    private RegionParser() {
    }

    /** Parse a string representation of an {@link ArcsSet}.
     * @param s string to parse
     * @return parsed region
     * @exception IOException if the string cannot be read
     * @exception ParseException if the string cannot be parsed
     */
    public static ArcsSet parseArcsSet(final String s)
        throws IOException, ParseException {
        final TreeBuilder<Sphere1D, S1Point, LimitAngle, SubLimitAngle> builder =
                new TreeBuilder<Sphere1D, S1Point, LimitAngle, SubLimitAngle>("ArcsSet", s) {

            /** {@inheritDoc} */
            @Override
            protected LimitAngle parseHyperplane() throws ParseException {
                return new LimitAngle(new S1Point(getNumber()), getBoolean(), getNumber());
            }

        };
        return new ArcsSet(builder.getTree(), builder.getTolerance());
    }

    /** Parse a string representation of a {@link SphericalPolygonsSet}.
     * @param s string to parse
     * @return parsed region
     * @exception IOException if the string cannot be read
     * @exception ParseException if the string cannot be parsed
     */
    public static SphericalPolygonsSet parseSphericalPolygonsSet(final String s)
        throws IOException, ParseException {
        final TreeBuilder<Sphere2D, S2Point, Circle, SubCircle> builder =
                new TreeBuilder<Sphere2D, S2Point, Circle, SubCircle>("SphericalPolygonsSet", s) {

            /** {@inheritDoc} */
            @Override
            public Circle parseHyperplane() {
                return new Circle(new Vector3D(getNumber(), getNumber(), getNumber()), getNumber());
            }

        };
        return new SphericalPolygonsSet(builder.getTree(), builder.getTolerance());
    }

    /** Parse a string representation of an {@link IntervalsSet}.
     * @param s string to parse
     * @return parsed region
     * @exception IOException if the string cannot be read
     * @exception ParseException if the string cannot be parsed
     */
    public static IntervalsSet parseIntervalsSet(final String s)
        throws IOException, ParseException {
        final TreeBuilder<Euclidean1D, Vector1D, OrientedPoint, SubOrientedPoint> builder =
                new TreeBuilder<Euclidean1D, Vector1D, OrientedPoint, SubOrientedPoint>("IntervalsSet", s) {

            /** {@inheritDoc} */
            @Override
            public OrientedPoint parseHyperplane() throws ParseException {
                return new OrientedPoint(new Vector1D(getNumber()), getBoolean(), getNumber());
            }

        };
        return new IntervalsSet(builder.getTree(), builder.getTolerance());
    }

    /** Parse a string representation of a {@link PolygonsSet}.
     * @param s string to parse
     * @return parsed region
     * @exception IOException if the string cannot be read
     * @exception ParseException if the string cannot be parsed
     */
    public static PolygonsSet parsePolygonsSet(final String s)
        throws IOException, ParseException {
        final TreeBuilder<Euclidean2D, Vector2D, Line, SubLine> builder =
                new TreeBuilder<Euclidean2D, Vector2D, Line, SubLine>("PolygonsSet", s) {

            /** {@inheritDoc} */
            @Override
            public Line parseHyperplane() {
                return new Line(new Vector2D(getNumber(), getNumber()), getNumber(), getNumber());
            }

        };
        return new PolygonsSet(builder.getTree(), builder.getTolerance());
    }

    /** Parse a string representation of a {@link PolyhedronsSet}.
     * @param s string to parse
     * @return parsed region
     * @exception IOException if the string cannot be read
     * @exception ParseException if the string cannot be parsed
     */
    public static PolyhedronsSet parsePolyhedronsSet(final String s)
        throws IOException, ParseException {
        final TreeBuilder<Euclidean3D, Vector3D, Plane, SubPlane> builder =
                new TreeBuilder<Euclidean3D, Vector3D, Plane, SubPlane>("PolyhedronsSet", s) {

            /** {@inheritDoc} */
            @Override
            public Plane parseHyperplane() {
                return new Plane(new Vector3D(getNumber(), getNumber(), getNumber()),
                                 new Vector3D(getNumber(), getNumber(), getNumber()),
                                 getNumber());
            }

        };
        return new PolyhedronsSet(builder.getTree(), builder.getTolerance());
    }

    /** Local class for building an {@link AbstractRegion} tree.
     * @param <S> Type of the space.
     * @param <P> Type of the points in space.
     * @param <H> Type of the hyperplane.
     * @param <I> Type of the sub-hyperplane.
     */
    private abstract static class TreeBuilder<S extends Space,
                                              P extends Point<S, P>,
                                              H extends Hyperplane<S, P, H, I>,
                                              I extends SubHyperplane<S, P, H, I>> {

        /** Keyword for tolerance. */
        private static final String TOLERANCE = "tolerance";

        /** Keyword for internal nodes. */
        private static final String INTERNAL  = "internal";

        /** Keyword for leaf nodes. */
        private static final String LEAF      = "leaf";

        /** Keyword for plus children trees. */
        private static final String PLUS      = "plus";

        /** Keyword for minus children trees. */
        private static final String MINUS     = "minus";

        /** Keyword for true flags. */
        private static final String TRUE      = "true";

        /** Keyword for false flags. */
        private static final String FALSE     = "false";

        /** Tree root. */
        private BSPTree<S, P, H, I> root;

        /** Tolerance. */
        private final double tolerance;

        /** Tokenizer parsing string representation. */
        private final StringTokenizer tokenizer;

        /** Simple constructor.
         * @param type type of the expected representation
         * @param s string representation
         * @exception IOException if the string cannot be read
         * @exception ParseException if the string cannot be parsed
         */
        public TreeBuilder(final String type, final String s)
            throws IOException, ParseException {
            root = null;
            tokenizer = new StringTokenizer(s);
            getWord(type);
            getWord(TOLERANCE);
            tolerance = getNumber();
            getWord(PLUS);
            root = new BSPTree<>();
            parseTree(root);
            if (tokenizer.hasMoreTokens()) {
                throw new ParseException("unexpected " + tokenizer.nextToken(), 0);
            }
        }

        /** Parse a tree.
         * @param node start node
         * @exception IOException if the string cannot be read
         * @exception ParseException if the string cannot be parsed
         */
        private void parseTree(final BSPTree<S, P, H, I> node)
            throws IOException, ParseException {
            if (INTERNAL.equals(getWord(INTERNAL, LEAF))) {
                // this is an internal node, it has a cut sub-hyperplane (stored as a whole hyperplane)
                // then a minus tree, then a plus tree
                node.insertCut(parseHyperplane());
                getWord(MINUS);
                parseTree(node.getMinus());
                getWord(PLUS);
                parseTree(node.getPlus());
            } else {
                // this is a leaf node, it has only an inside/outside flag
                node.setAttribute(getBoolean());
            }
        }

        /** Get next word.
         * @param allowed allowed values
         * @return parsed word
         * @exception ParseException if the string cannot be parsed
         */
        protected String getWord(final String ... allowed) throws ParseException {
            final String token = tokenizer.nextToken();
            for (final String a : allowed) {
                if (a.equals(token)) {
                    return token;
                }
            }
            throw new ParseException(token + " != " + allowed[0], 0);
        }

        /** Get next number.
         * @return parsed number
         * @exception NumberFormatException if the string cannot be parsed
         */
        protected double getNumber() throws NumberFormatException {
            return Double.parseDouble(tokenizer.nextToken());
        }

        /** Get next boolean.
         * @return parsed boolean
         * @exception ParseException if the string cannot be parsed
         */
        protected boolean getBoolean() throws ParseException {
            return getWord(TRUE, FALSE).equals(TRUE);
        }

        /** Get the built tree.
         * @return built tree
         */
        public BSPTree<S, P, H, I> getTree() {
            return root;
        }

        /** Get the tolerance.
         * @return tolerance
         */
        public double getTolerance() {
            return tolerance;
        }

        /** Parse an hyperplane.
         * @return next hyperplane from the stream
         * @exception IOException if the string cannot be read
         * @exception ParseException if the string cannot be parsed
         */
        protected abstract H parseHyperplane()
            throws IOException, ParseException;

    }

}
