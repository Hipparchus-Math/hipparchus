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
package org.hipparchus.geometry.euclidean.twod;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.geometry.euclidean.oned.Interval;
import org.hipparchus.geometry.euclidean.oned.IntervalsSet;
import org.hipparchus.geometry.euclidean.oned.Vector1D;
import org.hipparchus.geometry.partitioning.BSPTree;
import org.hipparchus.geometry.partitioning.BSPTreeVisitor;
import org.hipparchus.geometry.partitioning.BoundaryProjection;
import org.hipparchus.geometry.partitioning.Region;
import org.hipparchus.geometry.partitioning.Region.Location;
import org.hipparchus.geometry.partitioning.RegionFactory;
import org.hipparchus.random.RandomVectorGenerator;
import org.hipparchus.random.UncorrelatedRandomVectorGenerator;
import org.hipparchus.random.UniformRandomGenerator;
import org.hipparchus.random.Well1024a;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

class PolygonsSetTest {

    @Test
    void testSimplyConnected() {
        Vector2D[][] vertices = new Vector2D[][] {
            new Vector2D[] {
                new Vector2D(36.0, 22.0),
                new Vector2D(39.0, 32.0),
                new Vector2D(19.0, 32.0),
                new Vector2D( 6.0, 16.0),
                new Vector2D(31.0, 10.0),
                new Vector2D(42.0, 16.0),
                new Vector2D(34.0, 20.0),
                new Vector2D(29.0, 19.0),
                new Vector2D(23.0, 22.0),
                new Vector2D(33.0, 25.0)
            }
        };
        PolygonsSet set = buildSet(vertices);
        assertEquals(Region.Location.OUTSIDE, set.checkPoint(new Vector2D(50.0, 30.0)));
        checkPoints(Region.Location.INSIDE, set, new Vector2D[] {
            new Vector2D(30.0, 15.0),
            new Vector2D(15.0, 20.0),
            new Vector2D(24.0, 25.0),
            new Vector2D(35.0, 30.0),
            new Vector2D(19.0, 17.0)
        });
        checkPoints(Region.Location.OUTSIDE, set, new Vector2D[] {
            new Vector2D(50.0, 30.0),
            new Vector2D(30.0, 35.0),
            new Vector2D(10.0, 25.0),
            new Vector2D(10.0, 10.0),
            new Vector2D(40.0, 10.0),
            new Vector2D(50.0, 15.0),
            new Vector2D(30.0, 22.0)
        });
        checkPoints(Region.Location.BOUNDARY, set, new Vector2D[] {
            new Vector2D(30.0, 32.0),
            new Vector2D(34.0, 20.0)
        });
        checkVertices(set.getVertices(), vertices);
        assertEquals(Location.INSIDE, set.checkPoint(set.getInteriorPoint()));
    }

    @Test
    void testBox() {
        PolygonsSet box = new PolygonsSet(0, 2, -1, 1, 1.0e-10);
        assertEquals(4.0, box.getSize(), 1.0e-10);
        assertEquals(8.0, box.getBoundarySize(), 1.0e-10);
        assertEquals(Location.INSIDE, box.checkPoint(box.getInteriorPoint()));
    }

    @Test
    void testInfinite() {
        PolygonsSet box = new PolygonsSet(new BSPTree<>(Boolean.TRUE), 1.0e-10);
        assertTrue(Double.isInfinite(box.getSize()));
        assertEquals(Location.INSIDE, box.checkPoint(box.getInteriorPoint()));
    }

    @Test
    void testStair() {
        Vector2D[][] vertices = new Vector2D[][] {
            new Vector2D[] {
                new Vector2D( 0.0, 0.0),
                new Vector2D( 0.0, 2.0),
                new Vector2D(-0.1, 2.0),
                new Vector2D(-0.1, 1.0),
                new Vector2D(-0.3, 1.0),
                new Vector2D(-0.3, 1.5),
                new Vector2D(-1.3, 1.5),
                new Vector2D(-1.3, 2.0),
                new Vector2D(-1.8, 2.0),
                new Vector2D(-1.8 - 1.0 / FastMath.sqrt(2.0),
                            2.0 - 1.0 / FastMath.sqrt(2.0))
            }
        };

        PolygonsSet set = buildSet(vertices);
        checkVertices(set.getVertices(), vertices);

        assertEquals(1.1 + 0.95 * FastMath.sqrt(2.0), set.getSize(), 1.0e-10);
        assertEquals(Location.INSIDE, set.checkPoint(set.getInteriorPoint()));

    }

    @Test
    void testEmpty() {
        PolygonsSet empty = (PolygonsSet) new RegionFactory<Euclidean2D, Vector2D, Line, SubLine>().
                            getComplement(new PolygonsSet(1.0e-10));
        assertTrue(empty.isEmpty());
        assertEquals(0, empty.getVertices().length);
        assertEquals(0.0, empty.getBoundarySize(), 1.0e-10);
        assertEquals(0.0, empty.getSize(), 1.0e-10);
        for (double y = -1; y < 1; y += 0.1) {
            for (double x = -1; x < 1; x += 0.1) {
                assertEquals(Double.POSITIVE_INFINITY,
                                    empty.projectToBoundary(new Vector2D(x, y)).getOffset(),
                                    1.0e-10);
            }
        }
        assertNull(empty.getInteriorPoint());
    }

    @Test
    void testFull() {
        PolygonsSet full = new PolygonsSet(1.0e-10);
        assertFalse(full.isEmpty());
        assertEquals(0, full.getVertices().length);
        assertEquals(0.0, full.getBoundarySize(), 1.0e-10);
        assertEquals(Double.POSITIVE_INFINITY, full.getSize(), 1.0e-10);
        for (double y = -1; y < 1; y += 0.1) {
            for (double x = -1; x < 1; x += 0.1) {
                assertEquals(Double.NEGATIVE_INFINITY,
                                    full.projectToBoundary(new Vector2D(x, y)).getOffset(),
                                    1.0e-10);
            }
        }
        assertEquals(Location.INSIDE, full.checkPoint(full.getInteriorPoint()));
    }

    @Test
    void testHole() {
        Vector2D[][] vertices = new Vector2D[][] {
            new Vector2D[] {
                new Vector2D(0.0, 0.0),
                new Vector2D(3.0, 0.0),
                new Vector2D(3.0, 3.0),
                new Vector2D(0.0, 3.0)
            }, new Vector2D[] {
                new Vector2D(1.0, 2.0),
                new Vector2D(2.0, 2.0),
                new Vector2D(2.0, 1.0),
                new Vector2D(1.0, 1.0)
            }
        };
        PolygonsSet set = buildSet(vertices);
        checkPoints(Region.Location.INSIDE, set, new Vector2D[] {
            new Vector2D(0.5, 0.5),
            new Vector2D(1.5, 0.5),
            new Vector2D(2.5, 0.5),
            new Vector2D(0.5, 1.5),
            new Vector2D(2.5, 1.5),
            new Vector2D(0.5, 2.5),
            new Vector2D(1.5, 2.5),
            new Vector2D(2.5, 2.5),
            new Vector2D(0.5, 1.0)
        });
        checkPoints(Region.Location.OUTSIDE, set, new Vector2D[] {
            new Vector2D(1.5, 1.5),
            new Vector2D(3.5, 1.0),
            new Vector2D(4.0, 1.5),
            new Vector2D(6.0, 6.0)
        });
        checkPoints(Region.Location.BOUNDARY, set, new Vector2D[] {
            new Vector2D(1.0, 1.0),
            new Vector2D(1.5, 0.0),
            new Vector2D(1.5, 1.0),
            new Vector2D(1.5, 2.0),
            new Vector2D(1.5, 3.0),
            new Vector2D(3.0, 3.0)
        });
        checkVertices(set.getVertices(), vertices);

        for (double x = -0.999; x < 3.999; x += 0.11) {
            Vector2D v = new Vector2D(x, x + 0.5);
            BoundaryProjection<Euclidean2D, Vector2D> projection = set.projectToBoundary(v);
            assertSame(projection.getOriginal(), v);
            Vector2D p = projection.getProjected();
            if (x < -0.5) {
                assertEquals(0.0,      p.getX(), 1.0e-10);
                assertEquals(0.0,      p.getY(), 1.0e-10);
                assertEquals(+v.distance(Vector2D.ZERO), projection.getOffset(), 1.0e-10);
            } else if (x < 0.5) {
                assertEquals(0.0,      p.getX(), 1.0e-10);
                assertEquals(v.getY(), p.getY(), 1.0e-10);
                assertEquals(-v.getX(), projection.getOffset(), 1.0e-10);
            } else if (x < 1.25) {
                assertEquals(1.0,      p.getX(), 1.0e-10);
                assertEquals(v.getY(), p.getY(), 1.0e-10);
                assertEquals(v.getX() - 1.0, projection.getOffset(), 1.0e-10);
            } else if (x < 2.0) {
                assertEquals(v.getX(), p.getX(), 1.0e-10);
                assertEquals(2.0,      p.getY(), 1.0e-10);
                assertEquals(2.0 - v.getY(), projection.getOffset(), 1.0e-10);
            } else if (x < 3.0) {
                assertEquals(v.getX(), p.getX(), 1.0e-10);
                assertEquals(3.0,      p.getY(), 1.0e-10);
                assertEquals(v.getY() - 3.0, projection.getOffset(), 1.0e-10);
            } else {
                assertEquals(3.0,      p.getX(), 1.0e-10);
                assertEquals(3.0,      p.getY(), 1.0e-10);
                assertEquals(+v.distance(new Vector2D(3, 3)), projection.getOffset(), 1.0e-10);
            }

        }

        assertEquals(Location.INSIDE, set.checkPoint(set.getInteriorPoint()));

    }

    @Test
    void testDisjointPolygons() {
        Vector2D[][] vertices = new Vector2D[][] {
            new Vector2D[] {
                new Vector2D(0.0, 1.0),
                new Vector2D(2.0, 1.0),
                new Vector2D(1.0, 2.0)
            }, new Vector2D[] {
                new Vector2D(4.0, 0.0),
                new Vector2D(5.0, 1.0),
                new Vector2D(3.0, 1.0)
            }
        };
        PolygonsSet set = buildSet(vertices);
        assertEquals(Region.Location.INSIDE, set.checkPoint(new Vector2D(1.0, 1.5)));
        checkPoints(Region.Location.INSIDE, set, new Vector2D[] {
            new Vector2D(1.0, 1.5),
            new Vector2D(4.5, 0.8)
        });
        checkPoints(Region.Location.OUTSIDE, set, new Vector2D[] {
            new Vector2D(1.0, 0.0),
            new Vector2D(3.5, 1.2),
            new Vector2D(2.5, 1.0),
            new Vector2D(3.0, 4.0)
        });
        checkPoints(Region.Location.BOUNDARY, set, new Vector2D[] {
            new Vector2D(1.0, 1.0),
            new Vector2D(3.5, 0.5),
            new Vector2D(0.0, 1.0)
        });
        checkVertices(set.getVertices(), vertices);
        assertEquals(Location.INSIDE, set.checkPoint(set.getInteriorPoint()));

    }

    @Test
    void testOppositeHyperplanes() {
        Vector2D[][] vertices = new Vector2D[][] {
            new Vector2D[] {
                new Vector2D(1.0, 0.0),
                new Vector2D(2.0, 1.0),
                new Vector2D(3.0, 1.0),
                new Vector2D(2.0, 2.0),
                new Vector2D(1.0, 1.0),
                new Vector2D(0.0, 1.0)
            }
        };
        PolygonsSet set = buildSet(vertices);
        checkVertices(set.getVertices(), vertices);
        assertEquals(Location.INSIDE, set.checkPoint(set.getInteriorPoint()));

    }

    @Test
    void testSingularPoint() {
        Vector2D[][] vertices = new Vector2D[][] {
            new Vector2D[] {
                new Vector2D( 0.0,  0.0),
                new Vector2D( 1.0,  0.0),
                new Vector2D( 1.0,  1.0),
                new Vector2D( 0.0,  1.0),
                new Vector2D( 0.0,  0.0),
                new Vector2D(-1.0,  0.0),
                new Vector2D(-1.0, -1.0),
                new Vector2D( 0.0, -1.0)
            }
        };
        PolygonsSet set = buildSet(vertices);
        checkVertices(set.getVertices(), vertices);
        assertEquals(Location.INSIDE, set.checkPoint(set.getInteriorPoint()));

    }

    @Test
    void testLineIntersection() {
        Vector2D[][] vertices = new Vector2D[][] {
            new Vector2D[] {
                new Vector2D( 0.0,  0.0),
                new Vector2D( 2.0,  0.0),
                new Vector2D( 2.0,  1.0),
                new Vector2D( 3.0,  1.0),
                new Vector2D( 3.0,  3.0),
                new Vector2D( 1.0,  3.0),
                new Vector2D( 1.0,  2.0),
                new Vector2D( 0.0,  2.0)
            }
        };
        PolygonsSet set = buildSet(vertices);

        Line l1 = new Line(new Vector2D(-1.5, 0.0), FastMath.PI / 4, 1.0e-10);
        SubLine s1 = set.intersection(l1.wholeHyperplane());
        List<Interval> i1 = ((IntervalsSet) s1.getRemainingRegion()).asList();
        assertEquals(2, i1.size());
        Interval v10 = i1.get(0);
        Vector2D p10Lower = l1.toSpace(new Vector1D(v10.getInf()));
        assertEquals(0.0, p10Lower.getX(), 1.0e-10);
        assertEquals(1.5, p10Lower.getY(), 1.0e-10);
        Vector2D p10Upper = l1.toSpace(new Vector1D(v10.getSup()));
        assertEquals(0.5, p10Upper.getX(), 1.0e-10);
        assertEquals(2.0, p10Upper.getY(), 1.0e-10);
        Interval v11 = i1.get(1);
        Vector2D p11Lower = l1.toSpace(new Vector1D(v11.getInf()));
        assertEquals(1.0, p11Lower.getX(), 1.0e-10);
        assertEquals(2.5, p11Lower.getY(), 1.0e-10);
        Vector2D p11Upper = l1.toSpace(new Vector1D(v11.getSup()));
        assertEquals(1.5, p11Upper.getX(), 1.0e-10);
        assertEquals(3.0, p11Upper.getY(), 1.0e-10);

        Line l2 = new Line(new Vector2D(-1.0, 2.0), 0, 1.0e-10);
        SubLine s2 = set.intersection(l2.wholeHyperplane());
        List<Interval> i2 = ((IntervalsSet) s2.getRemainingRegion()).asList();
        assertEquals(1, i2.size());
        Interval v20 = i2.get(0);
        Vector2D p20Lower = l2.toSpace(new Vector1D(v20.getInf()));
        assertEquals(1.0, p20Lower.getX(), 1.0e-10);
        assertEquals(2.0, p20Lower.getY(), 1.0e-10);
        Vector2D p20Upper = l2.toSpace(new Vector1D(v20.getSup()));
        assertEquals(3.0, p20Upper.getX(), 1.0e-10);
        assertEquals(2.0, p20Upper.getY(), 1.0e-10);
        assertEquals(Location.INSIDE, set.checkPoint(set.getInteriorPoint()));

    }

    @Test
    void testUnlimitedSubHyperplane() {
        Vector2D[][] vertices1 = new Vector2D[][] {
            new Vector2D[] {
                new Vector2D(0.0, 0.0),
                new Vector2D(4.0, 0.0),
                new Vector2D(1.4, 1.5),
                new Vector2D(0.0, 3.5)
            }
        };
        PolygonsSet set1 = buildSet(vertices1);
        Vector2D[][] vertices2 = new Vector2D[][] {
            new Vector2D[] {
                new Vector2D(1.4,  0.2),
                new Vector2D(2.8, -1.2),
                new Vector2D(2.5,  0.6)
            }
        };
        PolygonsSet set2 = buildSet(vertices2);

        PolygonsSet set =
            (PolygonsSet) new RegionFactory<Euclidean2D, Vector2D, Line, SubLine>().union(set1.copySelf(), set2.copySelf());
        checkVertices(set1.getVertices(), vertices1);
        checkVertices(set2.getVertices(), vertices2);
        checkVertices(set.getVertices(), new Vector2D[][] {
            new Vector2D[] {
                new Vector2D(0.0,  0.0),
                new Vector2D(1.6,  0.0),
                new Vector2D(2.8, -1.2),
                new Vector2D(2.6,  0.0),
                new Vector2D(4.0,  0.0),
                new Vector2D(1.4,  1.5),
                new Vector2D(0.0,  3.5)
            }
        });
        assertEquals(Location.INSIDE, set.checkPoint(set.getInteriorPoint()));

    }

    @Test
    void testUnion() {
        Vector2D[][] vertices1 = new Vector2D[][] {
            new Vector2D[] {
                new Vector2D( 0.0,  0.0),
                new Vector2D( 2.0,  0.0),
                new Vector2D( 2.0,  2.0),
                new Vector2D( 0.0,  2.0)
            }
        };
        PolygonsSet set1 = buildSet(vertices1);
        Vector2D[][] vertices2 = new Vector2D[][] {
            new Vector2D[] {
                new Vector2D( 1.0,  1.0),
                new Vector2D( 3.0,  1.0),
                new Vector2D( 3.0,  3.0),
                new Vector2D( 1.0,  3.0)
            }
        };
        PolygonsSet set2 = buildSet(vertices2);
        PolygonsSet set  = (PolygonsSet) new RegionFactory<Euclidean2D, Vector2D, Line, SubLine>().union(set1.copySelf(), set2.copySelf());
        checkVertices(set1.getVertices(), vertices1);
        checkVertices(set2.getVertices(), vertices2);
        checkVertices(set.getVertices(), new Vector2D[][] {
            new Vector2D[] {
                new Vector2D( 0.0,  0.0),
                new Vector2D( 2.0,  0.0),
                new Vector2D( 2.0,  1.0),
                new Vector2D( 3.0,  1.0),
                new Vector2D( 3.0,  3.0),
                new Vector2D( 1.0,  3.0),
                new Vector2D( 1.0,  2.0),
                new Vector2D( 0.0,  2.0)
            }
        });
        checkPoints(Region.Location.INSIDE, set, new Vector2D[] {
            new Vector2D(1.0, 1.0),
            new Vector2D(0.5, 0.5),
            new Vector2D(2.0, 2.0),
            new Vector2D(2.5, 2.5),
            new Vector2D(0.5, 1.5),
            new Vector2D(1.5, 1.5),
            new Vector2D(1.5, 0.5),
            new Vector2D(1.5, 2.5),
            new Vector2D(2.5, 1.5),
            new Vector2D(2.5, 2.5)
        });
        checkPoints(Region.Location.OUTSIDE, set, new Vector2D[] {
            new Vector2D(-0.5, 0.5),
            new Vector2D( 0.5, 2.5),
            new Vector2D( 2.5, 0.5),
            new Vector2D( 3.5, 2.5)
        });
        checkPoints(Region.Location.BOUNDARY, set, new Vector2D[] {
            new Vector2D(0.0, 0.0),
            new Vector2D(0.5, 2.0),
            new Vector2D(2.0, 0.5),
            new Vector2D(2.5, 1.0),
            new Vector2D(3.0, 2.5)
        });

        assertEquals(Location.INSIDE, set.checkPoint(set.getInteriorPoint()));

    }

    @Test
    void testIntersection() {
        Vector2D[][] vertices1 = new Vector2D[][] {
            new Vector2D[] {
                new Vector2D( 0.0,  0.0),
                new Vector2D( 2.0,  0.0),
                new Vector2D( 2.0,  2.0),
                new Vector2D( 0.0,  2.0)
            }
        };
        PolygonsSet set1 = buildSet(vertices1);
        Vector2D[][] vertices2 = new Vector2D[][] {
            new Vector2D[] {
                new Vector2D( 1.0,  1.0),
                new Vector2D( 3.0,  1.0),
                new Vector2D( 3.0,  3.0),
                new Vector2D( 1.0,  3.0)
            }
        };
        PolygonsSet set2 = buildSet(vertices2);
        PolygonsSet set  = (PolygonsSet) new RegionFactory<Euclidean2D, Vector2D, Line, SubLine>().
                           intersection(set1.copySelf(), set2.copySelf());
        checkVertices(set1.getVertices(), vertices1);
        checkVertices(set2.getVertices(), vertices2);
        checkVertices(set.getVertices(), new Vector2D[][] {
            new Vector2D[] {
                new Vector2D( 1.0,  1.0),
                new Vector2D( 2.0,  1.0),
                new Vector2D( 2.0,  2.0),
                new Vector2D( 1.0,  2.0)
            }
        });
        checkPoints(Region.Location.INSIDE, set, new Vector2D[] {
            new Vector2D(1.5, 1.5)
        });
        checkPoints(Region.Location.OUTSIDE, set, new Vector2D[] {
            new Vector2D(0.5, 1.5),
            new Vector2D(2.5, 1.5),
            new Vector2D(1.5, 0.5),
            new Vector2D(0.5, 0.5)
        });
        checkPoints(Region.Location.BOUNDARY, set, new Vector2D[] {
            new Vector2D(1.0, 1.0),
            new Vector2D(2.0, 2.0),
            new Vector2D(1.0, 1.5),
            new Vector2D(1.5, 2.0)
        });
        assertEquals(Location.INSIDE, set.checkPoint(set.getInteriorPoint()));
    }

    @Test
    void testXor() {
        Vector2D[][] vertices1 = new Vector2D[][] {
            new Vector2D[] {
                new Vector2D( 0.0,  0.0),
                new Vector2D( 2.0,  0.0),
                new Vector2D( 2.0,  2.0),
                new Vector2D( 0.0,  2.0)
            }
        };
        PolygonsSet set1 = buildSet(vertices1);
        Vector2D[][] vertices2 = new Vector2D[][] {
            new Vector2D[] {
                new Vector2D( 1.0,  1.0),
                new Vector2D( 3.0,  1.0),
                new Vector2D( 3.0,  3.0),
                new Vector2D( 1.0,  3.0)
            }
        };
        PolygonsSet set2 = buildSet(vertices2);
        PolygonsSet set  = (PolygonsSet) new RegionFactory<Euclidean2D, Vector2D, Line, SubLine>().
                           xor(set1.copySelf(), set2.copySelf());
        checkVertices(set1.getVertices(), vertices1);
        checkVertices(set2.getVertices(), vertices2);
        checkVertices(set.getVertices(), new Vector2D[][] {
            new Vector2D[] {
                new Vector2D( 0.0,  0.0),
                new Vector2D( 2.0,  0.0),
                new Vector2D( 2.0,  1.0),
                new Vector2D( 3.0,  1.0),
                new Vector2D( 3.0,  3.0),
                new Vector2D( 1.0,  3.0),
                new Vector2D( 1.0,  2.0),
                new Vector2D( 0.0,  2.0)
            },
            new Vector2D[] {
                new Vector2D( 1.0,  1.0),
                new Vector2D( 1.0,  2.0),
                new Vector2D( 2.0,  2.0),
                new Vector2D( 2.0,  1.0)
            }
        });
        checkPoints(Region.Location.INSIDE, set, new Vector2D[] {
            new Vector2D(0.5, 0.5),
            new Vector2D(2.5, 2.5),
            new Vector2D(0.5, 1.5),
            new Vector2D(1.5, 0.5),
            new Vector2D(1.5, 2.5),
            new Vector2D(2.5, 1.5),
            new Vector2D(2.5, 2.5)
        });
        checkPoints(Region.Location.OUTSIDE, set, new Vector2D[] {
            new Vector2D(-0.5, 0.5),
            new Vector2D( 0.5, 2.5),
            new Vector2D( 2.5, 0.5),
            new Vector2D( 1.5, 1.5),
            new Vector2D( 3.5, 2.5)
        });
        checkPoints(Region.Location.BOUNDARY, set, new Vector2D[] {
            new Vector2D(1.0, 1.0),
            new Vector2D(2.0, 2.0),
            new Vector2D(1.5, 1.0),
            new Vector2D(2.0, 1.5),
            new Vector2D(0.0, 0.0),
            new Vector2D(0.5, 2.0),
            new Vector2D(2.0, 0.5),
            new Vector2D(2.5, 1.0),
            new Vector2D(3.0, 2.5)
        });
        assertEquals(Location.INSIDE, set.checkPoint(set.getInteriorPoint()));
    }

    @Test
    void testDifference() {
        Vector2D[][] vertices1 = new Vector2D[][] {
            new Vector2D[] {
                new Vector2D( 0.0,  0.0),
                new Vector2D( 2.0,  0.0),
                new Vector2D( 2.0,  2.0),
                new Vector2D( 0.0,  2.0)
            }
        };
        PolygonsSet set1 = buildSet(vertices1);
        Vector2D[][] vertices2 = new Vector2D[][] {
            new Vector2D[] {
                new Vector2D( 1.0,  1.0),
                new Vector2D( 3.0,  1.0),
                new Vector2D( 3.0,  3.0),
                new Vector2D( 1.0,  3.0)
            }
        };
        PolygonsSet set2 = buildSet(vertices2);
        PolygonsSet set  = (PolygonsSet) new RegionFactory<Euclidean2D, Vector2D, Line, SubLine>().
                           difference(set1.copySelf(), set2.copySelf());
        checkVertices(set1.getVertices(), vertices1);
        checkVertices(set2.getVertices(), vertices2);
        checkVertices(set.getVertices(), new Vector2D[][] {
            new Vector2D[] {
                new Vector2D( 0.0,  0.0),
                new Vector2D( 2.0,  0.0),
                new Vector2D( 2.0,  1.0),
                new Vector2D( 1.0,  1.0),
                new Vector2D( 1.0,  2.0),
                new Vector2D( 0.0,  2.0)
            }
        });
        checkPoints(Region.Location.INSIDE, set, new Vector2D[] {
            new Vector2D(0.5, 0.5),
            new Vector2D(0.5, 1.5),
            new Vector2D(1.5, 0.5)
        });
        checkPoints(Region.Location.OUTSIDE, set, new Vector2D[] {
            new Vector2D( 2.5, 2.5),
            new Vector2D(-0.5, 0.5),
            new Vector2D( 0.5, 2.5),
            new Vector2D( 2.5, 0.5),
            new Vector2D( 1.5, 1.5),
            new Vector2D( 3.5, 2.5),
            new Vector2D( 1.5, 2.5),
            new Vector2D( 2.5, 1.5),
            new Vector2D( 2.0, 1.5),
            new Vector2D( 2.0, 2.0),
            new Vector2D( 2.5, 1.0),
            new Vector2D( 2.5, 2.5),
            new Vector2D( 3.0, 2.5)
        });
        checkPoints(Region.Location.BOUNDARY, set, new Vector2D[] {
            new Vector2D(1.0, 1.0),
            new Vector2D(1.5, 1.0),
            new Vector2D(0.0, 0.0),
            new Vector2D(0.5, 2.0),
            new Vector2D(2.0, 0.5)
        });
        assertEquals(Location.INSIDE, set.checkPoint(set.getInteriorPoint()));
    }

    @Test
    void testEmptyDifference() {
        Vector2D[][] vertices1 = new Vector2D[][] {
            new Vector2D[] {
                new Vector2D( 0.5, 3.5),
                new Vector2D( 0.5, 4.5),
                new Vector2D(-0.5, 4.5),
                new Vector2D(-0.5, 3.5)
            }
        };
        PolygonsSet set1 = buildSet(vertices1);
        Vector2D[][] vertices2 = new Vector2D[][] {
            new Vector2D[] {
                new Vector2D( 1.0, 2.0),
                new Vector2D( 1.0, 8.0),
                new Vector2D(-1.0, 8.0),
                new Vector2D(-1.0, 2.0)
            }
        };
        PolygonsSet set2 = buildSet(vertices2);
        assertTrue(new RegionFactory<Euclidean2D, Vector2D, Line, SubLine>().
                   difference(set1.copySelf(), set2.copySelf()).isEmpty());
    }

    @Test
    void testChoppedHexagon() {
        double pi6   = FastMath.PI / 6.0;
        double sqrt3 = FastMath.sqrt(3.0);
        SubLine[] hyp = {
            new Line(new Vector2D(   0.0, 1.0),  5 * pi6, 1.0e-10).wholeHyperplane(),
            new Line(new Vector2D(-sqrt3, 1.0),  7 * pi6, 1.0e-10).wholeHyperplane(),
            new Line(new Vector2D(-sqrt3, 1.0),  9 * pi6, 1.0e-10).wholeHyperplane(),
            new Line(new Vector2D(-sqrt3, 0.0), 11 * pi6, 1.0e-10).wholeHyperplane(),
            new Line(new Vector2D(   0.0, 0.0), 13 * pi6, 1.0e-10).wholeHyperplane(),
            new Line(new Vector2D(   0.0, 1.0),  3 * pi6, 1.0e-10).wholeHyperplane(),
            new Line(new Vector2D(-5.0 * sqrt3 / 6.0, 0.0), 9 * pi6, 1.0e-10).wholeHyperplane()
        };
        hyp[1] = hyp[1].split(hyp[0].getHyperplane()).getMinus();
        hyp[2] = hyp[2].split(hyp[1].getHyperplane()).getMinus();
        hyp[3] = hyp[3].split(hyp[2].getHyperplane()).getMinus();
        hyp[4] = hyp[4].split(hyp[3].getHyperplane()).getMinus().split(hyp[0].getHyperplane()).getMinus();
        hyp[5] = hyp[5].split(hyp[4].getHyperplane()).getMinus().split(hyp[0].getHyperplane()).getMinus();
        hyp[6] = hyp[6].split(hyp[3].getHyperplane()).getMinus().split(hyp[1].getHyperplane()).getMinus();
        BSPTree<Euclidean2D, Vector2D, Line, SubLine> tree = new BSPTree<>(Boolean.TRUE);
        for (int i = hyp.length - 1; i >= 0; --i) {
            tree = new BSPTree<>(hyp[i], new BSPTree<>(Boolean.FALSE), tree, null);
        }
        PolygonsSet set = new PolygonsSet(tree, 1.0e-10);
        SubLine splitter =
            new Line(new Vector2D(-2.0 * sqrt3 / 3.0, 0.0), 9 * pi6, 1.0e-10).wholeHyperplane();
        PolygonsSet slice =
            new PolygonsSet(new BSPTree<>(splitter,
                                          set.getTree(false).split(splitter).getPlus(),
                                          new BSPTree<>(Boolean.FALSE), null),
                            1.0e-10);
        assertEquals(Region.Location.OUTSIDE,
                            slice.checkPoint(new Vector2D(0.1, 0.5)));
        assertEquals(11.0 / 3.0, slice.getBoundarySize(), 1.0e-10);
        assertEquals(Location.INSIDE, slice.checkPoint(slice.getInteriorPoint()));

    }

    @Test
    void testConcentric() {
        double h = FastMath.sqrt(3.0) / 2.0;
        Vector2D[][] vertices1 = new Vector2D[][] {
            new Vector2D[] {
                new Vector2D( 0.00, 0.1 * h),
                new Vector2D( 0.05, 0.1 * h),
                new Vector2D( 0.10, 0.2 * h),
                new Vector2D( 0.05, 0.3 * h),
                new Vector2D(-0.05, 0.3 * h),
                new Vector2D(-0.10, 0.2 * h),
                new Vector2D(-0.05, 0.1 * h)
            }
        };
        PolygonsSet set1 = buildSet(vertices1);
        Vector2D[][] vertices2 = new Vector2D[][] {
            new Vector2D[] {
                new Vector2D( 0.00, 0.0 * h),
                new Vector2D( 0.10, 0.0 * h),
                new Vector2D( 0.20, 0.2 * h),
                new Vector2D( 0.10, 0.4 * h),
                new Vector2D(-0.10, 0.4 * h),
                new Vector2D(-0.20, 0.2 * h),
                new Vector2D(-0.10, 0.0 * h)
            }
        };
        PolygonsSet set2 = buildSet(vertices2);
        assertTrue(set2.contains(set1));
    }

    @Test
    void testBug20040520() {
        BSPTree<Euclidean2D, Vector2D, Line, SubLine> a0 =
            new BSPTree<>(buildSegment(new Vector2D(0.85, -0.05),
                                       new Vector2D(0.90, -0.10)),
                          new BSPTree<>(Boolean.FALSE),
                          new BSPTree<>(Boolean.TRUE),
                          null);
        BSPTree<Euclidean2D, Vector2D, Line, SubLine> a1 =
            new BSPTree<>(buildSegment(new Vector2D(0.85, -0.10),
                                       new Vector2D(0.90, -0.10)),
                          new BSPTree<>(Boolean.FALSE), a0, null);
        BSPTree<Euclidean2D, Vector2D, Line, SubLine> a2 =
            new BSPTree<>(buildSegment(new Vector2D(0.90, -0.05),
                                       new Vector2D(0.85, -0.05)),
                          new BSPTree<>(Boolean.FALSE), a1, null);
        BSPTree<Euclidean2D, Vector2D, Line, SubLine> a3 =
            new BSPTree<>(buildSegment(new Vector2D(0.82, -0.05),
                                       new Vector2D(0.82, -0.08)),
                          new BSPTree<>(Boolean.FALSE),
                          new BSPTree<>(Boolean.TRUE),
                          null);
        BSPTree<Euclidean2D, Vector2D, Line, SubLine> a4 =
            new BSPTree<>(buildHalfLine(new Vector2D(0.85, -0.05),
                                        new Vector2D(0.80, -0.05),
                                        false),
                          new BSPTree<>(Boolean.FALSE), a3, null);
        BSPTree<Euclidean2D, Vector2D, Line, SubLine> a5 =
            new BSPTree<>(buildSegment(new Vector2D(0.82, -0.08),
                                       new Vector2D(0.82, -0.18)),
                          new BSPTree<>(Boolean.FALSE),
                          new BSPTree<>(Boolean.TRUE),
                          null);
        BSPTree<Euclidean2D, Vector2D, Line, SubLine> a6 =
            new BSPTree<>(buildHalfLine(new Vector2D(0.82, -0.18),
                                        new Vector2D(0.85, -0.15),
                                        true),
                          new BSPTree<>(Boolean.FALSE), a5, null);
        BSPTree<Euclidean2D, Vector2D, Line, SubLine> a7 =
            new BSPTree<>(buildHalfLine(new Vector2D(0.85, -0.05),
                                        new Vector2D(0.82, -0.08),
                                        false),
                          a4, a6, null);
        BSPTree<Euclidean2D, Vector2D, Line, SubLine> a8 =
            new BSPTree<>(buildLine(new Vector2D(0.85, -0.25),
                                    new Vector2D(0.85,  0.05)),
                          a2, a7, null);
        BSPTree<Euclidean2D, Vector2D, Line, SubLine> a9 =
            new BSPTree<>(buildLine(new Vector2D(0.90,  0.05),
                                    new Vector2D(0.90, -0.50)),
                          a8, new BSPTree<>(Boolean.FALSE), null);

        BSPTree<Euclidean2D, Vector2D, Line, SubLine> b0 =
            new BSPTree<>(buildSegment(new Vector2D(0.92, -0.12),
                                       new Vector2D(0.92, -0.08)),
                          new BSPTree<>(Boolean.FALSE), new BSPTree<>(Boolean.TRUE),
                          null);
        BSPTree<Euclidean2D, Vector2D, Line, SubLine> b1 =
            new BSPTree<>(buildHalfLine(new Vector2D(0.92, -0.08),
                                        new Vector2D(0.90, -0.10),
                                        true),
                          new BSPTree<>(Boolean.FALSE), b0, null);
        BSPTree<Euclidean2D, Vector2D, Line, SubLine> b2 =
            new BSPTree<>(buildSegment(new Vector2D(0.92, -0.18),
                                       new Vector2D(0.92, -0.12)),
                          new BSPTree<>(Boolean.FALSE), new BSPTree<>(Boolean.TRUE),
                          null);
        BSPTree<Euclidean2D, Vector2D, Line, SubLine> b3 =
            new BSPTree<>(buildSegment(new Vector2D(0.85, -0.15),
                                       new Vector2D(0.90, -0.20)),
                          new BSPTree<>(Boolean.FALSE), b2, null);
        BSPTree<Euclidean2D, Vector2D, Line, SubLine> b4 =
            new BSPTree<>(buildSegment(new Vector2D(0.95, -0.15),
                                       new Vector2D(0.85, -0.05)),
                          b1, b3, null);
        BSPTree<Euclidean2D, Vector2D, Line, SubLine> b5 =
            new BSPTree<>(buildHalfLine(new Vector2D(0.85, -0.05),
                                        new Vector2D(0.85, -0.25),
                                        true),
                          new BSPTree<>(Boolean.FALSE), b4, null);
        BSPTree<Euclidean2D, Vector2D, Line, SubLine> b6 =
            new BSPTree<>(buildLine(new Vector2D(0.0, -1.10),
                                    new Vector2D(1.0, -0.10)),
                          new BSPTree<>(Boolean.FALSE), b5, null);

        PolygonsSet c =
            (PolygonsSet) new RegionFactory<Euclidean2D, Vector2D, Line, SubLine>().
                     union(new PolygonsSet(a9, 1.0e-10), new PolygonsSet(b6, 1.0e-10));

        checkPoints(Region.Location.INSIDE, c, new Vector2D[] {
            new Vector2D(0.83, -0.06),
            new Vector2D(0.83, -0.15),
            new Vector2D(0.88, -0.15),
            new Vector2D(0.88, -0.09),
            new Vector2D(0.88, -0.07),
            new Vector2D(0.91, -0.18),
            new Vector2D(0.91, -0.10)
        });

        checkPoints(Region.Location.OUTSIDE, c, new Vector2D[] {
            new Vector2D(0.80, -0.10),
            new Vector2D(0.83, -0.50),
            new Vector2D(0.83, -0.20),
            new Vector2D(0.83, -0.02),
            new Vector2D(0.87, -0.50),
            new Vector2D(0.87, -0.20),
            new Vector2D(0.87, -0.02),
            new Vector2D(0.91, -0.20),
            new Vector2D(0.91, -0.08),
            new Vector2D(0.93, -0.15)
        });

        checkVertices(c.getVertices(),
                      new Vector2D[][] {
            new Vector2D[] {
                new Vector2D(0.85, -0.15),
                new Vector2D(0.90, -0.20),
                new Vector2D(0.92, -0.18),
                new Vector2D(0.92, -0.08),
                new Vector2D(0.90, -0.10),
                new Vector2D(0.90, -0.05),
                new Vector2D(0.82, -0.05),
                new Vector2D(0.82, -0.18),
            }
        });

    }

    @Test
    void testBug20041003() {

        Line[] l = {
            new Line(new Vector2D(0.0, 0.625000007541172),
                     new Vector2D(1.0, 0.625000007541172), 1.0e-10),
            new Line(new Vector2D(-0.19204433621902645, 0.0),
                     new Vector2D(-0.19204433621902645, 1.0), 1.0e-10),
            new Line(new Vector2D(-0.40303524786887,  0.4248364535319128),
                     new Vector2D(-1.12851149797877, -0.2634107480798909), 1.0e-10),
            new Line(new Vector2D(0.0, 2.0),
                     new Vector2D(1.0, 2.0), 1.0e-10)
        };

        BSPTree<Euclidean2D, Vector2D, Line, SubLine> node1 =
            new BSPTree<>(new SubLine(l[0],
                                      new IntervalsSet(intersectionAbscissa(l[0], l[1]),
                                                       intersectionAbscissa(l[0], l[2]),
                                                       1.0e-10)),
                          new BSPTree<>(Boolean.TRUE),
                          new BSPTree<>(Boolean.FALSE),
                          null);
        BSPTree<Euclidean2D, Vector2D, Line, SubLine> node2 =
            new BSPTree<>(new SubLine(l[1],
                                      new IntervalsSet(intersectionAbscissa(l[1], l[2]),
                                                       intersectionAbscissa(l[1], l[3]),
                                                       1.0e-10)),
                          node1,
                          new BSPTree<>(Boolean.FALSE),
                          null);
        BSPTree<Euclidean2D, Vector2D, Line, SubLine> node3 =
            new BSPTree<>(new SubLine(l[2],
                                      new IntervalsSet(intersectionAbscissa(l[2], l[3]),
                                                       Double.POSITIVE_INFINITY,
                                                       1.0e-10)),
                                     node2,
                          new BSPTree<>(Boolean.FALSE),
                          null);
        BSPTree<Euclidean2D, Vector2D, Line, SubLine> node4 =
            new BSPTree<>(l[3].wholeHyperplane(),
                          node3,
                          new BSPTree<>(Boolean.FALSE),
                          null);

        PolygonsSet set = new PolygonsSet(node4, 1.0e-10);
        assertEquals(0, set.getVertices().length);

    }

    @Test
    void testSqueezedHexa() {
        PolygonsSet set = new PolygonsSet(1.0e-10,
                                          new Vector2D(-6, -4), new Vector2D(-8, -8), new Vector2D(  8, -8),
                                          new Vector2D( 6, -4), new Vector2D(10,  4), new Vector2D(-10,  4));
        assertEquals(Location.OUTSIDE, set.checkPoint(new Vector2D(0, 6)));
    }

    @Test
    void testIssue880Simplified() {

        Vector2D[] vertices1 = new Vector2D[] {
            new Vector2D( 90.13595870833188,  38.33604606376991),
            new Vector2D( 90.14047850603913,  38.34600084496253),
            new Vector2D( 90.11045289492762,  38.36801537312368),
            new Vector2D( 90.10871471476526,  38.36878044144294),
            new Vector2D( 90.10424901707671,  38.374300101757),
            new Vector2D( 90.0979455456843,   38.373578376172475),
            new Vector2D( 90.09081227075944,  38.37526295920463),
            new Vector2D( 90.09081378927135,  38.375193883266434)
        };
        PolygonsSet set1 = new PolygonsSet(1.0e-10, vertices1);
        assertEquals(Location.OUTSIDE, set1.checkPoint(new Vector2D(90.12,  38.32)));
        assertEquals(Location.OUTSIDE, set1.checkPoint(new Vector2D(90.135, 38.355)));

    }

    @Test
    void testIssue880Complete() {
        Vector2D[] vertices1 = new Vector2D[] {
                new Vector2D( 90.08714908223715,  38.370299337260235),
                new Vector2D( 90.08709517675004,  38.3702895991413),
                new Vector2D( 90.08401538704919,  38.368849330127944),
                new Vector2D( 90.08258210430711,  38.367634558585564),
                new Vector2D( 90.08251455106665,  38.36763409247078),
                new Vector2D( 90.08106599752608,  38.36761621664249),
                new Vector2D( 90.08249585300035,  38.36753627557965),
                new Vector2D( 90.09075743352184,  38.35914647644972),
                new Vector2D( 90.09099945896571,  38.35896264724079),
                new Vector2D( 90.09269383800086,  38.34595756121246),
                new Vector2D( 90.09638631543191,  38.3457988093121),
                new Vector2D( 90.09666417351019,  38.34523360999418),
                new Vector2D( 90.1297082145872,  38.337670454923625),
                new Vector2D( 90.12971687748956,  38.337669827794684),
                new Vector2D( 90.1240820219179,  38.34328502001131),
                new Vector2D( 90.13084259656404,  38.34017811765017),
                new Vector2D( 90.13378567942857,  38.33860579180606),
                new Vector2D( 90.13519557833206,  38.33621054663689),
                new Vector2D( 90.13545616732307,  38.33614965452864),
                new Vector2D( 90.13553111202748,  38.33613962818305),
                new Vector2D( 90.1356903436448,  38.33610227127048),
                new Vector2D( 90.13576283227428,  38.33609255422783),
                new Vector2D( 90.13595870833188,  38.33604606376991),
                new Vector2D( 90.1361556630693,  38.3360024198866),
                new Vector2D( 90.13622408795709,  38.335987048115726),
                new Vector2D( 90.13696189099994,  38.33581914328681),
                new Vector2D( 90.13746655304897,  38.33616706665265),
                new Vector2D( 90.13845973716064,  38.33650776167099),
                new Vector2D( 90.13950901827667,  38.3368469456463),
                new Vector2D( 90.14393814424852,  38.337591835857495),
                new Vector2D( 90.14483839716831,  38.337076122362475),
                new Vector2D( 90.14565474433601,  38.33769000964429),
                new Vector2D( 90.14569421179482,  38.3377117256905),
                new Vector2D( 90.14577067124333,  38.33770883625908),
                new Vector2D( 90.14600350631684,  38.337714326520995),
                new Vector2D( 90.14600355139731,  38.33771435193319),
                new Vector2D( 90.14600369112401,  38.33771443882085),
                new Vector2D( 90.14600382486884,  38.33771453466096),
                new Vector2D( 90.14600395205912,  38.33771463904344),
                new Vector2D( 90.14600407214999,  38.337714751520764),
                new Vector2D( 90.14600418462749,  38.337714871611695),
                new Vector2D( 90.14600422249327,  38.337714915811034),
                new Vector2D( 90.14867838361471,  38.34113888210675),
                new Vector2D( 90.14923750157374,  38.341582537502575),
                new Vector2D( 90.14877083250991,  38.34160685841391),
                new Vector2D( 90.14816667319519,  38.34244232585684),
                new Vector2D( 90.14797696744586,  38.34248455284745),
                new Vector2D( 90.14484318014337,  38.34385573215269),
                new Vector2D( 90.14477919958296,  38.3453797747614),
                new Vector2D( 90.14202393306448,  38.34464324839456),
                new Vector2D( 90.14198920640195,  38.344651155237216),
                new Vector2D( 90.14155207025175,  38.34486424263724),
                new Vector2D( 90.1415196143314,  38.344871730519),
                new Vector2D( 90.14128611910814,  38.34500196593859),
                new Vector2D( 90.14047850603913,  38.34600084496253),
                new Vector2D( 90.14045907000337,  38.34601860032171),
                new Vector2D( 90.14039496493928,  38.346223030432384),
                new Vector2D( 90.14037626063737,  38.346240203360026),
                new Vector2D( 90.14030005823724,  38.34646920000705),
                new Vector2D( 90.13799164754806,  38.34903093011013),
                new Vector2D( 90.11045289492762,  38.36801537312368),
                new Vector2D( 90.10871471476526,  38.36878044144294),
                new Vector2D( 90.10424901707671,  38.374300101757),
                new Vector2D( 90.10263482039932,  38.37310041316073),
                new Vector2D( 90.09834601753448,  38.373615053823414),
                new Vector2D( 90.0979455456843,  38.373578376172475),
                new Vector2D( 90.09086514328669,  38.37527884194668),
                new Vector2D( 90.09084931407364,  38.37590801712463),
                new Vector2D( 90.09081227075944,  38.37526295920463),
                new Vector2D( 90.09081378927135,  38.375193883266434)
        };
        PolygonsSet set1 = new PolygonsSet(1.0e-8, vertices1);
        assertEquals(Location.OUTSIDE, set1.checkPoint(new Vector2D(90.0905,  38.3755)));
        assertEquals(Location.INSIDE,  set1.checkPoint(new Vector2D(90.09084, 38.3755)));
        assertEquals(Location.OUTSIDE, set1.checkPoint(new Vector2D(90.0913,  38.3755)));
        assertEquals(Location.INSIDE,  set1.checkPoint(new Vector2D(90.1042,  38.3739)));
        assertEquals(Location.INSIDE,  set1.checkPoint(new Vector2D(90.1111,  38.3673)));
        assertEquals(Location.OUTSIDE, set1.checkPoint(new Vector2D(90.0959,  38.3457)));

        Vector2D[] vertices2 = new Vector2D[] {
                new Vector2D( 90.13067558880044,  38.36977255037573),
                new Vector2D( 90.12907570488,  38.36817308242706),
                new Vector2D( 90.1342774136516,  38.356886880294724),
                new Vector2D( 90.13090330629757,  38.34664392676211),
                new Vector2D( 90.13078571364593,  38.344904617518466),
                new Vector2D( 90.1315602208914,  38.3447185040846),
                new Vector2D( 90.1316336226821,  38.34470643148342),
                new Vector2D( 90.134020944832,  38.340936644972885),
                new Vector2D( 90.13912536387306,  38.335497255122334),
                new Vector2D( 90.1396178806582,  38.334878075552126),
                new Vector2D( 90.14083049696671,  38.33316530644106),
                new Vector2D( 90.14145252901329,  38.33152722916191),
                new Vector2D( 90.1404779335565,  38.32863516047786),
                new Vector2D( 90.14282712131586,  38.327504432532066),
                new Vector2D( 90.14616669875488,  38.3237354115015),
                new Vector2D( 90.14860976050608,  38.315714862457924),
                new Vector2D( 90.14999277782437,  38.3164932507504),
                new Vector2D( 90.15005207194997,  38.316534677663356),
                new Vector2D( 90.15508513859612,  38.31878731691609),
                new Vector2D( 90.15919938519221,  38.31852743183782),
                new Vector2D( 90.16093758658837,  38.31880662005153),
                new Vector2D( 90.16099420184912,  38.318825953291594),
                new Vector2D( 90.1665411125756,  38.31859497874757),
                new Vector2D( 90.16999653861313,  38.32505772048029),
                new Vector2D( 90.17475243391698,  38.32594398441148),
                new Vector2D( 90.17940844844992,  38.327427213761325),
                new Vector2D( 90.20951909541378,  38.330616833491774),
                new Vector2D( 90.2155400467941,  38.331746223670336),
                new Vector2D( 90.21559881391778,  38.33175551425302),
                new Vector2D( 90.21916646426041,  38.332584299620805),
                new Vector2D( 90.23863749852285,  38.34778978875795),
                new Vector2D( 90.25459855175802,  38.357790570608984),
                new Vector2D( 90.25964298227257,  38.356918010203174),
                new Vector2D( 90.26024593994703,  38.361692743151366),
                new Vector2D( 90.26146187570015,  38.36311080550837),
                new Vector2D( 90.26614159359622,  38.36510808579902),
                new Vector2D( 90.26621342936448,  38.36507942500333),
                new Vector2D( 90.26652190211962,  38.36494042196722),
                new Vector2D( 90.26621240678867,  38.365113172030874),
                new Vector2D( 90.26614057102057,  38.365141832826794),
                new Vector2D( 90.26380080055299,  38.3660381760273),
                new Vector2D( 90.26315345241,  38.36670658276421),
                new Vector2D( 90.26251574942881,  38.367490323488084),
                new Vector2D( 90.26247873448426,  38.36755266444749),
                new Vector2D( 90.26234628016698,  38.36787989125406),
                new Vector2D( 90.26214559424784,  38.36945909356126),
                new Vector2D( 90.25861728442555,  38.37200753430875),
                new Vector2D( 90.23905557537864,  38.375405314295904),
                new Vector2D( 90.22517251874075,  38.38984691662256),
                new Vector2D( 90.22549955153215,  38.3911564273979),
                new Vector2D( 90.22434386063355,  38.391476432092134),
                new Vector2D( 90.22147729457276,  38.39134652252034),
                new Vector2D( 90.22142070120117,  38.391349167741964),
                new Vector2D( 90.20665060751588,  38.39475580900313),
                new Vector2D( 90.20042268367109,  38.39842558622888),
                new Vector2D( 90.17423771242085,  38.402727751805344),
                new Vector2D( 90.16756796257476,  38.40913898597597),
                new Vector2D( 90.16728283954308,  38.411255399912875),
                new Vector2D( 90.16703538220418,  38.41136059866693),
                new Vector2D( 90.16725865657685,  38.41013618805954),
                new Vector2D( 90.16746107640665,  38.40902614307544),
                new Vector2D( 90.16122795307462,  38.39773101873203)
        };
        PolygonsSet set2 = new PolygonsSet(1.0e-8, vertices2);
        PolygonsSet set  = (PolygonsSet) new RegionFactory<Euclidean2D, Vector2D, Line, SubLine>().
                           difference(set1.copySelf(), set2.copySelf());

        Vector2D[][] vertices = set.getVertices();
        assertNotNull(vertices[0][0]);
        assertEquals(1, vertices.length);
    }

    @Test
    void testTooThinBox() {
        assertEquals(0.0,
                            new PolygonsSet(0.0, 0.0, 0.0, 10.3206397147574, 1.0e-10).getSize(),
                            1.0e-10);
    }

    @Test
    void testWrongUsage() {
        // the following is a wrong usage of the constructor.
        // as explained in the javadoc, the failure is NOT detected at construction
        // time but occurs later on
        PolygonsSet ps = new PolygonsSet(new BSPTree<>(), 1.0e-10);
        assertNotNull(ps);
        try {
            ps.getSize();
            fail("an exception should have been thrown");
        } catch (NullPointerException npe) {
            // this is expected
        }
    }

    @Test
    void testIssue1162() {
        PolygonsSet p = new PolygonsSet(1.0e-10,
                                                new Vector2D(4.267199999996532, -11.928637756014894),
                                                new Vector2D(4.267200000026445, -14.12360595809307),
                                                new Vector2D(9.144000000273694, -14.12360595809307),
                                                new Vector2D(9.144000000233383, -11.928637756020067));

        PolygonsSet w = new PolygonsSet(1.0e-10,
                                                new Vector2D(2.56735636510452512E-9, -11.933116461089332),
                                                new Vector2D(2.56735636510452512E-9, -12.393225665247766),
                                                new Vector2D(2.56735636510452512E-9, -27.785625665247778),
                                                new Vector2D(4.267200000030211,      -27.785625665247778),
                                                new Vector2D(4.267200000030211,      -11.933116461089332));

        assertFalse(p.contains(w));

    }

    @Test
    void testThinRectangle() {

        RegionFactory<Euclidean2D, Vector2D, Line, SubLine> factory = new RegionFactory<>();
        Vector2D pA = new Vector2D(0.0,        1.0);
        Vector2D pB = new Vector2D(0.0,        0.0);
        Vector2D pC = new Vector2D(1.0 / 64.0, 0.0);
        Vector2D pD = new Vector2D(1.0 / 64.0, 1.0);

        // if tolerance is smaller than rectangle width, the rectangle is computed accurately
        Line[] h1 = new Line[] {
            new Line(pA, pB, 1.0 / 256),
            new Line(pB, pC, 1.0 / 256),
            new Line(pC, pD, 1.0 / 256),
            new Line(pD, pA, 1.0 / 256)
        };
        Region<Euclidean2D, Vector2D, Line, SubLine> accuratePolygon = factory.buildConvex(h1);
        assertEquals(1.0 / 64.0, accuratePolygon.getSize(), 1.0e-10);
        assertTrue(Double.isInfinite(new RegionFactory<Euclidean2D, Vector2D, Line, SubLine>().
                                     getComplement(accuratePolygon).getSize()));
        assertEquals(2 * (1.0 + 1.0 / 64.0), accuratePolygon.getBoundarySize(), 1.0e-10);

        // if tolerance is larger than rectangle width, the rectangle degenerates
        // as of 3.3, its two long edges cannot be distinguished anymore and this part of the test did fail
        // this has been fixed in 3.4 (issue MATH-1174)
        Line[] h2 = new Line[] {
            new Line(pA, pB, 1.0 / 16),
            new Line(pB, pC, 1.0 / 16),
            new Line(pC, pD, 1.0 / 16),
            new Line(pD, pA, 1.0 / 16)
        };
        Region<Euclidean2D, Vector2D, Line, SubLine> degeneratedPolygon = factory.buildConvex(h2);
        assertEquals(0.0, degeneratedPolygon.getSize(), 1.0e-10);
        assertTrue(degeneratedPolygon.isEmpty());

    }

    @Test
    void testInconsistentHyperplanes() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            double tolerance = 1.0e-10;
            new RegionFactory<Euclidean2D, Vector2D, Line, SubLine>().
                    buildConvex(new Line(new Vector2D(0, 0), new Vector2D(0, 1), tolerance),
                new Line(new Vector2D(1, 1), new Vector2D(1, 0), tolerance));
        });
    }

    @Test
    void testBoundarySimplification() {

        // a simple square will result in a 4 cuts and 5 leafs tree
        PolygonsSet square = new PolygonsSet(1.0e-10,
                                             new Vector2D(0, 0),
                                             new Vector2D(1, 0),
                                             new Vector2D(1, 1),
                                             new Vector2D(0, 1));
        Vector2D[][] squareBoundary = square.getVertices();
        assertEquals(1, squareBoundary.length);
        assertEquals(4, squareBoundary[0].length);
        Counter squareCount = new Counter();
        squareCount.count(square);
        assertEquals(4, squareCount.getInternalNodes());
        assertEquals(5, squareCount.getLeafNodes());

        // splitting the square in two halves increases the BSP tree
        // with 3 more cuts and 3 more leaf nodes
        SubLine cut = new Line(new Vector2D(0.5, 0.5), 0.0, square.getTolerance()).wholeHyperplane();
        PolygonsSet splitSquare = new PolygonsSet(square.getTree(false).split(cut),
                                                  square.getTolerance());
        Counter splitSquareCount = new Counter();
        splitSquareCount.count(splitSquare);
        assertEquals(squareCount.getInternalNodes() + 3, splitSquareCount.getInternalNodes());
        assertEquals(squareCount.getLeafNodes()     + 3, splitSquareCount.getLeafNodes());

        // the number of vertices should not change, as the intermediate vertices
        // at (0.0, 0.5) and (1.0, 0.5) induced by the top level horizontal split
        // should be removed during the boundary extraction process
        Vector2D[][] splitBoundary = splitSquare.getVertices();
        assertEquals(1, splitBoundary.length);
        assertEquals(4, splitBoundary[0].length);

    }

    @Test
    void testOppositeEdges() {
        PolygonsSet polygon = new PolygonsSet(1.0e-6,
                                              new Vector2D(+1, -2),
                                              new Vector2D(+1,  0),
                                              new Vector2D(+2,  0),
                                              new Vector2D(-1, +2),
                                              new Vector2D(-1,  0),
                                              new Vector2D(-2,  0));
        assertEquals(6.0, polygon.getSize(), 1.0e-10);
        assertEquals(2.0 * FastMath.sqrt(13.0) + 6.0, polygon.getBoundarySize(), 1.0e-10);
    }

    @Test
    void testInfiniteQuadrant() {
        final double tolerance = 1.0e-10;
        BSPTree<Euclidean2D, Vector2D, Line, SubLine> bsp = new BSPTree<>();
        bsp.insertCut(new Line(Vector2D.ZERO, 0.0, tolerance));
        bsp.getPlus().setAttribute(Boolean.FALSE);
        bsp.getMinus().insertCut(new Line(Vector2D.ZERO, MathUtils.SEMI_PI, tolerance));
        bsp.getMinus().getPlus().setAttribute(Boolean.FALSE);
        bsp.getMinus().getMinus().setAttribute(Boolean.TRUE);
        PolygonsSet polygons = new PolygonsSet(bsp, tolerance);
        assertEquals(Double.POSITIVE_INFINITY, polygons.getSize(), 1.0e-10);
    }

    @Test
    void testZigZagBoundaryOversampledIssue46() {
        final double tol = 1.0e-4;
        // sample region, non-convex, not too big, not too small
        final List<Vector2D> vertices = Arrays.asList(
                new Vector2D(-0.12630940610562444e1, (0.8998192093789258 - 0.89) * 100),
                new Vector2D(-0.12731320182988207e1, (0.8963735568774486 - 0.89) * 100),
                new Vector2D(-0.1351107624622557e1, (0.8978258663483273 - 0.89) * 100),
                new Vector2D(-0.13545331405131725e1, (0.8966781238246179 - 0.89) * 100),
                new Vector2D(-0.14324883017454967e1, (0.8981309629283796 - 0.89) * 100),
                new Vector2D(-0.14359875625524995e1, (0.896983965573036 - 0.89) * 100),
                new Vector2D(-0.14749650541159384e1, (0.8977109994666864 - 0.89) * 100),
                new Vector2D(-0.14785037758231825e1, (0.8965644005442432 - 0.89) * 100),
                new Vector2D(-0.15369807257448784e1, (0.8976550608135502 - 0.89) * 100),
                new Vector2D(-0.1526225554339386e1, (0.9010934265410458 - 0.89) * 100),
                new Vector2D(-0.14679028466684121e1, (0.9000043396997698 - 0.89) * 100),
                new Vector2D(-0.14643807494172612e1, (0.9011511073761742 - 0.89) * 100),
                new Vector2D(-0.1386609051963748e1, (0.8996991539048602 - 0.89) * 100),
                new Vector2D(-0.13831601655974668e1, (0.9008466623902937 - 0.89) * 100),
                new Vector2D(-0.1305365419828323e1, (0.8993961857946309 - 0.89) * 100),
                new Vector2D(-0.1301989630405964e1, (0.9005444294061787 - 0.89) * 100));
        Collections.reverse(vertices);
        PolygonsSet expected = new PolygonsSet(tol, vertices.toArray(new Vector2D[0]));
        // sample high resolution boundary
        List<Vector2D> points = new ArrayList<>();
        final Vector2D[] boundary = expected.getVertices()[0];
        double step = tol / 10;
        for (int i = 1; i < boundary.length; i++) {
            Segment edge = new Segment(boundary[i - 1], boundary[i], tol);
            final double length = edge.getLength();
            final double x0 = edge.getLine().toSubSpace(edge.getStart()).getX();
            int n = (int) (length / step);
            for (int j = 0; j < n; j++) {
                points.add(edge.getLine().getPointAt(new Vector1D(x0 + j * step), 0));
            }
        }
        // create zone from high resolution boundary
        PolygonsSet zone = new PolygonsSet(tol, points.toArray(new Vector2D[0]));
        assertEquals(expected.getSize(), zone.getSize(), tol);
        assertEquals(0, expected.getBarycenter().distance(zone.getBarycenter()), tol);
        assertEquals(Location.INSIDE, zone.checkPoint(zone.getBarycenter()), "" + expected.getBarycenter() +  zone.getBarycenter());
        // extra tolerance at corners due to SPS tolerance being a hyperplaneThickness
        // a factor of 3.1 corresponds to a edge intersection angle of ~19 degrees
        final double cornerTol = 3.1 * tol;
        for (Vector2D vertex : vertices) {
            // check original points are on the boundary
            assertEquals(Location.BOUNDARY, zone.checkPoint(vertex), "" + vertex);
            double offset = FastMath.abs(zone.projectToBoundary(vertex).getOffset());
            assertEquals(0, offset, cornerTol, vertex + " offset: " + offset);
        }
    }

    @Test
    void testPositiveQuadrantByVerticesDetailIssue46() {
        double tol = 0.01;
        Line x = new Line(Vector2D.ZERO, new Vector2D(1, 0), tol);
        Line y = new Line(Vector2D.ZERO, new Vector2D(0, 1), tol);
        double length = 1;
        double step = tol / 10;
        // sample high resolution boundary
        int n = (int) (length / step);
        List<Vector2D> points = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            double t = i * step;
            points.add(y.getPointAt(new Vector1D(t), 0));
        }
        for (int i = 0; i < n; i++) {
            double t = i * step;
            points.add(x.getPointAt(new Vector1D(t), 0).add(new Vector2D(0, 1)));
        }
        for (int i = 0; i < n; i++) {
            double t = i * step;
            points.add(new Vector2D(1, 1).subtract(y.getPointAt(new Vector1D(t), 0)));
        }
        Collections.reverse(points);
        PolygonsSet set = new PolygonsSet(tol, points.toArray(new Vector2D[0]));
        RandomVectorGenerator random =
                new UncorrelatedRandomVectorGenerator(2, new UniformRandomGenerator(new Well1024a(0xb8fc5acc91044308L)));
        /* Where exactly the boundaries fall depends on which points are kept from
         * decimation, which can vary by up to tol. So a point up to 2*tol away from a
         * input point may be on the boundary. All input points are guaranteed to be on
         * the boundary, just not the center of the boundary.
         */
        for (int i = 0; i < 1000; ++i) {
            Vector2D v = new Vector2D(random.nextVector());
            final Location actual = set.checkPoint(v);
            if ((v.getX() > tol) && (v.getY() > tol) && (v.getX() < 1 - tol) && (v.getY() < 1 - tol)) {
                if ((v.getX() > 2 * tol) && (v.getY() > 2 * tol) && (v.getX() < 1 - 2 * tol) && (v.getY() < 1 - 2 * tol)) {
                    // certainly inside
                    assertEquals(Location.INSIDE, actual, "" + v);
                } else {
                    // may be inside or boundary
                    assertNotEquals(Location.OUTSIDE, actual, "" + v);
                }
            } else if ((v.getX() < 0) || (v.getY() < 0) || (v.getX() > 1) || (v.getY() > 1)) {
                if ((v.getX() < -tol) || (v.getY() < -tol) || (v.getX() > 1 + tol) || (v.getY() > 1 + tol)) {
                    // certainly outside
                    assertEquals(Location.OUTSIDE, actual);
                } else {
                    // may be outside or boundary
                    assertNotEquals(Location.INSIDE, actual);
                }
            } else {
                // certainly boundary
                assertEquals(Location.BOUNDARY, actual);
            }
        }
        // all input points are on the boundary
        for (Vector2D point : points) {
            assertEquals(Location.BOUNDARY, set.checkPoint(point), "" + point);
        }
    }

    @Test
    void testOversampledCircleIssue64() {
        // setup
        double tol = 1e-2;
        double length = FastMath.PI * 2;
        int n = (int) (length / tol);
        double step = length / n;
        List<Vector2D> points = new ArrayList<>(n);
        for (int i = 0; i < n; i++) {
            double angle = i * step;
            points.add(new Vector2D(FastMath.cos(angle), FastMath.sin(angle)));
        }

        // action
        PolygonsSet set = new PolygonsSet(tol, points.toArray(new Vector2D[0]));

        // verify
        assertEquals(FastMath.PI, set.getSize(), tol);
        assertEquals(0, set.getBarycenter().distance(Vector2D.ZERO), tol);
        // each segment is shorter than boundary, so error builds up
        assertEquals(length, set.getBoundarySize(), 2 * tol);
        for (Vector2D point : points) {
            assertEquals(Location.BOUNDARY, set.checkPoint(point));
        }
    }

    @Test
    public void testWholeSpace() {
        PolygonsSet set = new PolygonsSet(1.0e-10, new Vector2D[0]);
        assertEquals(Double.POSITIVE_INFINITY, set.getSize());
        assertTrue(set.isFull());
        assertFalse(set.isEmpty());
    }

    @Test
    public void testOpenLoops() {
        PolygonsSet negativeXHalfPlane =
                new PolygonsSet(Collections.singleton(new SubLine(new Line(Vector2D.ZERO, MathUtils.SEMI_PI, 1.0e-10),
                                                                  new IntervalsSet(1.0e-10))),
                                1.0e-10);
        PolygonsSet positiveYHalfPlane =
                new PolygonsSet(Collections.singleton(new SubLine(new Line(Vector2D.ZERO, 0.0, 1.0e-10),
                                                                  new IntervalsSet(1.0e-10))),
                                1.0e-10);
        RegionFactory<Euclidean2D, Vector2D, Line, SubLine> factory = new RegionFactory<>();
        PolygonsSet xNegYposQuadrant = (PolygonsSet) factory.intersection(negativeXHalfPlane, positiveYHalfPlane);
        assertEquals(Location.OUTSIDE, xNegYposQuadrant.checkPoint(new Vector2D( 1,  1)));
        assertEquals(Location.INSIDE,  xNegYposQuadrant.checkPoint(new Vector2D(-1,  1)));
        assertEquals(Location.OUTSIDE, xNegYposQuadrant.checkPoint(new Vector2D(-1, -1)));
        assertEquals(Location.OUTSIDE, xNegYposQuadrant.checkPoint(new Vector2D( 1, -1)));

        Vector2D[][] vertices = xNegYposQuadrant.getVertices();
        assertEquals(1, vertices.length);
        assertEquals(5, vertices[0].length);
        assertNull(vertices[0][0]);
        assertEquals(-1.0, vertices[0][1].getX(), 1.0e-15);
        assertEquals( 0.0, vertices[0][1].getY(), 1.0e-15);
        assertEquals( 0.0, vertices[0][2].getX(), 1.0e-15);
        assertEquals( 0.0, vertices[0][2].getY(), 1.0e-15);
        assertEquals( 0.0, vertices[0][3].getX(), 1.0e-15);
        assertEquals( 1.0, vertices[0][3].getY(), 1.0e-15);
        assertNull(vertices[0][4]);
    }

    private static class Counter {

        private int internalNodes;
        private int leafNodes;

        public void count(PolygonsSet polygonsSet) {
            leafNodes     = 0;
            internalNodes = 0;
            polygonsSet.getTree(false).visit(new BSPTreeVisitor<Euclidean2D, Vector2D, Line, SubLine>() {
                public Order visitOrder(BSPTree<Euclidean2D, Vector2D, Line, SubLine> node) {
                    return Order.SUB_PLUS_MINUS;
                }
                public void visitInternalNode(BSPTree<Euclidean2D, Vector2D, Line, SubLine> node) {
                    ++internalNodes;
                }
                public void visitLeafNode(BSPTree<Euclidean2D, Vector2D, Line, SubLine> node) {
                    ++leafNodes;
                }

            });
        }

        public int getInternalNodes() {
            return internalNodes;
        }

        public int getLeafNodes() {
            return leafNodes;
        }

    }

    private PolygonsSet buildSet(Vector2D[][] vertices) {
        ArrayList<SubLine> edges = new ArrayList<>();
        for (final Vector2D[] vertex : vertices) {
            int l = vertex.length;
            for (int j = 0; j < l; ++j) {
                edges.add(buildSegment(vertex[j], vertex[(j + 1) % l]));
            }
        }
        return new PolygonsSet(edges, 1.0e-10);
    }

    private SubLine buildLine(Vector2D start, Vector2D end) {
        return new Line(start, end, 1.0e-10).wholeHyperplane();
    }

    private double intersectionAbscissa(Line l0, Line l1) {
        Vector2D p = l0.intersection(l1);
        return (l0.toSubSpace(p)).getX();
    }

    private SubLine buildHalfLine(Vector2D start, Vector2D end, boolean startIsVirtual) {
        Line   line  = new Line(start, end, 1.0e-10);
        double lower = startIsVirtual ? Double.NEGATIVE_INFINITY : (line.toSubSpace(start)).getX();
        double upper = startIsVirtual ? (line.toSubSpace(end)).getX() : Double.POSITIVE_INFINITY;
        return new SubLine(line, new IntervalsSet(lower, upper, 1.0e-10));
    }

    private SubLine buildSegment(Vector2D start, Vector2D end) {
        Line   line  = new Line(start, end, 1.0e-10);
        double lower = (line.toSubSpace(start)).getX();
        double upper = (line.toSubSpace(end)).getX();
        return new SubLine(line, new IntervalsSet(lower, upper, 1.0e-10));
    }

    private void checkPoints(Region.Location expected, PolygonsSet set, Vector2D[] points) {
        for (final Vector2D point : points) {
            assertEquals(expected, set.checkPoint(point));
        }
    }

    private boolean checkInSegment(Vector2D p,
                                   Vector2D p1, Vector2D p2,
                                   double tolerance) {
        Line line = new Line(p1, p2, 1.0e-10);
        if (line.getOffset(p) < tolerance) {
            double x  = (line.toSubSpace(p)).getX();
            double x1 = (line.toSubSpace(p1)).getX();
            double x2 = (line.toSubSpace(p2)).getX();
            return (((x - x1) * (x - x2) <= 0.0)
                    || (p1.distance(p) < tolerance)
                    || (p2.distance(p) < tolerance));
        } else {
            return false;
        }
    }

    private void checkVertices(Vector2D[][] rebuiltVertices,
                               Vector2D[][] vertices) {

        // each rebuilt vertex should be in a segment joining two original vertices
        for (final Vector2D[] rebuiltVertex : rebuiltVertices) {
            for (final Vector2D vertex : rebuiltVertex) {
                boolean  inSegment = false;
                for (Vector2D[] loop : vertices) {
                    int length = loop.length;
                    for (int l = 0; (!inSegment) && (l < length); ++l) {
                        inSegment = checkInSegment(vertex, loop[l], loop[(l + 1) % length], 1.0e-10);
                    }
                }
                assertTrue(inSegment);
            }
        }

        // each original vertex should have a corresponding rebuilt vertex
        for (final Vector2D[] vertex : vertices) {
            for (final Vector2D vector2D : vertex) {
                double min = Double.POSITIVE_INFINITY;
                for (final Vector2D[] rebuiltVertex : rebuiltVertices) {
                    for (final Vector2D d : rebuiltVertex) {
                        min = FastMath.min(vector2D.distance(d), min);
                    }
                }
                assertEquals(0.0, min, 1.0e-10);
            }
        }

    }

}
