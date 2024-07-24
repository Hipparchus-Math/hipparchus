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
package org.hipparchus.linear;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.fraction.Fraction;
import org.hipparchus.fraction.FractionField;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;


/**
 * Test cases for the {@link SparseFieldVector} class.
 *
 */
class SparseFieldVectorTest {

    //
    protected Fraction[][] ma1 = {{new Fraction(1), new Fraction(2), new Fraction(3)}, {new Fraction(4), new Fraction(5), new Fraction(6)}, {new Fraction(7), new Fraction(8), new Fraction(9)}};
    protected Fraction[] vec1 = {new Fraction(1), new Fraction(2), new Fraction(3)};
    protected Fraction[] vec2 = {new Fraction(4), new Fraction(5), new Fraction(6)};
    protected Fraction[] vec3 = {new Fraction(7), new Fraction(8), new Fraction(9)};
    protected Fraction[] vec4 = {new Fraction(1), new Fraction(2), new Fraction(3), new Fraction(4), new Fraction(5), new Fraction(6), new Fraction(7), new Fraction(8), new Fraction(9)};
    protected Fraction[] vec_null = {new Fraction(0), new Fraction(0), new Fraction(0)};
    protected Fraction[] dvec1 = {new Fraction(1), new Fraction(2), new Fraction(3), new Fraction(4), new Fraction(5), new Fraction(6), new Fraction(7), new Fraction(8),new Fraction(9)};
    protected Fraction[][] mat1 = {{new Fraction(1), new Fraction(2), new Fraction(3)}, {new Fraction(4), new Fraction(5), new Fraction(6)},{ new Fraction(7), new Fraction(8), new Fraction(9)}};

    // tolerances
    protected double entryTolerance = 10E-16;
    protected double normTolerance = 10E-14;

    protected FractionField field = FractionField.getInstance();

    @Test
    void testMapFunctions() {
        SparseFieldVector<Fraction> v1 = new SparseFieldVector<Fraction>(field,vec1);

        //octave =  v1 .+ 2.0
        FieldVector<Fraction> v_mapAdd = v1.mapAdd(new Fraction(2));
        Fraction[] result_mapAdd = {new Fraction(3), new Fraction(4), new Fraction(5)};
        assertArrayEquals(result_mapAdd,v_mapAdd.toArray(),"compare vectors");

        //octave =  v1 .+ 2.0
        FieldVector<Fraction> v_mapAddToSelf = v1.copy();
        v_mapAddToSelf.mapAddToSelf(new Fraction(2));
        Fraction[] result_mapAddToSelf = {new Fraction(3), new Fraction(4), new Fraction(5)};
        assertArrayEquals(result_mapAddToSelf,v_mapAddToSelf.toArray(),"compare vectors");

        //octave =  v1 .- 2.0
        FieldVector<Fraction> v_mapSubtract = v1.mapSubtract(new Fraction(2));
        Fraction[] result_mapSubtract = {new Fraction(-1), new Fraction(0), new Fraction(1)};
        assertArrayEquals(result_mapSubtract,v_mapSubtract.toArray(),"compare vectors");

        //octave =  v1 .- 2.0
        FieldVector<Fraction> v_mapSubtractToSelf = v1.copy();
        v_mapSubtractToSelf.mapSubtractToSelf(new Fraction(2));
        Fraction[] result_mapSubtractToSelf = {new Fraction(-1), new Fraction(0), new Fraction(1)};
        assertArrayEquals(result_mapSubtractToSelf,v_mapSubtractToSelf.toArray(),"compare vectors");

        //octave =  v1 .* 2.0
        FieldVector<Fraction> v_mapMultiply = v1.mapMultiply(new Fraction(2));
        Fraction[] result_mapMultiply = {new Fraction(2), new Fraction(4), new Fraction(6)};
        assertArrayEquals(result_mapMultiply,v_mapMultiply.toArray(),"compare vectors");

        //octave =  v1 .* 2.0
        FieldVector<Fraction> v_mapMultiplyToSelf = v1.copy();
        v_mapMultiplyToSelf.mapMultiplyToSelf(new Fraction(2));
        Fraction[] result_mapMultiplyToSelf = {new Fraction(2), new Fraction(4), new Fraction(6)};
        assertArrayEquals(result_mapMultiplyToSelf,v_mapMultiplyToSelf.toArray(),"compare vectors");

        //octave =  v1 ./ 2.0
        FieldVector<Fraction> v_mapDivide = v1.mapDivide(new Fraction(2));
        Fraction[] result_mapDivide = {new Fraction(.5d), new Fraction(1), new Fraction(1.5d)};
        assertArrayEquals(result_mapDivide,v_mapDivide.toArray(),"compare vectors");

        //octave =  v1 ./ 2.0
        FieldVector<Fraction> v_mapDivideToSelf = v1.copy();
        v_mapDivideToSelf.mapDivideToSelf(new Fraction(2));
        Fraction[] result_mapDivideToSelf = {new Fraction(.5d), new Fraction(1), new Fraction(1.5d)};
        assertArrayEquals(result_mapDivideToSelf,v_mapDivideToSelf.toArray(),"compare vectors");

        //octave =  v1 .^-1
        FieldVector<Fraction> v_mapInv = v1.mapInv();
        Fraction[] result_mapInv = {new Fraction(1),new Fraction(0.5d),new Fraction(3.333333333333333e-01d)};
        assertArrayEquals(result_mapInv,v_mapInv.toArray(),"compare vectors");

        //octave =  v1 .^-1
        FieldVector<Fraction> v_mapInvToSelf = v1.copy();
        v_mapInvToSelf.mapInvToSelf();
        Fraction[] result_mapInvToSelf = {new Fraction(1),new Fraction(0.5d),new Fraction(3.333333333333333e-01d)};
        assertArrayEquals(result_mapInvToSelf,v_mapInvToSelf.toArray(),"compare vectors");


    }

    @Test
    void testBasicFunctions() {
        SparseFieldVector<Fraction> v1 = (SparseFieldVector<Fraction>) new SparseFieldVector<>(field).append(new SparseFieldVector<>(field,vec1));
        SparseFieldVector<Fraction> v2 = (SparseFieldVector<Fraction>) new SparseFieldVector<>(field).append(new ArrayFieldVector<>(vec2));

        assertSame(field, v1.getField());
        FieldVector<Fraction> v2_t = new ArrayFieldVectorTest.FieldVectorTestImpl<Fraction>(vec2);

        //octave =  v1 + v2
        FieldVector<Fraction> v_add = v1.add(v2);
        Fraction[] result_add = {new Fraction(5), new Fraction(7), new Fraction(9)};
        assertArrayEquals(v_add.toArray(),result_add,"compare vect");

        FieldVector<Fraction> vt2 = new ArrayFieldVectorTest.FieldVectorTestImpl<Fraction>(vec2);
        FieldVector<Fraction> v_add_i = v1.add(vt2);
        Fraction[] result_add_i = {new Fraction(5), new Fraction(7), new Fraction(9)};
        assertArrayEquals(v_add_i.toArray(),result_add_i,"compare vect");

        //octave =  v1 - v2
        SparseFieldVector<Fraction> v_subtract = v1.subtract(v2);
        Fraction[] result_subtract = {new Fraction(-3), new Fraction(-3), new Fraction(-3)};
        customAssertClose("compare vect" , v_subtract.toArray(), result_subtract, normTolerance);

        FieldVector<Fraction> v_subtract_i = v1.subtract(vt2);
        Fraction[] result_subtract_i = {new Fraction(-3), new Fraction(-3), new Fraction(-3)};
        customAssertClose("compare vect" , v_subtract_i.toArray(), result_subtract_i, normTolerance);

        // octave v1 .* v2
        FieldVector<Fraction>  v_ebeMultiply = v1.ebeMultiply(v2);
        Fraction[] result_ebeMultiply = {new Fraction(4), new Fraction(10), new Fraction(18)};
        customAssertClose("compare vect" , v_ebeMultiply.toArray(), result_ebeMultiply, normTolerance);

        FieldVector<Fraction>  v_ebeMultiply_2 = v1.ebeMultiply(v2_t);
        Fraction[] result_ebeMultiply_2 = {new Fraction(4), new Fraction(10), new Fraction(18)};
        customAssertClose("compare vect" , v_ebeMultiply_2.toArray(), result_ebeMultiply_2, normTolerance);

        // octave v1 ./ v2
        FieldVector<Fraction>  v_ebeDivide = v1.ebeDivide(v2);
        Fraction[] result_ebeDivide = {new Fraction(0.25d), new Fraction(0.4d), new Fraction(0.5d)};
        customAssertClose("compare vect" , v_ebeDivide.toArray(), result_ebeDivide, normTolerance);

        FieldVector<Fraction>  v_ebeDivide_2 = v1.ebeDivide(v2_t);
        Fraction[] result_ebeDivide_2 = {new Fraction(0.25d), new Fraction(0.4d), new Fraction(0.5d)};
        customAssertClose("compare vect" , v_ebeDivide_2.toArray(), result_ebeDivide_2, normTolerance);

        // octave  dot(v1,v2)
        Fraction dot =  v1.dotProduct(v2);
        assertEquals(new Fraction(32), dot, "compare val ");

        // octave  dot(v1,v2_t)
        Fraction dot_2 =  v1.dotProduct(v2_t);
        assertEquals(new Fraction(32), dot_2, "compare val ");

        FieldMatrix<Fraction> m_outerProduct = v1.outerProduct(v2);
        assertEquals(new Fraction(4), m_outerProduct.getEntry(0,0), "compare val ");

        FieldMatrix<Fraction> m_outerProduct_2 = v1.outerProduct(v2_t);
        assertEquals(new Fraction(4), m_outerProduct_2.getEntry(0,0), "compare val ");

    }

    @Test
    void testOuterProduct() {
        final SparseFieldVector<Fraction> u
            = (SparseFieldVector<Fraction>) new SparseFieldVector<>(FractionField.getInstance(),
                                                                    new Fraction[] { new Fraction(1),
                                                                                     new Fraction(2) }).
                                            append(new Fraction(-3));
        final SparseFieldVector<Fraction> v
            = new SparseFieldVector<Fraction>(FractionField.getInstance(),
                                              new Fraction[] {new Fraction(4),
                                                              new Fraction(-2)});

        final FieldMatrix<Fraction> uv = u.outerProduct(v);

        final double tol = Math.ulp(1d);
        assertEquals(new Fraction(4).doubleValue(), uv.getEntry(0, 0).doubleValue(), tol);
        assertEquals(new Fraction(-2).doubleValue(), uv.getEntry(0, 1).doubleValue(), tol);
        assertEquals(new Fraction(8).doubleValue(), uv.getEntry(1, 0).doubleValue(), tol);
        assertEquals(new Fraction(-4).doubleValue(), uv.getEntry(1, 1).doubleValue(), tol);
        assertEquals(new Fraction(-12).doubleValue(), uv.getEntry(2, 0).doubleValue(), tol);
        assertEquals(new Fraction(6).doubleValue(), uv.getEntry(2, 1).doubleValue(), tol);
    }

    @Test
    void testMisc() {
        SparseFieldVector<Fraction> v1 = new SparseFieldVector<Fraction>(field,vec1);

        String out1 = v1.toString();
        assertTrue(out1.length()!=0,  "some output ");
        try {
            v1.checkVectorDimensions(2);
            fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected behavior
        }


    }

    @Test
    void testPredicates() {

        SparseFieldVector<Fraction> v = new SparseFieldVector<Fraction>(field, new Fraction[] { new Fraction(0), new Fraction(1), new Fraction(2) });

        v.setEntry(0, field.getZero());
        assertEquals(v, new SparseFieldVector<Fraction>(field, new Fraction[] { new Fraction(0), new Fraction(1), new Fraction(2) }));
        assertNotSame(v, new SparseFieldVector<Fraction>(field, new Fraction[] { new Fraction(0), new Fraction(1), new Fraction(2), new Fraction(3) }));

    }

    /** verifies that two vectors are close (sup norm) */
    protected void customAssertEquals(String msg, Fraction[] m, Fraction[] n) {
        if (m.length != n.length) {
            fail("vectors have different lengths");
        }
        for (int i = 0; i < m.length; i++) {
            assertEquals(m[i],n[i],msg + " " +  i + " elements differ");
        }
    }

    /** verifies that two vectors are close (sup norm) */
    protected void customAssertClose(String msg, Fraction[] m, Fraction[] n, double tolerance) {
        if (m.length != n.length) {
            fail("vectors have different lengths");
        }
        for (int i = 0; i < m.length; i++) {
            assertEquals(m[i].doubleValue(),n[i].doubleValue(), tolerance, msg + " " +  i + " elements differ");
        }
    }

    /*
     * TESTS OF THE VISITOR PATTERN
     */

    /** The whole vector is visited. */
    @Test
    void testWalkInDefaultOrderPreservingVisitor1() {
        final Fraction[] data = new Fraction[] {
            Fraction.ZERO, Fraction.ONE, Fraction.ZERO,
            Fraction.ZERO, Fraction.TWO, Fraction.ZERO,
            Fraction.ZERO, Fraction.ZERO, new Fraction(3)
        };
        final SparseFieldVector<Fraction> v = new SparseFieldVector<Fraction>(field, data);
        final FieldVectorPreservingVisitor<Fraction> visitor;
        visitor = new FieldVectorPreservingVisitor<Fraction>() {

            private int expectedIndex;

            public void visit(final int actualIndex, final Fraction actualValue) {
                assertEquals(expectedIndex, actualIndex);
                assertEquals(data[actualIndex], actualValue, Integer.toString(actualIndex));
                ++expectedIndex;
            }

            public void start(final int actualSize, final int actualStart,
                              final int actualEnd) {
                assertEquals(data.length, actualSize);
                assertEquals(0, actualStart);
                assertEquals(data.length - 1, actualEnd);
                expectedIndex = 0;
            }

            public Fraction end() {
                return Fraction.ZERO;
            }
        };
        v.walkInDefaultOrder(visitor);
    }

    /** Visiting an invalid subvector. */
    @Test
    void testWalkInDefaultOrderPreservingVisitor2() {
        final SparseFieldVector<Fraction> v = create(5);
        final FieldVectorPreservingVisitor<Fraction> visitor;
        visitor = new FieldVectorPreservingVisitor<Fraction>() {

            public void visit(int index, Fraction value) {
                // Do nothing
            }

            public void start(int dimension, int start, int end) {
                // Do nothing
            }

            public Fraction end() {
                return Fraction.ZERO;
            }
        };
        try {
            v.walkInDefaultOrder(visitor, -1, 4);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 5, 4);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 0, -1);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 0, 5);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 4, 0);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
    }

    /** Visiting a valid subvector. */
    @Test
    void testWalkInDefaultOrderPreservingVisitor3() {
        final Fraction[] data = new Fraction[] {
            Fraction.ZERO, Fraction.ONE, Fraction.ZERO,
            Fraction.ZERO, Fraction.TWO, Fraction.ZERO,
            Fraction.ZERO, Fraction.ZERO, new Fraction(3)
        };
        final SparseFieldVector<Fraction> v = new SparseFieldVector<Fraction>(field, data);
        final int expectedStart = 2;
        final int expectedEnd = 7;
        final FieldVectorPreservingVisitor<Fraction> visitor;
        visitor = new FieldVectorPreservingVisitor<Fraction>() {

            private int expectedIndex;

            public void visit(final int actualIndex, final Fraction actualValue) {
                assertEquals(expectedIndex, actualIndex);
                assertEquals(data[actualIndex], actualValue, Integer.toString(actualIndex));
                ++expectedIndex;
            }

            public void start(final int actualSize, final int actualStart,
                              final int actualEnd) {
                assertEquals(data.length, actualSize);
                assertEquals(expectedStart, actualStart);
                assertEquals(expectedEnd, actualEnd);
                expectedIndex = expectedStart;
            }

            public Fraction end() {
                return Fraction.ZERO;
            }
        };
        v.walkInDefaultOrder(visitor, expectedStart, expectedEnd);
    }

    /** The whole vector is visited. */
    @Test
    void testWalkInOptimizedOrderPreservingVisitor1() {
        final Fraction[] data = new Fraction[] {
            Fraction.ZERO, Fraction.ONE, Fraction.ZERO,
            Fraction.ZERO, Fraction.TWO, Fraction.ZERO,
            Fraction.ZERO, Fraction.ZERO, new Fraction(3)
        };
        final SparseFieldVector<Fraction> v = new SparseFieldVector<Fraction>(field, data);
        final FieldVectorPreservingVisitor<Fraction> visitor;
        visitor = new FieldVectorPreservingVisitor<Fraction>() {
            private final boolean[] visited = new boolean[data.length];

            public void visit(final int actualIndex, final Fraction actualValue) {
                visited[actualIndex] = true;
                assertEquals(data[actualIndex], actualValue, Integer.toString(actualIndex));
            }

            public void start(final int actualSize, final int actualStart,
                              final int actualEnd) {
                assertEquals(data.length, actualSize);
                assertEquals(0, actualStart);
                assertEquals(data.length - 1, actualEnd);
                Arrays.fill(visited, false);
            }

            public Fraction end() {
                for (int i = 0; i < data.length; i++) {
                    assertTrue(visited[i],
                                      "entry " + i + "has not been visited");
                }
                return Fraction.ZERO;
            }
        };
        v.walkInOptimizedOrder(visitor);
    }

    /** Visiting an invalid subvector. */
    @Test
    void testWalkInOptimizedOrderPreservingVisitor2() {
        final SparseFieldVector<Fraction> v = create(5);
        final FieldVectorPreservingVisitor<Fraction> visitor;
        visitor = new FieldVectorPreservingVisitor<Fraction>() {

            public void visit(int index, Fraction value) {
                // Do nothing
            }

            public void start(int dimension, int start, int end) {
                // Do nothing
            }

            public Fraction end() {
                return Fraction.ZERO;
            }
        };
        try {
            v.walkInOptimizedOrder(visitor, -1, 4);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 5, 4);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 0, -1);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 0, 5);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 4, 0);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
    }

    /** Visiting a valid subvector. */
    @Test
    void testWalkInOptimizedOrderPreservingVisitor3() {
        final Fraction[] data = new Fraction[] {
            Fraction.ZERO, Fraction.ONE, Fraction.ZERO,
            Fraction.ZERO, Fraction.TWO, Fraction.ZERO,
            Fraction.ZERO, Fraction.ZERO, new Fraction(3)
        };
        final SparseFieldVector<Fraction> v = new SparseFieldVector<Fraction>(field, data);
        final int expectedStart = 2;
        final int expectedEnd = 7;
        final FieldVectorPreservingVisitor<Fraction> visitor;
        visitor = new FieldVectorPreservingVisitor<Fraction>() {
            private final boolean[] visited = new boolean[data.length];

            public void visit(final int actualIndex, final Fraction actualValue) {
                assertEquals(data[actualIndex], actualValue, Integer.toString(actualIndex));
                visited[actualIndex] = true;
            }

            public void start(final int actualSize, final int actualStart,
                              final int actualEnd) {
                assertEquals(data.length, actualSize);
                assertEquals(expectedStart, actualStart);
                assertEquals(expectedEnd, actualEnd);
                Arrays.fill(visited, true);
            }

            public Fraction end() {
                for (int i = expectedStart; i <= expectedEnd; i++) {
                    assertTrue(visited[i],
                                      "entry " + i + "has not been visited");
                }
                return Fraction.ZERO;
            }
        };
        v.walkInOptimizedOrder(visitor, expectedStart, expectedEnd);
    }

    /** The whole vector is visited. */
    @Test
    void testWalkInDefaultOrderChangingVisitor1() {
        final Fraction[] data = new Fraction[] {
            Fraction.ZERO, Fraction.ONE, Fraction.ZERO,
            Fraction.ZERO, Fraction.TWO, Fraction.ZERO,
            Fraction.ZERO, Fraction.ZERO, new Fraction(3)
        };
        final SparseFieldVector<Fraction> v = new SparseFieldVector<Fraction>(field, data);
        final FieldVectorChangingVisitor<Fraction> visitor;
        visitor = new FieldVectorChangingVisitor<Fraction>() {

            private int expectedIndex;

            public Fraction visit(final int actualIndex, final Fraction actualValue) {
                assertEquals(expectedIndex, actualIndex);
                assertEquals(data[actualIndex], actualValue, Integer.toString(actualIndex));
                ++expectedIndex;
                return actualValue.add(actualIndex);
            }

            public void start(final int actualSize, final int actualStart,
                              final int actualEnd) {
                assertEquals(data.length, actualSize);
                assertEquals(0, actualStart);
                assertEquals(data.length - 1, actualEnd);
                expectedIndex = 0;
            }

            public Fraction end() {
                return Fraction.ZERO;
            }
        };
        v.walkInDefaultOrder(visitor);
        for (int i = 0; i < data.length; i++) {
            assertEquals(data[i].add(i), v.getEntry(i), "entry " + i);
        }
    }

    /** Visiting an invalid subvector. */
    @Test
    void testWalkInDefaultOrderChangingVisitor2() {
        final SparseFieldVector<Fraction> v = create(5);
        final FieldVectorChangingVisitor<Fraction> visitor;
        visitor = new FieldVectorChangingVisitor<Fraction>() {

            public Fraction visit(int index, Fraction value) {
                return Fraction.ZERO;
            }

            public void start(int dimension, int start, int end) {
                // Do nothing
            }

            public Fraction end() {
                return Fraction.ZERO;
            }
        };
        try {
            v.walkInDefaultOrder(visitor, -1, 4);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 5, 4);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 0, -1);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 0, 5);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 4, 0);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
    }

    /** Visiting a valid subvector. */
    @Test
    void testWalkInDefaultOrderChangingVisitor3() {
        final Fraction[] data = new Fraction[] {
            Fraction.ZERO, Fraction.ONE, Fraction.ZERO,
            Fraction.ZERO, Fraction.TWO, Fraction.ZERO,
            Fraction.ZERO, Fraction.ZERO, new Fraction(3)
        };
        final SparseFieldVector<Fraction> v = new SparseFieldVector<Fraction>(field, data);
        final int expectedStart = 2;
        final int expectedEnd = 7;
        final FieldVectorChangingVisitor<Fraction> visitor;
        visitor = new FieldVectorChangingVisitor<Fraction>() {

            private int expectedIndex;

            public Fraction visit(final int actualIndex, final Fraction actualValue) {
                assertEquals(expectedIndex, actualIndex);
                assertEquals(data[actualIndex], actualValue, Integer.toString(actualIndex));
                ++expectedIndex;
                return actualValue.add(actualIndex);
            }

            public void start(final int actualSize, final int actualStart,
                              final int actualEnd) {
                assertEquals(data.length, actualSize);
                assertEquals(expectedStart, actualStart);
                assertEquals(expectedEnd, actualEnd);
                expectedIndex = expectedStart;
            }

            public Fraction end() {
                return Fraction.ZERO;
            }
        };
        v.walkInDefaultOrder(visitor, expectedStart, expectedEnd);
        for (int i = expectedStart; i <= expectedEnd; i++) {
            assertEquals(data[i].add(i), v.getEntry(i), "entry " + i);
        }
    }

    /** The whole vector is visited. */
    @Test
    void testWalkInOptimizedOrderChangingVisitor1() {
        final Fraction[] data = new Fraction[] {
            Fraction.ZERO, Fraction.ONE, Fraction.ZERO,
            Fraction.ZERO, Fraction.TWO, Fraction.ZERO,
            Fraction.ZERO, Fraction.ZERO, new Fraction(3)
        };
        final SparseFieldVector<Fraction> v = new SparseFieldVector<Fraction>(field, data);
        final FieldVectorChangingVisitor<Fraction> visitor;
        visitor = new FieldVectorChangingVisitor<Fraction>() {
            private final boolean[] visited = new boolean[data.length];

            public Fraction visit(final int actualIndex, final Fraction actualValue) {
                visited[actualIndex] = true;
                assertEquals(data[actualIndex], actualValue, Integer.toString(actualIndex));
                return actualValue.add(actualIndex);
            }

            public void start(final int actualSize, final int actualStart,
                              final int actualEnd) {
                assertEquals(data.length, actualSize);
                assertEquals(0, actualStart);
                assertEquals(data.length - 1, actualEnd);
                Arrays.fill(visited, false);
            }

            public Fraction end() {
                for (int i = 0; i < data.length; i++) {
                    assertTrue(visited[i],
                                      "entry " + i + "has not been visited");
                }
                return Fraction.ZERO;
            }
        };
        v.walkInOptimizedOrder(visitor);
        for (int i = 0; i < data.length; i++) {
            assertEquals(data[i].add(i), v.getEntry(i), "entry " + i);
        }
    }

    /** Visiting an invalid subvector. */
    @Test
    void testWalkInOptimizedOrderChangingVisitor2() {
        final SparseFieldVector<Fraction> v = create(5);
        final FieldVectorChangingVisitor<Fraction> visitor;
        visitor = new FieldVectorChangingVisitor<Fraction>() {

            public Fraction visit(int index, Fraction value) {
                return Fraction.ZERO;
            }

            public void start(int dimension, int start, int end) {
                // Do nothing
            }

            public Fraction end() {
                return Fraction.ZERO;
            }
        };
        try {
            v.walkInOptimizedOrder(visitor, -1, 4);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 5, 4);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 0, -1);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 0, 5);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 4, 0);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
    }

    /** Visiting a valid subvector. */
    @Test
    void testWalkInOptimizedOrderChangingVisitor3() {
        final Fraction[] data = new Fraction[] {
            Fraction.ZERO, Fraction.ONE, Fraction.ZERO,
            Fraction.ZERO, Fraction.TWO, Fraction.ZERO,
            Fraction.ZERO, Fraction.ZERO, new Fraction(3)
        };
        final SparseFieldVector<Fraction> v = new SparseFieldVector<Fraction>(field, data);
        final int expectedStart = 2;
        final int expectedEnd = 7;
        final FieldVectorChangingVisitor<Fraction> visitor;
        visitor = new FieldVectorChangingVisitor<Fraction>() {
            private final boolean[] visited = new boolean[data.length];

            public Fraction visit(final int actualIndex, final Fraction actualValue) {
                assertEquals(data[actualIndex], actualValue, Integer.toString(actualIndex));
                visited[actualIndex] = true;
                return actualValue.add(actualIndex);
            }

            public void start(final int actualSize, final int actualStart,
                              final int actualEnd) {
                assertEquals(data.length, actualSize);
                assertEquals(expectedStart, actualStart);
                assertEquals(expectedEnd, actualEnd);
                Arrays.fill(visited, true);
            }

            public Fraction end() {
                for (int i = expectedStart; i <= expectedEnd; i++) {
                    assertTrue(visited[i],
                                      "entry " + i + "has not been visited");
                }
                return Fraction.ZERO;
            }
        };
        v.walkInOptimizedOrder(visitor, expectedStart, expectedEnd);
        for (int i = expectedStart; i <= expectedEnd; i++) {
            assertEquals(data[i].add(i), v.getEntry(i), "entry " + i);
        }
    }

    private SparseFieldVector<Fraction> create(int n) {
        Fraction[] t = new Fraction[n];
        for (int i = 0; i < n; ++i) {
            t[i] = Fraction.ZERO;
        }
        return new SparseFieldVector<Fraction>(field, t);
    }
}
