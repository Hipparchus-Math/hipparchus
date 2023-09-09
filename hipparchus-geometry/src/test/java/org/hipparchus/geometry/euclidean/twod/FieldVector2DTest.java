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
package org.hipparchus.geometry.euclidean.twod;

import java.text.NumberFormat;
import java.util.Locale;

import org.hipparchus.Field;
import org.hipparchus.CalculusFieldElement;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.geometry.LocalizedGeometryFormats;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;
import org.hipparchus.util.SinCos;
import org.junit.Assert;
import org.junit.Test;

public class FieldVector2DTest {

    @Test
    public void testConstructors() {
        doTestConstructors(Binary64Field.getInstance());
    }

    @Test
    public void testConstants() {
        doTestConstants(Binary64Field.getInstance());
    }

    @Test
    public void testToMethods() {
        doTestToMethods(Binary64Field.getInstance());
    }

    @Test
    public void testNorms() {
        doTestNorms(Binary64Field.getInstance());
    }

    @Test
    public void testDistances() {
        doTestDistances(Binary64Field.getInstance());
    }

    @Test
    public void testAdd() {
        doTestAdd(Binary64Field.getInstance());
    }

    @Test
    public void testSubtract() {
        doTestSubtract(Binary64Field.getInstance());
    }

    @Test
    public void testNormalize() {
        doTestNormalize(Binary64Field.getInstance());
    }

    @Test
    public void testAngle() {
        doTestAngle(Binary64Field.getInstance());
    }

    @Test
    public void testNegate() {
        doTestNegate(Binary64Field.getInstance());
    }

    @Test
    public void testScalarMultiply() {
        doTestScalarMultiply(Binary64Field.getInstance());
    }

    @Test
    public void testIsNaN() {
        doTestIsNaN(Binary64Field.getInstance());
    }

    @Test
    public void testIsInfinite() {
        doTestIsInfinite(Binary64Field.getInstance());
    }

    @Test
    public void testEquals() {
        doTestEquals(Binary64Field.getInstance());
    }

    @Test
    public void testHashCode() {
        doTestHashCode(Binary64Field.getInstance());
    }

    @Test
    public void testCrossProduct() {
        doTestCrossProduct(Binary64Field.getInstance());
    }

    @Test
    public void testOrientation() {
        doTestOrientation(Binary64Field.getInstance());
    }

    private <T extends CalculusFieldElement<T>> void doTestConstructors(final Field<T> field) {
        final T p40 = field.getZero().add( 4.0);
        final T p20 = field.getZero().add( 2.0);
        final T p25 = field.getZero().add( 2.5);
        final T p10 = field.getOne();
        final T m05 = field.getZero().add(-0.5);
        final T m30 = field.getZero().add(-3.0);
        check(new FieldVector2D<>(p25, m05), 2.5, -0.5, 1.0e-15);
        final T[] a = MathArrays.buildArray(field, 2);
        a[0] = field.getOne();
        a[1] = field.getZero();
        check(new FieldVector2D<>(a), 1.0, 0.0, 1.0e-15);
        try {
            new FieldVector2D<>(MathArrays.buildArray(field, 3));
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            Assert.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
            Assert.assertEquals(3, ((Integer) miae.getParts()[0]).intValue());
            Assert.assertEquals(2, ((Integer) miae.getParts()[1]).intValue());
        }
        check(new FieldVector2D<>(p20, new FieldVector2D<>(p25, m05)), 5.0, -1.0, 1.0e-15);
        check(new FieldVector2D<>(p20, new Vector2D(2.5, -0.5)), 5.0, -1.0, 1.0e-15);
        check(new FieldVector2D<>(2, new FieldVector2D<>(p25, m05)), 5.0, -1.0, 1.0e-15);
        check(new FieldVector2D<>(p20, new FieldVector2D<>(p25, m05), m30, new FieldVector2D<>(m05, p40)),
              6.5, -13.0, 1.0e-15);
        check(new FieldVector2D<>(p20, new Vector2D(2.5, -0.5), m30, new Vector2D(-0.5, 4.0)),
              6.5, -13.0, 1.0e-15);
        check(new FieldVector2D<>(2.0, new FieldVector2D<>(p25, m05), -3.0, new FieldVector2D<>(m05, p40)),
              6.5, -13.0, 1.0e-15);
        check(new FieldVector2D<>(p20, new FieldVector2D<>(p25, m05), m30, new FieldVector2D<>(m05, p40),
                                  p40, new FieldVector2D<>(p25, m30)),
              16.5, -25.0, 1.0e-15);
        check(new FieldVector2D<>(p20, new Vector2D(2.5, -0.5), m30, new Vector2D(-0.5, 4.0),
                                  p40, new Vector2D(2.5, -3.0)),
              16.5, -25.0, 1.0e-15);
        check(new FieldVector2D<>(2.0, new FieldVector2D<>(p25, m05), -3.0, new FieldVector2D<>(m05, p40),
                                  4.0, new FieldVector2D<>(p25, m30)),
              16.5, -25.0, 1.0e-15);
        check(new FieldVector2D<>(p20, new FieldVector2D<>(p25, m05), m30, new FieldVector2D<>(m05, p40),
                                  p40, new FieldVector2D<>(p25, m30), p10, new FieldVector2D<>(p10, p10)),
              17.5, -24.0, 1.0e-15);
        check(new FieldVector2D<>(p20, new Vector2D(2.5, -0.5), m30, new Vector2D(-0.5, 4.0),
                                  p40, new Vector2D(2.5, -3.0), p10, new Vector2D(1.0, 1.0)),
              17.5, -24.0, 1.0e-15);
        check(new FieldVector2D<>(2.0, new FieldVector2D<>(p25, m05), -3.0, new FieldVector2D<>(m05, p40),
                                  4.0, new FieldVector2D<>(p25, m30),  1.0, new FieldVector2D<>(p10, p10)),
              17.5, -24.0, 1.0e-15);
    }

    private <T extends CalculusFieldElement<T>> void doTestConstants(final Field<T> field) {
        check(FieldVector2D.getZero(field),    0.0,  0.0, 1.0e-15);
        check(FieldVector2D.getPlusI(field),   1.0,  0.0, 1.0e-15);
        check(FieldVector2D.getMinusI(field), -1.0,  0.0, 1.0e-15);
        check(FieldVector2D.getPlusJ(field),   0.0,  1.0, 1.0e-15);
        check(FieldVector2D.getMinusJ(field),  0.0, -1.0, 1.0e-15);
        check(FieldVector2D.getPositiveInfinity(field), Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, 1.0e-15);
        check(FieldVector2D.getNegativeInfinity(field), Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, 1.0e-15);
        Assert.assertTrue(Double.isNaN(FieldVector2D.getNaN(field).getX().getReal()));
        Assert.assertTrue(Double.isNaN(FieldVector2D.getNaN(field).getY().getReal()));
    }

    private <T extends CalculusFieldElement<T>> void doTestToMethods(final Field<T> field) {
        final FieldVector2D<T> v = new FieldVector2D<>(field, new Vector2D(2.5, -0.5));
        Assert.assertEquals( 2,   v.toArray().length);
        Assert.assertEquals( 2.5, v.toArray()[0].getReal(), 1.0e-15);
        Assert.assertEquals(-0.5, v.toArray()[1].getReal(), 1.0e-15);
        Assert.assertEquals(new Vector2D(2.5, -0.5), v.toVector2D());
        Assert.assertEquals("{2.5; -0.5}", v.toString().replaceAll(",", "."));
        Assert.assertEquals("{2,5; -0,5}", v.toString(NumberFormat.getInstance(Locale.FRENCH)));
    }

    private <T extends CalculusFieldElement<T>> void doTestNorms(final Field<T> field) {
        final FieldVector2D<T> v = new FieldVector2D<>(field, new Vector2D(3.0, -4.0));
        Assert.assertEquals( 7.0, v.getNorm1().getReal(),   1.0e-15);
        Assert.assertEquals( 5.0, v.getNorm().getReal(),    1.0e-15);
        Assert.assertEquals(25.0, v.getNormSq().getReal(),  1.0e-15);
        Assert.assertEquals( 4.0, v.getNormInf().getReal(), 1.0e-15);
    }

    private <T extends CalculusFieldElement<T>> void doTestDistances(final Field<T> field) {
        final FieldVector2D<T> u = new FieldVector2D<>(field, new Vector2D( 2.0, -2.0));
        final FieldVector2D<T> v = new FieldVector2D<>(field, new Vector2D(-1.0,  2.0));
        Assert.assertEquals( 7.0, FieldVector2D.distance1(u, v).getReal(),   1.0e-15);
        Assert.assertEquals( 5.0, FieldVector2D.distance(u, v).getReal(),    1.0e-15);
        Assert.assertEquals(25.0, FieldVector2D.distanceSq(u, v).getReal(),  1.0e-15);
        Assert.assertEquals( 4.0, FieldVector2D.distanceInf(u, v).getReal(), 1.0e-15);
        Assert.assertEquals( 7.0, FieldVector2D.distance1(u, v.toVector2D()).getReal(),   1.0e-15);
        Assert.assertEquals( 5.0, FieldVector2D.distance(u, v.toVector2D()).getReal(),    1.0e-15);
        Assert.assertEquals(25.0, FieldVector2D.distanceSq(u, v.toVector2D()).getReal(),  1.0e-15);
        Assert.assertEquals( 4.0, FieldVector2D.distanceInf(u, v.toVector2D()).getReal(), 1.0e-15);
        Assert.assertEquals( 7.0, FieldVector2D.distance1(u.toVector2D(), v).getReal(),   1.0e-15);
        Assert.assertEquals( 5.0, FieldVector2D.distance(u.toVector2D(), v).getReal(),    1.0e-15);
        Assert.assertEquals(25.0, FieldVector2D.distanceSq(u.toVector2D(), v).getReal(),  1.0e-15);
        Assert.assertEquals( 4.0, FieldVector2D.distanceInf(u.toVector2D(), v).getReal(), 1.0e-15);
    }

    private <T extends CalculusFieldElement<T>> void doTestAdd(final Field<T> field) {
        final FieldVector2D<T> u = new FieldVector2D<>(field, new Vector2D( 2.0, -2.0));
        final FieldVector2D<T> v = new FieldVector2D<>(field, new Vector2D(-1.0,  2.0));
        check(u.add(v), 1.0, 0.0, 1.0e-15);
        check(u.add(v.toVector2D()), 1.0, 0.0, 1.0e-15);
        check(u.add(field.getZero().add(5), v), -3.0, 8.0, 1.0e-15);
        check(u.add(field.getZero().add(5), v.toVector2D()), -3.0, 8.0, 1.0e-15);
        check(u.add(5.0, v), -3.0, 8.0, 1.0e-15);
        check(u.add(5.0, v.toVector2D()), -3.0, 8.0, 1.0e-15);
    }

    private <T extends CalculusFieldElement<T>> void doTestSubtract(final Field<T> field) {
        final FieldVector2D<T> u = new FieldVector2D<>(field, new Vector2D( 2.0, -2.0));
        final FieldVector2D<T> v = new FieldVector2D<>(field, new Vector2D( 1.0, -2.0));
        check(u.subtract(v), 1.0, 0.0, 1.0e-15);
        check(u.subtract(v.toVector2D()), 1.0, 0.0, 1.0e-15);
        check(u.subtract(field.getZero().add(5), v), -3.0, 8.0, 1.0e-15);
        check(u.subtract(field.getZero().add(5), v.toVector2D()), -3.0, 8.0, 1.0e-15);
        check(u.subtract(5.0, v), -3.0, 8.0, 1.0e-15);
        check(u.subtract(5.0, v.toVector2D()), -3.0, 8.0, 1.0e-15);
    }

    private <T extends CalculusFieldElement<T>> void doTestNormalize(final Field<T> field) {
        try {
            FieldVector2D.getZero(field).normalize();
            Assert.fail("an exception should habe been thrown");
        } catch (MathRuntimeException mre) {
            Assert.assertEquals(LocalizedGeometryFormats.CANNOT_NORMALIZE_A_ZERO_NORM_VECTOR, mre.getSpecifier());
        }
        check(new FieldVector2D<>(field, new Vector2D(3, -4)).normalize(), 0.6, -0.8, 1.0e-15);
    }

    private <T extends CalculusFieldElement<T>> void doTestAngle(final Field<T> field) {
        try {
            FieldVector2D.angle(FieldVector2D.getZero(field), FieldVector2D.getPlusI(field));
            Assert.fail("an exception should habe been thrown");
        } catch (MathRuntimeException mre) {
            Assert.assertEquals(LocalizedCoreFormats.ZERO_NORM, mre.getSpecifier());
        }
        final double alpha = 0.01;
        final SinCos sc = FastMath.sinCos(alpha);
        Assert.assertEquals(alpha,
                            FieldVector2D.angle(new FieldVector2D<>(field, new Vector2D(sc.cos(), sc.sin())),
                                                FieldVector2D.getPlusI(field)).getReal(),
                            1.0e-15);
        Assert.assertEquals(FastMath.PI - alpha,
                            FieldVector2D.angle(new FieldVector2D<>(field, new Vector2D(-sc.cos(), sc.sin())),
                                                FieldVector2D.getPlusI(field)).getReal(),
                            1.0e-15);
        Assert.assertEquals(0.5 * FastMath.PI - alpha,
                            FieldVector2D.angle(new FieldVector2D<>(field, new Vector2D(sc.sin(), sc.cos())),
                                                FieldVector2D.getPlusI(field)).getReal(),
                            1.0e-15);
        Assert.assertEquals(0.5 * FastMath.PI + alpha,
                            FieldVector2D.angle(new FieldVector2D<>(field, new Vector2D(-sc.sin(), sc.cos())),
                                                FieldVector2D.getPlusI(field)).getReal(),
                            1.0e-15);
        try {
            FieldVector2D.angle(FieldVector2D.getZero(field), Vector2D.PLUS_I);
            Assert.fail("an exception should habe been thrown");
        } catch (MathRuntimeException mre) {
            Assert.assertEquals(LocalizedCoreFormats.ZERO_NORM, mre.getSpecifier());
        }
        Assert.assertEquals(alpha,
                            FieldVector2D.angle(new FieldVector2D<>(field, new Vector2D(sc.cos(), sc.sin())),
                                                Vector2D.PLUS_I).getReal(),
                            1.0e-15);
        Assert.assertEquals(FastMath.PI - alpha,
                            FieldVector2D.angle(new FieldVector2D<>(field, new Vector2D(-sc.cos(), sc.sin())),
                                                Vector2D.PLUS_I).getReal(),
                            1.0e-15);
        Assert.assertEquals(0.5 * FastMath.PI - alpha,
                            FieldVector2D.angle(new FieldVector2D<>(field, new Vector2D(sc.sin(), sc.cos())),
                                                Vector2D.PLUS_I).getReal(),
                            1.0e-15);
        Assert.assertEquals(0.5 * FastMath.PI + alpha,
                            FieldVector2D.angle(new FieldVector2D<>(field, new Vector2D(-sc.sin(), sc.cos())),
                                                Vector2D.PLUS_I).getReal(),
                            1.0e-15);
        try {
            FieldVector2D.angle(Vector2D.ZERO, FieldVector2D.getPlusI(field));
            Assert.fail("an exception should habe been thrown");
        } catch (MathRuntimeException mre) {
            Assert.assertEquals(LocalizedCoreFormats.ZERO_NORM, mre.getSpecifier());
        }
        Assert.assertEquals(alpha,
                            FieldVector2D.angle(new Vector2D(sc.cos(), sc.sin()),
                                                FieldVector2D.getPlusI(field)).getReal(),
                            1.0e-15);
        Assert.assertEquals(FastMath.PI - alpha,
                            FieldVector2D.angle(new Vector2D(-sc.cos(), sc.sin()),
                                                FieldVector2D.getPlusI(field)).getReal(),
                            1.0e-15);
        Assert.assertEquals(0.5 * FastMath.PI - alpha,
                            FieldVector2D.angle(new Vector2D(sc.sin(), sc.cos()),
                                                FieldVector2D.getPlusI(field)).getReal(),
                            1.0e-15);
        Assert.assertEquals(0.5 * FastMath.PI + alpha,
                            FieldVector2D.angle(new Vector2D(-sc.sin(), sc.cos()),
                                                FieldVector2D.getPlusI(field)).getReal(),
                            1.0e-15);
    }


    private <T extends CalculusFieldElement<T>> void doTestNegate(final Field<T> field) {
        check(new FieldVector2D<>(field, new Vector2D(3.0, -4.0)).negate(), -3.0, 4.0, 1.0e-15);
    }

    private <T extends CalculusFieldElement<T>> void doTestScalarMultiply(final Field<T> field) {
        check(new FieldVector2D<>(field, new Vector2D(3.0, -4.0)).scalarMultiply(2.0), 6.0, -8.0, 1.0e-15);
        check(new FieldVector2D<>(field, new Vector2D(3.0, -4.0)).scalarMultiply(field.getZero().add(2.0)), 6.0, -8.0, 1.0e-15);
    }

    private <T extends CalculusFieldElement<T>> void doTestIsNaN(final Field<T> field) {
        Assert.assertTrue(new FieldVector2D<>(field, new Vector2D(Double.NaN, 0.0)).isNaN());
        Assert.assertTrue(new FieldVector2D<>(field, new Vector2D(0.0, Double.NaN)).isNaN());
        Assert.assertTrue(new FieldVector2D<>(field, new Vector2D(Double.NaN, Double.NaN)).isNaN());
        Assert.assertTrue(FieldVector2D.getNaN(field).isNaN());
        Assert.assertFalse(new FieldVector2D<>(field, new Vector2D(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY)).isNaN());
        Assert.assertFalse(FieldVector2D.getMinusI(field).isNaN());
    }

    private <T extends CalculusFieldElement<T>> void doTestIsInfinite(final Field<T> field) {
        Assert.assertFalse(new FieldVector2D<>(field, new Vector2D(Double.NaN, 0.0)).isInfinite());
        Assert.assertTrue(new FieldVector2D<>(field, new Vector2D(Double.POSITIVE_INFINITY, 0.0)).isInfinite());
        Assert.assertTrue(new FieldVector2D<>(field, new Vector2D(0.0, Double.POSITIVE_INFINITY)).isInfinite());
        Assert.assertTrue(new FieldVector2D<>(field, new Vector2D(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY)).isInfinite());
        Assert.assertTrue(new FieldVector2D<>(field, new Vector2D(Double.NEGATIVE_INFINITY, 0.0)).isInfinite());
        Assert.assertTrue(new FieldVector2D<>(field, new Vector2D(0.0, Double.NEGATIVE_INFINITY)).isInfinite());
        Assert.assertTrue(new FieldVector2D<>(field, new Vector2D(Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY)).isInfinite());
        Assert.assertFalse(FieldVector2D.getNaN(field).isInfinite());
        Assert.assertFalse(FieldVector2D.getMinusI(field).isInfinite());
    }

    private <T extends CalculusFieldElement<T>> void doTestEquals(final Field<T> field) {
        final FieldVector2D<T> u1 = new FieldVector2D<>(field, Vector2D.PLUS_I);
        final FieldVector2D<T> u2 = new FieldVector2D<>(field, Vector2D.MINUS_I.negate());
        final FieldVector2D<T> v1 = new FieldVector2D<>(field, new Vector2D(1.0, 0.001));
        final FieldVector2D<T> v2 = new FieldVector2D<>(field, new Vector2D(0.001, 1.0));
        Assert.assertEquals(u1, u1);
        Assert.assertEquals(u1, u2);
        Assert.assertNotEquals(u1, v1);
        Assert.assertNotEquals(u1, v2);
        Assert.assertNotEquals(u1, Vector2D.PLUS_I);
        Assert.assertEquals(u1.toVector2D(), Vector2D.PLUS_I);
        Assert.assertEquals(new FieldVector2D<>(Double.NaN, u1), FieldVector2D.getNaN(field));
        Assert.assertNotEquals(u1, FieldVector2D.getNaN(field));
        Assert.assertNotEquals(FieldVector2D.getNaN(field), v2);
    }

    private <T extends CalculusFieldElement<T>> void doTestHashCode(final Field<T> field) {
        Assert.assertEquals(542, FieldVector2D.getNaN(field).hashCode());
        Assert.assertEquals(1325400064, new FieldVector2D<>(field, new Vector2D(1.5, -0.5)).hashCode());
    }

    private <T extends CalculusFieldElement<T>> void doTestCrossProduct(final Field<T> field) {
        final double epsilon = 1e-10;

        FieldVector2D<T> p1 = new FieldVector2D<>(field, new Vector2D(1, 1));
        FieldVector2D<T> p2 = new FieldVector2D<>(field, new Vector2D(2, 2));

        FieldVector2D<T> p3 = new FieldVector2D<>(field, new Vector2D(3, 3));
        Assert.assertEquals(0.0, p3.crossProduct(p1, p2).getReal(), epsilon);

        FieldVector2D<T> p4 = new FieldVector2D<>(field, new Vector2D(1, 2));
        Assert.assertEquals(1.0, p4.crossProduct(p1, p2).getReal(), epsilon);

        FieldVector2D<T> p5 = new FieldVector2D<>(field, new Vector2D(2, 1));
        Assert.assertEquals(-1.0, p5.crossProduct(p1, p2).getReal(), epsilon);
        Assert.assertEquals(-1.0, p5.crossProduct(p1.toVector2D(), p2.toVector2D()).getReal(), epsilon);
    }

    private <T extends CalculusFieldElement<T>> void doTestOrientation(final Field<T> field) {
        Assert.assertTrue(FieldVector2D.orientation(new FieldVector2D<>(field, new Vector2D(0, 0)),
                                                    new FieldVector2D<>(field, new Vector2D(1, 0)),
                                                    new FieldVector2D<>(field, new Vector2D(1, 1))).getReal() > 0);
        Assert.assertTrue(FieldVector2D.orientation(new FieldVector2D<>(field, new Vector2D(1, 0)),
                                                    new FieldVector2D<>(field, new Vector2D(0, 0)),
                                                    new FieldVector2D<>(field, new Vector2D(1, 1))).getReal() < 0);
        Assert.assertEquals(0.0,
                            FieldVector2D.orientation(new FieldVector2D<>(field, new Vector2D(0, 0)),
                                                      new FieldVector2D<>(field, new Vector2D(1, 0)),
                                                      new FieldVector2D<>(field, new Vector2D(1, 0))).getReal(),
                            1.0e-15);
        Assert.assertEquals(0.0,
                            FieldVector2D.orientation(new FieldVector2D<>(field, new Vector2D(0, 0)),
                                                      new FieldVector2D<>(field, new Vector2D(1, 0)),
                                                      new FieldVector2D<>(field, new Vector2D(2, 0))).getReal(),
                            1.0e-15);
    }

    private <T extends CalculusFieldElement<T>> void check(final FieldVector2D<T> v,
                                                       final double x, final double y, final double tol) {
        Assert.assertEquals(x, v.getX().getReal(), tol);
        Assert.assertEquals(y, v.getY().getReal(), tol);
    }

}
