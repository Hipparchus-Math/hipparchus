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

package org.hipparchus.geometry.euclidean.threed;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.geometry.LocalizedGeometryFormats;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

class RotationTest {

    @Test
    void testIssue304Cardan() {
        for (final RotationConvention convention : RotationConvention.values()) {
            for (final RotationOrder order : Arrays.asList(RotationOrder.XYZ,
                                                           RotationOrder.XZY,
                                                           RotationOrder.YXZ,
                                                           RotationOrder.YZX,
                                                           RotationOrder.ZXY,
                                                           RotationOrder.ZYX)) {

                // first singularity
                Rotation singularPlus = new Rotation(order, convention, 0.0, MathUtils.SEMI_PI, 0.125);
                assertEquals(0.0, singularPlus.getAngles(order, convention)[0], 1.0e-16);
                assertEquals(MathUtils.SEMI_PI, singularPlus.getAngles(order, convention)[1], 1.0e-16);
                assertEquals(0.125, singularPlus.getAngles(order, convention)[2], 1.0e-16);

                // second singularity
                Rotation singularMinus = new Rotation(order, convention, 0.0, -MathUtils.SEMI_PI, 0.125);
                assertEquals(0.0, singularMinus.getAngles(order, convention)[0], 1.0e-16);
                assertEquals(-MathUtils.SEMI_PI, singularMinus.getAngles(order, convention)[1], 1.0e-16);
                assertEquals(0.125, singularMinus.getAngles(order, convention)[2], 1.0e-16);

            }
        }
    }

    @Test
    void testIssue304Euler() {
        for (final RotationConvention convention : RotationConvention.values()) {
            for (final RotationOrder order : Arrays.asList(RotationOrder.XYX,
                                                           RotationOrder.XZX,
                                                           RotationOrder.YXY,
                                                           RotationOrder.YZY,
                                                           RotationOrder.ZXZ,
                                                           RotationOrder.ZYZ)) {

                // first singularity
                Rotation singularZero = new Rotation(order, convention, 0.125, 0.0, 0.0);
                assertEquals(0.125, singularZero.getAngles(order, convention)[0], 1.0e-16);
                assertEquals(0.0, singularZero.getAngles(order, convention)[1], 1.0e-16);
                assertEquals(0.0, singularZero.getAngles(order, convention)[2], 1.0e-16);

                // second singularity
                Rotation singularPi = new Rotation(order, convention, 0.125, FastMath.PI, 0.0);
                assertEquals(0.125, singularPi.getAngles(order, convention)[0], 1.0e-16);
                assertEquals(FastMath.PI, singularPi.getAngles(order, convention)[1], 1.0e-16);
                assertEquals(0.0, singularPi.getAngles(order, convention)[2], 1.0e-16);

            }
        }
    }

    @Test
    void testIdentity() {

    Rotation r = Rotation.IDENTITY;
    checkVector(r.applyTo(Vector3D.PLUS_I), Vector3D.PLUS_I);
    checkVector(r.applyTo(Vector3D.PLUS_J), Vector3D.PLUS_J);
    checkVector(r.applyTo(Vector3D.PLUS_K), Vector3D.PLUS_K);
    checkAngle(r.getAngle(), 0);

    r = new Rotation(-1, 0, 0, 0, false);
    checkVector(r.applyTo(Vector3D.PLUS_I), Vector3D.PLUS_I);
    checkVector(r.applyTo(Vector3D.PLUS_J), Vector3D.PLUS_J);
    checkVector(r.applyTo(Vector3D.PLUS_K), Vector3D.PLUS_K);
    checkAngle(r.getAngle(), 0);

    r = new Rotation(42, 0, 0, 0, true);
    checkVector(r.applyTo(Vector3D.PLUS_I), Vector3D.PLUS_I);
    checkVector(r.applyTo(Vector3D.PLUS_J), Vector3D.PLUS_J);
    checkVector(r.applyTo(Vector3D.PLUS_K), Vector3D.PLUS_K);
    checkAngle(r.getAngle(), 0);

  }

    @Test
    void testAxisAngleVectorOperator() throws MathIllegalArgumentException {

    Rotation r = new Rotation(new Vector3D(10, 10, 10), 2 * FastMath.PI / 3, RotationConvention.VECTOR_OPERATOR);
    checkVector(r.applyTo(Vector3D.PLUS_I), Vector3D.PLUS_J);
    checkVector(r.applyTo(Vector3D.PLUS_J), Vector3D.PLUS_K);
    checkVector(r.applyTo(Vector3D.PLUS_K), Vector3D.PLUS_I);
    double s = 1 / FastMath.sqrt(3);
    checkVector(r.getAxis(RotationConvention.VECTOR_OPERATOR), new Vector3D( s,  s,  s));
    checkVector(r.getAxis(RotationConvention.FRAME_TRANSFORM), new Vector3D(-s, -s, -s));
    checkAngle(r.getAngle(), 2 * FastMath.PI / 3);

    try {
      new Rotation(new Vector3D(0, 0, 0), 2 * FastMath.PI / 3, RotationConvention.VECTOR_OPERATOR);
      fail("an exception should have been thrown");
    } catch (MathIllegalArgumentException e) {
        assertEquals(LocalizedGeometryFormats.ZERO_NORM_FOR_ROTATION_AXIS, e.getSpecifier());
    }

    r = new Rotation(Vector3D.PLUS_K, 1.5 * FastMath.PI, RotationConvention.VECTOR_OPERATOR);
    checkVector(r.getAxis(RotationConvention.VECTOR_OPERATOR), new Vector3D(0, 0, -1));
    checkVector(r.getAxis(RotationConvention.FRAME_TRANSFORM), new Vector3D(0, 0, +1));
    checkAngle(r.getAngle(), MathUtils.SEMI_PI);

    r = new Rotation(Vector3D.PLUS_J, FastMath.PI, RotationConvention.VECTOR_OPERATOR);
    checkVector(r.getAxis(RotationConvention.VECTOR_OPERATOR), Vector3D.PLUS_J);
    checkVector(r.getAxis(RotationConvention.FRAME_TRANSFORM), Vector3D.MINUS_J);
    checkAngle(r.getAngle(), FastMath.PI);

    checkVector(Rotation.IDENTITY.getAxis(RotationConvention.VECTOR_OPERATOR), Vector3D.PLUS_I);
    checkVector(Rotation.IDENTITY.getAxis(RotationConvention.FRAME_TRANSFORM), Vector3D.MINUS_I);

  }

    @Test
    void testAxisAngleFrameTransform() throws MathIllegalArgumentException {

    Rotation r = new Rotation(new Vector3D(10, 10, 10), 2 * FastMath.PI / 3, RotationConvention.FRAME_TRANSFORM);
    checkVector(r.applyTo(Vector3D.PLUS_I), Vector3D.PLUS_K);
    checkVector(r.applyTo(Vector3D.PLUS_J), Vector3D.PLUS_I);
    checkVector(r.applyTo(Vector3D.PLUS_K), Vector3D.PLUS_J);
    double s = 1 / FastMath.sqrt(3);
    checkVector(r.getAxis(RotationConvention.FRAME_TRANSFORM), new Vector3D( s,  s,  s));
    checkVector(r.getAxis(RotationConvention.VECTOR_OPERATOR), new Vector3D(-s, -s, -s));
    checkAngle(r.getAngle(), 2 * FastMath.PI / 3);

    try {
      new Rotation(new Vector3D(0, 0, 0), 2 * FastMath.PI / 3, RotationConvention.FRAME_TRANSFORM);
      fail("an exception should have been thrown");
    } catch (MathIllegalArgumentException e) {
        assertEquals(LocalizedGeometryFormats.ZERO_NORM_FOR_ROTATION_AXIS, e.getSpecifier());
    }

    r = new Rotation(Vector3D.PLUS_K, 1.5 * FastMath.PI, RotationConvention.FRAME_TRANSFORM);
    checkVector(r.getAxis(RotationConvention.FRAME_TRANSFORM), new Vector3D(0, 0, -1));
    checkVector(r.getAxis(RotationConvention.VECTOR_OPERATOR), new Vector3D(0, 0, +1));
    checkAngle(r.getAngle(), MathUtils.SEMI_PI);

    r = new Rotation(Vector3D.PLUS_J, FastMath.PI, RotationConvention.FRAME_TRANSFORM);
    checkVector(r.getAxis(RotationConvention.FRAME_TRANSFORM), Vector3D.PLUS_J);
    checkVector(r.getAxis(RotationConvention.VECTOR_OPERATOR), Vector3D.MINUS_J);
    checkAngle(r.getAngle(), FastMath.PI);

    checkVector(Rotation.IDENTITY.getAxis(RotationConvention.FRAME_TRANSFORM), Vector3D.MINUS_I);
    checkVector(Rotation.IDENTITY.getAxis(RotationConvention.VECTOR_OPERATOR), Vector3D.PLUS_I);

  }

    @Test
    void testWrongMatrix() {
      checkWrongMatrix(new double[2][2]);
      checkWrongMatrix(new double[][] { new double[2], new double[3], new double[3]});
      checkWrongMatrix(new double[][] { new double[3], new double[2], new double[3]});
      checkWrongMatrix(new double[][] { new double[3], new double[3], new double[2]});
  }

  private void checkWrongMatrix(final double[][] m) {
      try {
          new Rotation(m, 0.001);
          fail("an exception should have been thrown");
      } catch (MathIllegalArgumentException miae) {
          assertEquals(LocalizedGeometryFormats.ROTATION_MATRIX_DIMENSIONS, miae.getSpecifier());
      }
  }

    @Test
    void testRevertVectorOperator() {
    Rotation r = new Rotation(0.001, 0.36, 0.48, 0.8, true);
    Rotation reverted = r.revert();
    checkRotation(r.compose(reverted, RotationConvention.VECTOR_OPERATOR), 1, 0, 0, 0);
    checkRotation(reverted.compose(r, RotationConvention.VECTOR_OPERATOR), 1, 0, 0, 0);
    assertEquals(r.getAngle(), reverted.getAngle(), 1.0e-12);
    assertEquals(-1,
                        Vector3D.dotProduct(r.getAxis(RotationConvention.VECTOR_OPERATOR),
                                           reverted.getAxis(RotationConvention.VECTOR_OPERATOR)),
                        1.0e-12);
  }

    @Test
    void testRevertFrameTransform() {
    Rotation r = new Rotation(0.001, 0.36, 0.48, 0.8, true);
    Rotation reverted = r.revert();
    checkRotation(r.compose(reverted, RotationConvention.FRAME_TRANSFORM), 1, 0, 0, 0);
    checkRotation(reverted.compose(r, RotationConvention.FRAME_TRANSFORM), 1, 0, 0, 0);
    assertEquals(r.getAngle(), reverted.getAngle(), 1.0e-12);
    assertEquals(-1,
                        Vector3D.dotProduct(r.getAxis(RotationConvention.FRAME_TRANSFORM),
                                           reverted.getAxis(RotationConvention.FRAME_TRANSFORM)),
                        1.0e-12);
  }

    @Test
    void testVectorOnePair() throws MathRuntimeException {

    Vector3D u = new Vector3D(3, 2, 1);
    Vector3D v = new Vector3D(-4, 2, 2);
    Rotation r = new Rotation(u, v);
    checkVector(r.applyTo(u.scalarMultiply(v.getNorm())), v.scalarMultiply(u.getNorm()));

    checkAngle(new Rotation(u, u.negate()).getAngle(), FastMath.PI);

    try {
        new Rotation(u, Vector3D.ZERO);
        fail("an exception should have been thrown");
    } catch (MathRuntimeException e) {
        // expected behavior
    }

  }

    @Test
    void testVectorTwoPairs() throws MathRuntimeException {

    Vector3D u1 = new Vector3D(3, 0, 0);
    Vector3D u2 = new Vector3D(0, 5, 0);
    Vector3D v1 = new Vector3D(0, 0, 2);
    Vector3D v2 = new Vector3D(-2, 0, 2);
    Rotation r = new Rotation(u1, u2, v1, v2);
    checkVector(r.applyTo(Vector3D.PLUS_I), Vector3D.PLUS_K);
    checkVector(r.applyTo(Vector3D.PLUS_J), Vector3D.MINUS_I);

    r = new Rotation(u1, u2, u1.negate(), u2.negate());
    Vector3D axis = r.getAxis(RotationConvention.VECTOR_OPERATOR);
    if (Vector3D.dotProduct(axis, Vector3D.PLUS_K) > 0) {
      checkVector(axis, Vector3D.PLUS_K);
    } else {
      checkVector(axis, Vector3D.MINUS_K);
    }
    checkAngle(r.getAngle(), FastMath.PI);

    double sqrt = FastMath.sqrt(2) / 2;
    r = new Rotation(Vector3D.PLUS_I,  Vector3D.PLUS_J,
                     new Vector3D(0.5, 0.5,  sqrt),
                     new Vector3D(0.5, 0.5, -sqrt));
    checkRotation(r, sqrt, 0.5, 0.5, 0);

    r = new Rotation(u1, u2, u1, Vector3D.crossProduct(u1, u2));
    checkRotation(r, sqrt, -sqrt, 0, 0);

    checkRotation(new Rotation(u1, u2, u1, u2), 1, 0, 0, 0);

    try {
        new Rotation(u1, u2, Vector3D.ZERO, v2);
        fail("an exception should have been thrown");
    } catch (MathRuntimeException e) {
      // expected behavior
    }

  }

    @Test
    void testMatrix()
        throws MathIllegalArgumentException {

    try {
      new Rotation(new double[][] {
                     { 0.0, 1.0, 0.0 },
                     { 1.0, 0.0, 0.0 }
                   }, 1.0e-7);
      fail("Expecting MathIllegalArgumentException");
    } catch (MathIllegalArgumentException nrme) {
      // expected behavior
    }

    try {
      new Rotation(new double[][] {
                     {  0.445888,  0.797184, -0.407040 },
                     {  0.821760, -0.184320,  0.539200 },
                     { -0.354816,  0.574912,  0.737280 }
                   }, 1.0e-7);
      fail("Expecting MathIllegalArgumentException");
    } catch (MathIllegalArgumentException nrme) {
      // expected behavior
    }

    try {
        new Rotation(new double[][] {
                       {  0.4,  0.8, -0.4 },
                       { -0.4,  0.6,  0.7 },
                       {  0.8, -0.2,  0.5 }
                     }, 1.0e-15);
        fail("Expecting MathIllegalArgumentException");
      } catch (MathIllegalArgumentException nrme) {
        // expected behavior
      }

    checkRotation(new Rotation(new double[][] {
                                 {  0.445888,  0.797184, -0.407040 },
                                 { -0.354816,  0.574912,  0.737280 },
                                 {  0.821760, -0.184320,  0.539200 }
                               }, 1.0e-10),
                  0.8, 0.288, 0.384, 0.36);

    checkRotation(new Rotation(new double[][] {
                                 {  0.539200,  0.737280,  0.407040 },
                                 {  0.184320, -0.574912,  0.797184 },
                                 {  0.821760, -0.354816, -0.445888 }
                              }, 1.0e-10),
                  0.36, 0.8, 0.288, 0.384);

    checkRotation(new Rotation(new double[][] {
                                 { -0.445888,  0.797184, -0.407040 },
                                 {  0.354816,  0.574912,  0.737280 },
                                 {  0.821760,  0.184320, -0.539200 }
                               }, 1.0e-10),
                  0.384, 0.36, 0.8, 0.288);

    checkRotation(new Rotation(new double[][] {
                                 { -0.539200,  0.737280,  0.407040 },
                                 { -0.184320, -0.574912,  0.797184 },
                                 {  0.821760,  0.354816,  0.445888 }
                               }, 1.0e-10),
                  0.288, 0.384, 0.36, 0.8);

    double[][] m1 = { { 0.0, 1.0, 0.0 },
                      { 0.0, 0.0, 1.0 },
                      { 1.0, 0.0, 0.0 } };
    Rotation r = new Rotation(m1, 1.0e-7);
    checkVector(r.applyTo(Vector3D.PLUS_I), Vector3D.PLUS_K);
    checkVector(r.applyTo(Vector3D.PLUS_J), Vector3D.PLUS_I);
    checkVector(r.applyTo(Vector3D.PLUS_K), Vector3D.PLUS_J);

    double[][] m2 = { { 0.83203, -0.55012, -0.07139 },
                      { 0.48293,  0.78164, -0.39474 },
                      { 0.27296,  0.29396,  0.91602 } };
    r = new Rotation(m2, 1.0e-12);

    double[][] m3 = r.getMatrix();
    double d00 = m2[0][0] - m3[0][0];
    double d01 = m2[0][1] - m3[0][1];
    double d02 = m2[0][2] - m3[0][2];
    double d10 = m2[1][0] - m3[1][0];
    double d11 = m2[1][1] - m3[1][1];
    double d12 = m2[1][2] - m3[1][2];
    double d20 = m2[2][0] - m3[2][0];
    double d21 = m2[2][1] - m3[2][1];
    double d22 = m2[2][2] - m3[2][2];

    assertTrue(FastMath.abs(d00) < 6.0e-6);
    assertTrue(FastMath.abs(d01) < 6.0e-6);
    assertTrue(FastMath.abs(d02) < 6.0e-6);
    assertTrue(FastMath.abs(d10) < 6.0e-6);
    assertTrue(FastMath.abs(d11) < 6.0e-6);
    assertTrue(FastMath.abs(d12) < 6.0e-6);
    assertTrue(FastMath.abs(d20) < 6.0e-6);
    assertTrue(FastMath.abs(d21) < 6.0e-6);
    assertTrue(FastMath.abs(d22) < 6.0e-6);

    assertTrue(FastMath.abs(d00) > 4.0e-7);
    assertTrue(FastMath.abs(d01) > 4.0e-7);
    assertTrue(FastMath.abs(d02) > 4.0e-7);
    assertTrue(FastMath.abs(d10) > 4.0e-7);
    assertTrue(FastMath.abs(d11) > 4.0e-7);
    assertTrue(FastMath.abs(d12) > 4.0e-7);
    assertTrue(FastMath.abs(d20) > 4.0e-7);
    assertTrue(FastMath.abs(d21) > 4.0e-7);
    assertTrue(FastMath.abs(d22) > 4.0e-7);

    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        double m3tm3 = m3[i][0] * m3[j][0]
                     + m3[i][1] * m3[j][1]
                     + m3[i][2] * m3[j][2];
        if (i == j) {
          assertTrue(FastMath.abs(m3tm3 - 1.0) < 1.0e-10);
        } else {
          assertTrue(FastMath.abs(m3tm3) < 1.0e-10);
        }
      }
    }

    checkVector(r.applyTo(Vector3D.PLUS_I),
                new Vector3D(m3[0][0], m3[1][0], m3[2][0]));
    checkVector(r.applyTo(Vector3D.PLUS_J),
                new Vector3D(m3[0][1], m3[1][1], m3[2][1]));
    checkVector(r.applyTo(Vector3D.PLUS_K),
                new Vector3D(m3[0][2], m3[1][2], m3[2][2]));

    double[][] m4 = { { 1.0,  0.0,  0.0 },
                      { 0.0, -1.0,  0.0 },
                      { 0.0,  0.0, -1.0 } };
    r = new Rotation(m4, 1.0e-7);
    checkAngle(r.getAngle(), FastMath.PI);

    try {
      double[][] m5 = { { 0.0, 0.0, 1.0 },
                        { 0.0, 1.0, 0.0 },
                        { 1.0, 0.0, 0.0 } };
      r = new Rotation(m5, 1.0e-7);
      fail("got " + r + ", should have caught an exception");
    } catch (MathIllegalArgumentException e) {
      // expected
    }

  }

    @Test
    void testAngles()
        throws MathIllegalStateException {

      for (RotationConvention convention : RotationConvention.values()) {
          RotationOrder[] CardanOrders = {
              RotationOrder.XYZ, RotationOrder.XZY, RotationOrder.YXZ,
              RotationOrder.YZX, RotationOrder.ZXY, RotationOrder.ZYX
          };

          for (final RotationOrder cardanOrder : CardanOrders) {
              for (double alpha1 = 0.1; alpha1 < 6.2; alpha1 += 0.3) {
                  for (double alpha2 = -1.55; alpha2 < 1.55; alpha2 += 0.3) {
                      for (double alpha3 = 0.1; alpha3 < 6.2; alpha3 += 0.3) {
                          Rotation r      = new Rotation(cardanOrder, convention, alpha1, alpha2, alpha3);
                          double[] angles = r.getAngles(cardanOrder, convention);
                          checkAngle(angles[0], alpha1);
                          checkAngle(angles[1], alpha2);
                          checkAngle(angles[2], alpha3);
                      }
                  }
              }
          }

          RotationOrder[] EulerOrders = {
              RotationOrder.XYX, RotationOrder.XZX, RotationOrder.YXY,
              RotationOrder.YZY, RotationOrder.ZXZ, RotationOrder.ZYZ
          };

          for (final RotationOrder eulerOrder : EulerOrders) {
              for (double alpha1 = 0.1; alpha1 < 6.2; alpha1 += 0.3) {
                  for (double alpha2 = 0.05; alpha2 < 3.1; alpha2 += 0.3) {
                      for (double alpha3 = 0.1; alpha3 < 6.2; alpha3 += 0.3) {
                          Rotation r = new Rotation(eulerOrder, convention,
                                                    alpha1, alpha2, alpha3);
                          double[] angles = r.getAngles(eulerOrder, convention);
                          checkAngle(angles[0], alpha1);
                          checkAngle(angles[1], alpha2);
                          checkAngle(angles[2], alpha3);
                      }
                  }
              }
          }
      }

  }

    @Test
    void testSingularities() {

      for (RotationConvention convention : RotationConvention.values()) {
          RotationOrder[] CardanOrders = {
              RotationOrder.XYZ, RotationOrder.XZY, RotationOrder.YXZ,
              RotationOrder.YZX, RotationOrder.ZXY, RotationOrder.ZYX
          };

          double[] singularCardanAngle = {
              -FastMath.PI / 2, -FastMath.PI / 2 + 1.0e-12, -FastMath.PI / 2 + 1.0e-10,
               FastMath.PI / 2 - 1.0e-10, FastMath.PI / 2 - 1.0e-12, FastMath.PI / 2
          };
          for (final RotationOrder cardanOrder : CardanOrders) {
              for (final double v : singularCardanAngle) {
                  Rotation r = new Rotation(cardanOrder, convention, 0.1, v, 0.3);
                  assertEquals(v, r.getAngles(cardanOrder, convention)[1], 4.5e-16);
              }
          }

          RotationOrder[] EulerOrders = {
              RotationOrder.XYX, RotationOrder.XZX, RotationOrder.YXY,
              RotationOrder.YZY, RotationOrder.ZXZ, RotationOrder.ZYZ
          };

          double[] singularEulerAngle = { 0, 1.0e-12, 1.0e-10, FastMath.PI - 1.0e-10, FastMath.PI - 1.0e-12, FastMath.PI };
          for (final RotationOrder eulerOrder : EulerOrders) {
              for (final double v : singularEulerAngle) {
                  Rotation r = new Rotation(eulerOrder, convention, 0.1, v, 0.3);
                  assertEquals(v, r.getAngles(eulerOrder, convention)[1], 1.0e-24);
              }
          }
      }


  }

    @Test
    void testQuaternion() throws MathIllegalArgumentException {

    Rotation r1 = new Rotation(new Vector3D(2, -3, 5), 1.7, RotationConvention.VECTOR_OPERATOR);
    double n = 23.5;
    Rotation r2 = new Rotation(n * r1.getQ0(), n * r1.getQ1(),
                               n * r1.getQ2(), n * r1.getQ3(),
                               true);
    for (double x = -0.9; x < 0.9; x += 0.2) {
      for (double y = -0.9; y < 0.9; y += 0.2) {
        for (double z = -0.9; z < 0.9; z += 0.2) {
          Vector3D u = new Vector3D(x, y, z);
          checkVector(r2.applyTo(u), r1.applyTo(u));
        }
      }
    }

    r1 = new Rotation( 0.288,  0.384,  0.36,  0.8, false);
    checkRotation(r1, -r1.getQ0(), -r1.getQ1(), -r1.getQ2(), -r1.getQ3());

  }

    @Test
    void testApplyTo() throws MathIllegalArgumentException {

    Rotation r1 = new Rotation(new Vector3D(2, -3, 5), 1.7, RotationConvention.VECTOR_OPERATOR);
    Rotation r2 = new Rotation(new Vector3D(-1, 3, 2), 0.3, RotationConvention.VECTOR_OPERATOR);
    Rotation r3 = r2.applyTo(r1);

    for (double x = -0.9; x < 0.9; x += 0.2) {
      for (double y = -0.9; y < 0.9; y += 0.2) {
        for (double z = -0.9; z < 0.9; z += 0.2) {
          Vector3D u = new Vector3D(x, y, z);
          checkVector(r2.applyTo(r1.applyTo(u)), r3.applyTo(u));
        }
      }
    }

  }

    @Test
    void testComposeVectorOperator() throws MathIllegalArgumentException {

    Rotation r1 = new Rotation(new Vector3D(2, -3, 5), 1.7, RotationConvention.VECTOR_OPERATOR);
    Rotation r2 = new Rotation(new Vector3D(-1, 3, 2), 0.3, RotationConvention.VECTOR_OPERATOR);
    Rotation r3 = r2.compose(r1, RotationConvention.VECTOR_OPERATOR);

    for (double x = -0.9; x < 0.9; x += 0.2) {
      for (double y = -0.9; y < 0.9; y += 0.2) {
        for (double z = -0.9; z < 0.9; z += 0.2) {
          Vector3D u = new Vector3D(x, y, z);
          checkVector(r2.applyTo(r1.applyTo(u)), r3.applyTo(u));
        }
      }
    }

  }

    @Test
    void testComposeFrameTransform() throws MathIllegalArgumentException {

    Rotation r1 = new Rotation(new Vector3D(2, -3, 5), 1.7, RotationConvention.FRAME_TRANSFORM);
    Rotation r2 = new Rotation(new Vector3D(-1, 3, 2), 0.3, RotationConvention.FRAME_TRANSFORM);
    Rotation r3 = r2.compose(r1, RotationConvention.FRAME_TRANSFORM);
    Rotation r4 = r1.compose(r2, RotationConvention.VECTOR_OPERATOR);
    assertEquals(0.0, Rotation.distance(r3, r4), 1.0e-15);

    for (double x = -0.9; x < 0.9; x += 0.2) {
      for (double y = -0.9; y < 0.9; y += 0.2) {
        for (double z = -0.9; z < 0.9; z += 0.2) {
          Vector3D u = new Vector3D(x, y, z);
          checkVector(r1.applyTo(r2.applyTo(u)), r3.applyTo(u));
        }
      }
    }

  }

    @Test
    void testApplyInverseToRotation() throws MathIllegalArgumentException {

    Rotation r1 = new Rotation(new Vector3D(2, -3, 5), 1.7, RotationConvention.VECTOR_OPERATOR);
    Rotation r2 = new Rotation(new Vector3D(-1, 3, 2), 0.3, RotationConvention.VECTOR_OPERATOR);
    Rotation r3 = r2.applyInverseTo(r1);

    for (double x = -0.9; x < 0.9; x += 0.2) {
      for (double y = -0.9; y < 0.9; y += 0.2) {
        for (double z = -0.9; z < 0.9; z += 0.2) {
          Vector3D u = new Vector3D(x, y, z);
          checkVector(r2.applyInverseTo(r1.applyTo(u)), r3.applyTo(u));
        }
      }
    }

  }

    @Test
    void testComposeInverseVectorOperator() throws MathIllegalArgumentException {

    Rotation r1 = new Rotation(new Vector3D(2, -3, 5), 1.7, RotationConvention.VECTOR_OPERATOR);
    Rotation r2 = new Rotation(new Vector3D(-1, 3, 2), 0.3, RotationConvention.VECTOR_OPERATOR);
    Rotation r3 = r2.composeInverse(r1, RotationConvention.VECTOR_OPERATOR);

    for (double x = -0.9; x < 0.9; x += 0.2) {
      for (double y = -0.9; y < 0.9; y += 0.2) {
        for (double z = -0.9; z < 0.9; z += 0.2) {
          Vector3D u = new Vector3D(x, y, z);
          checkVector(r2.applyInverseTo(r1.applyTo(u)), r3.applyTo(u));
        }
      }
    }

  }

    @Test
    void testComposeInverseFrameTransform() throws MathIllegalArgumentException {

    Rotation r1 = new Rotation(new Vector3D(2, -3, 5), 1.7, RotationConvention.FRAME_TRANSFORM);
    Rotation r2 = new Rotation(new Vector3D(-1, 3, 2), 0.3, RotationConvention.FRAME_TRANSFORM);
    Rotation r3 = r2.composeInverse(r1, RotationConvention.FRAME_TRANSFORM);
    Rotation r4 = r1.revert().composeInverse(r2.revert(), RotationConvention.VECTOR_OPERATOR);
    assertEquals(0.0, Rotation.distance(r3, r4), 1.0e-15);

    for (double x = -0.9; x < 0.9; x += 0.2) {
      for (double y = -0.9; y < 0.9; y += 0.2) {
        for (double z = -0.9; z < 0.9; z += 0.2) {
          Vector3D u = new Vector3D(x, y, z);
          checkVector(r1.applyTo(r2.applyInverseTo(u)), r3.applyTo(u));
        }
      }
    }

  }

    @Test
    void testArray() throws MathIllegalArgumentException {

      Rotation r = new Rotation(new Vector3D(2, -3, 5), 1.7, RotationConvention.VECTOR_OPERATOR);

      for (double x = -0.9; x < 0.9; x += 0.2) {
          for (double y = -0.9; y < 0.9; y += 0.2) {
              for (double z = -0.9; z < 0.9; z += 0.2) {
                  Vector3D u = new Vector3D(x, y, z);
                  Vector3D v = r.applyTo(u);
                  double[] inOut = new double[] { x, y, z };
                  r.applyTo(inOut, inOut);
                  assertEquals(v.getX(), inOut[0], 1.0e-10);
                  assertEquals(v.getY(), inOut[1], 1.0e-10);
                  assertEquals(v.getZ(), inOut[2], 1.0e-10);
                  r.applyInverseTo(inOut, inOut);
                  assertEquals(u.getX(), inOut[0], 1.0e-10);
                  assertEquals(u.getY(), inOut[1], 1.0e-10);
                  assertEquals(u.getZ(), inOut[2], 1.0e-10);
              }
          }
      }

  }

    @Test
    void testApplyInverseTo() throws MathIllegalArgumentException {

    Rotation r = new Rotation(new Vector3D(2, -3, 5), 1.7, RotationConvention.VECTOR_OPERATOR);
    for (double lambda = 0; lambda < 6.2; lambda += 0.2) {
      for (double phi = -1.55; phi < 1.55; phi += 0.2) {
          Vector3D u = new Vector3D(FastMath.cos(lambda) * FastMath.cos(phi),
                                    FastMath.sin(lambda) * FastMath.cos(phi),
                                    FastMath.sin(phi));
          r.applyInverseTo(r.applyTo(u));
          checkVector(u, r.applyInverseTo(r.applyTo(u)));
          checkVector(u, r.applyTo(r.applyInverseTo(u)));
      }
    }

    r = Rotation.IDENTITY;
    for (double lambda = 0; lambda < 6.2; lambda += 0.2) {
      for (double phi = -1.55; phi < 1.55; phi += 0.2) {
          Vector3D u = new Vector3D(FastMath.cos(lambda) * FastMath.cos(phi),
                                    FastMath.sin(lambda) * FastMath.cos(phi),
                                    FastMath.sin(phi));
          checkVector(u, r.applyInverseTo(r.applyTo(u)));
          checkVector(u, r.applyTo(r.applyInverseTo(u)));
      }
    }

    r = new Rotation(Vector3D.PLUS_K, FastMath.PI, RotationConvention.VECTOR_OPERATOR);
    for (double lambda = 0; lambda < 6.2; lambda += 0.2) {
      for (double phi = -1.55; phi < 1.55; phi += 0.2) {
          Vector3D u = new Vector3D(FastMath.cos(lambda) * FastMath.cos(phi),
                                    FastMath.sin(lambda) * FastMath.cos(phi),
                                    FastMath.sin(phi));
          checkVector(u, r.applyInverseTo(r.applyTo(u)));
          checkVector(u, r.applyTo(r.applyInverseTo(u)));
      }
    }

  }

    @Test
    void testIssue639() throws MathRuntimeException{
      Vector3D u1 = new Vector3D(-1321008684645961.0 /  268435456.0,
                                 -5774608829631843.0 /  268435456.0,
                                 -3822921525525679.0 / 4294967296.0);
      Vector3D u2 =new Vector3D( -5712344449280879.0 /    2097152.0,
                                 -2275058564560979.0 /    1048576.0,
                                  4423475992255071.0 /      65536.0);
      Rotation rot = new Rotation(u1, u2, Vector3D.PLUS_I,Vector3D.PLUS_K);
      assertEquals( 0.6228370359608200639829222, rot.getQ0(), 1.0e-15);
      assertEquals( 0.0257707621456498790029987, rot.getQ1(), 1.0e-15);
      assertEquals(-0.0000000002503012255839931, rot.getQ2(), 1.0e-15);
      assertEquals(-0.7819270390861109450724902, rot.getQ3(), 1.0e-15);
  }

    @Test
    void testIssue801() throws MathRuntimeException {
      Vector3D u1 = new Vector3D(0.9999988431610581, -0.0015210774290851095, 0.0);
      Vector3D u2 = new Vector3D(0.0, 0.0, 1.0);

      Vector3D v1 = new Vector3D(0.9999999999999999, 0.0, 0.0);
      Vector3D v2 = new Vector3D(0.0, 0.0, -1.0);

      Rotation quat = new Rotation(u1, u2, v1, v2);
      double q2 = quat.getQ0() * quat.getQ0() +
                  quat.getQ1() * quat.getQ1() +
                  quat.getQ2() * quat.getQ2() +
                  quat.getQ3() * quat.getQ3();
      assertEquals(1.0, q2, 1.0e-14);
      assertEquals(0.0, Vector3D.angle(v1, quat.applyTo(u1)), 1.0e-14);
      assertEquals(0.0, Vector3D.angle(v2, quat.applyTo(u2)), 1.0e-14);

  }

    @Test
    void testGithubPullRequest22A() {
      final RotationOrder order = RotationOrder.ZYX;
      final double xRotation = FastMath.toDegrees(30);
      final double yRotation = FastMath.toDegrees(20);
      final double zRotation = FastMath.toDegrees(10);
      final Vector3D startingVector = Vector3D.PLUS_I;
      Vector3D appliedIndividually = startingVector;
      appliedIndividually = new Rotation(order, RotationConvention.FRAME_TRANSFORM, zRotation, 0, 0).applyTo(appliedIndividually);
      appliedIndividually = new Rotation(order, RotationConvention.FRAME_TRANSFORM, 0, yRotation, 0).applyTo(appliedIndividually);
      appliedIndividually = new Rotation(order, RotationConvention.FRAME_TRANSFORM, 0, 0, xRotation).applyTo(appliedIndividually);

      final Vector3D bad = new Rotation(order, RotationConvention.FRAME_TRANSFORM, zRotation, yRotation, xRotation).applyTo(startingVector);

      assertEquals(bad.getX(), appliedIndividually.getX(), 1e-12);
      assertEquals(bad.getY(), appliedIndividually.getY(), 1e-12);
      assertEquals(bad.getZ(), appliedIndividually.getZ(), 1e-12);
  }

    @Test
    void testGithubPullRequest22B() {
      final RotationOrder order = RotationOrder.ZYX;
      final double xRotation = FastMath.toDegrees(30);
      final double yRotation = FastMath.toDegrees(20);
      final double zRotation = FastMath.toDegrees(10);
      final Vector3D startingVector = Vector3D.PLUS_I;
      Vector3D appliedIndividually = startingVector;
      appliedIndividually = new Rotation(order, RotationConvention.FRAME_TRANSFORM, zRotation, 0, 0).applyTo(appliedIndividually);
      appliedIndividually = new Rotation(order, RotationConvention.FRAME_TRANSFORM, 0, yRotation, 0).applyTo(appliedIndividually);
      appliedIndividually = new Rotation(order, RotationConvention.FRAME_TRANSFORM, 0, 0, xRotation).applyTo(appliedIndividually);

      final Rotation r1 = new Rotation(order.getA1(), zRotation, RotationConvention.FRAME_TRANSFORM);
      final Rotation r2 = new Rotation(order.getA2(), yRotation, RotationConvention.FRAME_TRANSFORM);
      final Rotation r3 = new Rotation(order.getA3(), xRotation, RotationConvention.FRAME_TRANSFORM);
      final Rotation composite = r1.compose(r2.compose(r3,
                                                       RotationConvention.FRAME_TRANSFORM),
                                            RotationConvention.FRAME_TRANSFORM);
      final Vector3D good = composite.applyTo(startingVector);

      assertEquals(good.getX(), appliedIndividually.getX(), 1e-12);
      assertEquals(good.getY(), appliedIndividually.getY(), 1e-12);
      assertEquals(good.getZ(), appliedIndividually.getZ(), 1e-12);
  }

  private void checkVector(Vector3D v1, Vector3D v2) {
    assertTrue(v1.subtract(v2).getNorm() < 1.0e-10);
  }

  private void checkAngle(double a1, double a2) {
    assertEquals(a1, MathUtils.normalizeAngle(a2, a1), 1.0e-10);
  }

  private void checkRotation(Rotation r, double q0, double q1, double q2, double q3) {
    assertEquals(0, Rotation.distance(r, new Rotation(q0, q1, q2, q3, false)), 1.0e-12);
  }

}
