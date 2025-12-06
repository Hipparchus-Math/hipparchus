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
package org.hipparchus.complex;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;
import org.junit.jupiter.api.Test;

import java.util.Random;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

class QuaternionTest {
    /** Epsilon for double comparison. */
    private static final double EPS = Math.ulp(1d);
    /** Epsilon for double comparison. */
    private static final double COMPARISON_EPS = 1e-14;

    @Test
    final void testAccessors1() {
        final double q0 = 2;
        final double q1 = 5.4;
        final double q2 = 17;
        final double q3 = 0.0005;
        final Quaternion q = new Quaternion(q0, q1, q2, q3);

        assertEquals(q0, q.getQ0(), 0);
        assertEquals(q1, q.getQ1(), 0);
        assertEquals(q2, q.getQ2(), 0);
        assertEquals(q3, q.getQ3(), 0);
    }

    @Test
    final void testAccessors2() {
        final double q0 = 2;
        final double q1 = 5.4;
        final double q2 = 17;
        final double q3 = 0.0005;
        final Quaternion q = new Quaternion(q0, q1, q2, q3);

        final double sP = q.getScalarPart();
        final double[] vP = q.getVectorPart();

        assertEquals(q0, sP, 0);
        assertEquals(q1, vP[0], 0);
        assertEquals(q2, vP[1], 0);
        assertEquals(q3, vP[2], 0);
    }

    @Test
    final void testAccessors3() {
        final double q0 = 2;
        final double q1 = 5.4;
        final double q2 = 17;
        final double q3 = 0.0005;
        final Quaternion q = new Quaternion(q0, new double[] { q1, q2, q3 });

        final double sP = q.getScalarPart();
        final double[] vP = q.getVectorPart();

        assertEquals(q0, sP, 0);
        assertEquals(q1, vP[0], 0);
        assertEquals(q2, vP[1], 0);
        assertEquals(q3, vP[2], 0);
    }

    @Test
    void testWrongDimension() {
        assertThrows(MathIllegalArgumentException.class, () -> new Quaternion(new double[]{1, 2}));
    }

    @Test
    final void testConjugate() {
        final double q0 = 2;
        final double q1 = 5.4;
        final double q2 = 17;
        final double q3 = 0.0005;
        final Quaternion q = new Quaternion(q0, q1, q2, q3);

        final Quaternion qConjugate = q.getConjugate();

        assertEquals(q0, qConjugate.getQ0(), 0);
        assertEquals(-q1, qConjugate.getQ1(), 0);
        assertEquals(-q2, qConjugate.getQ2(), 0);
        assertEquals(-q3, qConjugate.getQ3(), 0);
    }

    @Test
    final void testProductQuaternionQuaternion() {

        // Case : analytic test case

        final Quaternion qA = new Quaternion(1, 0.5, -3, 4);
        final Quaternion qB = new Quaternion(6, 2, 1, -9);
        final Quaternion qResult = Quaternion.multiply(qA, qB);

        assertEquals(44, qResult.getQ0(), EPS);
        assertEquals(28, qResult.getQ1(), EPS);
        assertEquals(-4.5, qResult.getQ2(), EPS);
        assertEquals(21.5, qResult.getQ3(), EPS);

        // Conjugate of the product of two quaternions and product of their conjugates :
        // Conj(qA * qB) = Conj(qB) * Conj(qA)

        final Quaternion conjugateOfProduct = qB.getConjugate().multiply(qA.getConjugate());
        final Quaternion productOfConjugate = (qA.multiply(qB)).getConjugate();

        assertEquals(conjugateOfProduct.getQ0(), productOfConjugate.getQ0(), EPS);
        assertEquals(conjugateOfProduct.getQ1(), productOfConjugate.getQ1(), EPS);
        assertEquals(conjugateOfProduct.getQ2(), productOfConjugate.getQ2(), EPS);
        assertEquals(conjugateOfProduct.getQ3(), productOfConjugate.getQ3(), EPS);
    }

    @Test
    final void testProductQuaternionVector() {

        // Case : Product between a vector and a quaternion : QxV

        final Quaternion quaternion = new Quaternion(4, 7, -1, 2);
        final double[] vector = {2.0, 1.0, 3.0};
        final Quaternion qResultQxV = Quaternion.multiply(quaternion, new Quaternion(vector));

        assertEquals(-19, qResultQxV.getQ0(), EPS);
        assertEquals(3, qResultQxV.getQ1(), EPS);
        assertEquals(-13, qResultQxV.getQ2(), EPS);
        assertEquals(21, qResultQxV.getQ3(), EPS);

        // comparison with the result given by the formula :
        // qResult = (- vectorQ . vector) + (scalarQ * vector + vectorQ ^ vector)

        final double[] vectorQ = quaternion.getVectorPart();

        final double scalarPartRefQxV = -MathArrays.linearCombination(vectorQ, vector);
        assertEquals(scalarPartRefQxV, qResultQxV.getScalarPart(), EPS);

        // Case : Product between a vector and a quaternion : VxQ

        final Quaternion qResultVxQ = Quaternion.multiply(new Quaternion(vector), quaternion);

        assertEquals(-19, qResultVxQ.getQ0(), EPS);
        assertEquals(13, qResultVxQ.getQ1(), EPS);
        assertEquals(21, qResultVxQ.getQ2(), EPS);
        assertEquals(3, qResultVxQ.getQ3(), EPS);

        // comparison with the result given by the formula :
        // qResult = (- vector . vectorQ) + (scalarQ * vector + vector ^ vectorQ)

        final double scalarPartRefVxQ = -MathArrays.linearCombination(vectorQ, vector);
        assertEquals(scalarPartRefVxQ, qResultVxQ.getScalarPart(), EPS);

    }

    @Test
    final void testDotProductQuaternionQuaternion() {
        // expected output
        final double expected = -6.;
        // inputs
        final Quaternion q1 = new Quaternion(1, 2, 2, 1);
        final Quaternion q2 = new Quaternion(3, -2, -1, -3);

        final double actual1 = Quaternion.dotProduct(q1, q2);
        final double actual2 = q1.dotProduct(q2);

        assertEquals(expected, actual1, EPS);
        assertEquals(expected, actual2, EPS);
    }

    @Test
    final void testScalarMultiplyDouble() {
        // expected outputs
        final double w = 1.6;
        final double x = -4.8;
        final double y = 11.20;
        final double z = 2.56;
        // inputs
        final Quaternion q1 = new Quaternion(0.5, -1.5, 3.5, 0.8);
        final double a = 3.2;

        final Quaternion q = q1.multiply(a);

        assertEquals(w, q.getQ0(), COMPARISON_EPS);
        assertEquals(x, q.getQ1(), COMPARISON_EPS);
        assertEquals(y, q.getQ2(), COMPARISON_EPS);
        assertEquals(z, q.getQ3(), COMPARISON_EPS);
    }

    @Test
    final void testAddQuaternionQuaternion() {
        // expected outputs
        final double w = 4;
        final double x = -1;
        final double y = 2;
        final double z = -4;
        // inputs
        final Quaternion q1 = new Quaternion(1., 2., -2., -1.);
        final Quaternion q2 = new Quaternion(3., -3., 4., -3.);

        final Quaternion qa = Quaternion.add(q1, q2);
        final Quaternion qb = q1.add(q2);

        assertEquals(w, qa.getQ0(), EPS);
        assertEquals(x, qa.getQ1(), EPS);
        assertEquals(y, qa.getQ2(), EPS);
        assertEquals(z, qa.getQ3(), EPS);

        assertEquals(w, qb.getQ0(), EPS);
        assertEquals(x, qb.getQ1(), EPS);
        assertEquals(y, qb.getQ2(), EPS);
        assertEquals(z, qb.getQ3(), EPS);
    }

    @Test
    final void testSubtractQuaternionQuaternion() {
        // expected outputs
        final double w = -2.;
        final double x = 5.;
        final double y = -6.;
        final double z = 2.;
        // inputs
        final Quaternion q1 = new Quaternion(1., 2., -2., -1.);
        final Quaternion q2 = new Quaternion(3., -3., 4., -3.);

        final Quaternion qa = Quaternion.subtract(q1, q2);
        final Quaternion qb = q1.subtract(q2);

        assertEquals(w, qa.getQ0(), EPS);
        assertEquals(x, qa.getQ1(), EPS);
        assertEquals(y, qa.getQ2(), EPS);
        assertEquals(z, qa.getQ3(), EPS);

        assertEquals(w, qb.getQ0(), EPS);
        assertEquals(x, qb.getQ1(), EPS);
        assertEquals(y, qb.getQ2(), EPS);
        assertEquals(z, qb.getQ3(), EPS);
}

    @Test
    final void testNorm() {

        final double q0 = 2;
        final double q1 = 1;
        final double q2 = -4;
        final double q3 = 3;
        final Quaternion q = new Quaternion(q0, q1, q2, q3);

        final double norm = q.getNorm();

        assertEquals(FastMath.sqrt(30), norm, 0);

        final double normSquareRef = Quaternion.multiply(q, q.getConjugate()).getScalarPart();
        assertEquals(FastMath.sqrt(normSquareRef), norm, 0);
    }

    @Test
    final void testNormalize() {

        final Quaternion q = new Quaternion(2, 1, -4, -2);

        final Quaternion versor = q.normalize();

        assertEquals(2.0 / 5.0, versor.getQ0(), 0);
        assertEquals(1.0 / 5.0, versor.getQ1(), 0);
        assertEquals(-4.0 / 5.0, versor.getQ2(), 0);
        assertEquals(-2.0 / 5.0, versor.getQ3(), 0);

        assertEquals(1, versor.getNorm(), 0);
    }

    @Test
    final void testNormalizeFail() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final Quaternion zeroQ = new Quaternion(0, 0, 0, 0);
            zeroQ.normalize();
        });
    }

    @Test
    final void testObjectEquals() {
        final double one = 1;
        final Quaternion q1 = new Quaternion(one, one, one, one);
        assertEquals(q1, q1);

        final Quaternion q2 = new Quaternion(one, one, one, one);
        assertEquals(q2, q1);

        final Quaternion q3 = new Quaternion(one, FastMath.nextUp(one), one, one);
        assertNotEquals(q3, q1);
    }

    @Test
    final void testQuaternionEquals() {
        final double inc = 1e-5;
        final Quaternion q1 = new Quaternion(2, 1, -4, -2);
        final Quaternion q2 = new Quaternion(q1.getQ0() + inc, q1.getQ1(), q1.getQ2(), q1.getQ3());
        final Quaternion q3 = new Quaternion(q1.getQ0(), q1.getQ1() + inc, q1.getQ2(), q1.getQ3());
        final Quaternion q4 = new Quaternion(q1.getQ0(), q1.getQ1(), q1.getQ2() + inc, q1.getQ3());
        final Quaternion q5 = new Quaternion(q1.getQ0(), q1.getQ1(), q1.getQ2(), q1.getQ3() + inc);

        assertFalse(q1.equals(q2, 0.9 * inc));
        assertFalse(q1.equals(q3, 0.9 * inc));
        assertFalse(q1.equals(q4, 0.9 * inc));
        assertFalse(q1.equals(q5, 0.9 * inc));

        assertTrue(q1.equals(q2, 1.1 * inc));
        assertTrue(q1.equals(q3, 1.1 * inc));
        assertTrue(q1.equals(q4, 1.1 * inc));
        assertTrue(q1.equals(q5, 1.1 * inc));
    }

    @Test
    final void testQuaternionEquals2() {
        final Quaternion q1 = new Quaternion(1, 4, 2, 3);
        final double gap = 1e-5;
        final Quaternion q2 = new Quaternion(1 + gap, 4 + gap, 2 + gap, 3 + gap);

        assertTrue(q1.equals(q2, 10 * gap));
        assertFalse(q1.equals(q2, gap));
        assertFalse(q1.equals(q2, gap / 10));
    }

    @Test
    final void testIsUnitQuaternion() {
        final Random r = new Random(48);
        final int numberOfTrials = 1000;
        for (int i = 0; i < numberOfTrials; i++) {
            final Quaternion q1 = new Quaternion(r.nextDouble(), r.nextDouble(), r.nextDouble(), r.nextDouble());
            final Quaternion q2 = q1.normalize();
            assertTrue(q2.isUnitQuaternion(COMPARISON_EPS));
        }

        final Quaternion q = new Quaternion(1, 1, 1, 1);
        assertFalse(q.isUnitQuaternion(COMPARISON_EPS));
    }

    @Test
    final void testIsPureQuaternion() {
        final Quaternion q1 = new Quaternion(0, 5, 4, 8);
        assertTrue(q1.isPureQuaternion(EPS));

        final Quaternion q2 = new Quaternion(0 - EPS, 5, 4, 8);
        assertTrue(q2.isPureQuaternion(EPS));

        final Quaternion q3 = new Quaternion(0 - 1.1 * EPS, 5, 4, 8);
        assertFalse(q3.isPureQuaternion(EPS));

        final Random r = new Random(48);
        final double[] v = {r.nextDouble(), r.nextDouble(), r.nextDouble()};
        final Quaternion q4 = new Quaternion(v);
        assertTrue(q4.isPureQuaternion(0));

        final Quaternion q5 = new Quaternion(0, v);
        assertTrue(q5.isPureQuaternion(0));
    }

    @Test
    final void testGetInverse() {
        final Quaternion q = new Quaternion(1.5, 4, 2, -2.5);

        final Quaternion inverseQ = q.getInverse();
        assertEquals(1.5 / 28.5, inverseQ.getQ0(), 0);
        assertEquals(-4.0 / 28.5, inverseQ.getQ1(), 0);
        assertEquals(-2.0 / 28.5, inverseQ.getQ2(), 0);
        assertEquals(2.5 / 28.5, inverseQ.getQ3(), 0);

        final Quaternion product = Quaternion.multiply(inverseQ, q);
        assertEquals(1, product.getQ0(), EPS);
        assertEquals(0, product.getQ1(), EPS);
        assertEquals(0, product.getQ2(), EPS);
        assertEquals(0, product.getQ3(), EPS);

        final Quaternion qNul = new Quaternion(0, 0, 0, 0);
        try {
            final Quaternion inverseQNul = qNul.getInverse();
            fail("expecting MathIllegalArgumentException but got : " + inverseQNul);
        } catch (MathIllegalArgumentException ex) {
            // expected
        }
    }

    @Test
    final void testToString() {
        final Quaternion q = new Quaternion(1, 2, 3, 4);
        assertEquals("[1.0 2.0 3.0 4.0]", q.toString());
    }
}
