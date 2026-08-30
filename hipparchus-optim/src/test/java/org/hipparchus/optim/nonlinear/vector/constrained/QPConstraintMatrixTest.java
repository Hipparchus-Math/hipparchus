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
package org.hipparchus.optim.nonlinear.vector.constrained;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

public class QPConstraintMatrixTest {

    @Test
    public void testVariables() {
        try {
            new QPConstraintMatrix(-6,
                                   MatrixUtils.createRealMatrix(3, 3),
                                   MatrixUtils.createRealMatrix(3, 3),
                                   MatrixUtils.createRealMatrix(3, 3));
            Assertions.fail("an exception should have been thrown");
        }  catch (final MathIllegalArgumentException e) {
            Assertions.assertEquals(LocalizedCoreFormats.NEGATIVE_VALUE, e.getSpecifier());
            Assertions.assertEquals(-6, (Integer) e.getParts()[0]);
        }
    }

    @Test
    public void testDimensions() {
        try {
            new QPConstraintMatrix(4,
                                   MatrixUtils.createRealMatrix(3, 3),
                                   MatrixUtils.createRealMatrix(3, 3),
                                   MatrixUtils.createRealMatrix(3, 3));
            Assertions.fail("an exception should have been thrown");
        }  catch (final MathIllegalArgumentException e) {
            Assertions.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, e.getSpecifier());
            Assertions.assertEquals(3, (Integer) e.getParts()[0]);
            Assertions.assertEquals(4, (Integer) e.getParts()[1]);
        }
    }

    @Test
    public void testGetters() {
        final QPConstraintMatrix cm =
                new QPConstraintMatrix(2,
                                       MatrixUtils.createRealMatrix(new double[][] {
                                           {  1.0,  2.0 }
                                       }),
                                       MatrixUtils.createRealMatrix(new double[][] {
                                           {  3.0,  4.0 },
                                           {  5.0,  6.0 }
                                       }),
                                       MatrixUtils.createRealMatrix(new double[][] {
                                           {  7.0,  8.0 },
                                           {  9.0, 10.0 },
                                           { 11.0, 12.0 }
                                       }));
        Assertions.assertEquals(   2, cm.getColumnDimension());
        Assertions.assertEquals(   9, cm.getRowDimension());

        // equality
        Assertions.assertEquals( 1.0, cm.getEntry( 0,  0), 1.0e-12);
        Assertions.assertEquals( 2.0, cm.getEntry( 0,  1), 1.0e-12);
        Assertions.assertArrayEquals(new double[] { 1.0, 2.0 }, cm.getRow( 0), 1.0e-12);

        // inequality
        Assertions.assertEquals( 3.0, cm.getEntry( 1,  0), 1.0e-12);
        Assertions.assertEquals( 4.0, cm.getEntry( 1,  1), 1.0e-12);
        Assertions.assertEquals( 5.0, cm.getEntry( 2,  0), 1.0e-12);
        Assertions.assertEquals( 6.0, cm.getEntry( 2,  1), 1.0e-12);
        Assertions.assertArrayEquals(new double[] { 3.0, 4.0 }, cm.getRow( 1), 1.0e-12);
        Assertions.assertArrayEquals(new double[] { 5.0, 6.0 }, cm.getRow( 2), 1.0e-12);

        // lower bound
        Assertions.assertEquals( 7.0, cm.getEntry( 3,  0), 1.0e-12);
        Assertions.assertEquals( 8.0, cm.getEntry( 3,  1), 1.0e-12);
        Assertions.assertEquals( 9.0, cm.getEntry( 4,  0), 1.0e-12);
        Assertions.assertEquals(10.0, cm.getEntry( 4,  1), 1.0e-12);
        Assertions.assertEquals(11.0, cm.getEntry( 5,  0), 1.0e-12);
        Assertions.assertEquals(12.0, cm.getEntry( 5,  1), 1.0e-12);
        Assertions.assertArrayEquals(new double[] {  7.0,  8.0 }, cm.getRow( 3), 1.0e-12);
        Assertions.assertArrayEquals(new double[] {  9.0, 10.0 }, cm.getRow( 4), 1.0e-12);
        Assertions.assertArrayEquals(new double[] { 11.0, 12.0 }, cm.getRow( 5), 1.0e-12);

        // upper bound
        Assertions.assertEquals( -7.0, cm.getEntry( 6,  0), 1.0e-12);
        Assertions.assertEquals( -8.0, cm.getEntry( 6,  1), 1.0e-12);
        Assertions.assertEquals( -9.0, cm.getEntry( 7,  0), 1.0e-12);
        Assertions.assertEquals(-10.0, cm.getEntry( 7,  1), 1.0e-12);
        Assertions.assertEquals(-11.0, cm.getEntry( 8,  0), 1.0e-12);
        Assertions.assertEquals(-12.0, cm.getEntry( 8,  1), 1.0e-12);
        Assertions.assertArrayEquals(new double[] {  -7.0,  -8.0 }, cm.getRow( 6), 1.0e-12);
        Assertions.assertArrayEquals(new double[] {  -9.0, -10.0 }, cm.getRow( 7), 1.0e-12);
        Assertions.assertArrayEquals(new double[] { -11.0, -12.0 }, cm.getRow( 8), 1.0e-12);

    }

    @Test
    public void testSetEntry() {
        try {
            final QPConstraintMatrix cm =
                new QPConstraintMatrix(2,
                                       MatrixUtils.createRealMatrix(new double[][] {
                                           {  1.0,  2.0 }
                                       }),
                                       MatrixUtils.createRealMatrix(new double[][] {
                                           {  3.0,  4.0 },
                                           {  5.0,  6.0 }
                                       }),
                                       MatrixUtils.createRealMatrix(new double[][] {
                                           {  7.0,  8.0 },
                                           {  9.0, 10.0 },
                                           { 11.0, 12.0 }
                                       }));
            cm.setEntry(1, 2, 3.5);
            Assertions.fail("an exception should have been thrown");
        } catch (UnsupportedOperationException upe) {
            // expected
        }
    }

    @Test
    public void testIndependentMatrix() {
        final QPConstraintMatrix cm =
            new QPConstraintMatrix(2,
                                   MatrixUtils.createRealMatrix(new double[][] {
                                       {  1.0,  2.0 }
                                   }),
                                   MatrixUtils.createRealMatrix(new double[][] {
                                       {  3.0,  4.0 },
                                       {  5.0,  6.0 }
                                   }),
                                   MatrixUtils.createRealMatrix(new double[][] {
                                       {  7.0,  8.0 },
                                       {  9.0, 10.0 },
                                       { 11.0, 12.0 }
                                   }));
        final RealMatrix independent = cm.createMatrix(9, 2);
        Assertions.assertEquals(10.0, cm.getEntry(4, 1), 1.0e-12);
        Assertions.assertEquals( 0.0, independent.getEntry(4, 1), 1.0e-12);
    }

    @Test
    public void testCopy() {
        final QPConstraintMatrix cm =
            new QPConstraintMatrix(2,
                                   MatrixUtils.createRealMatrix(new double[][] {
                                       {  1.0,  2.0 }
                                   }),
                                   MatrixUtils.createRealMatrix(new double[][] {
                                       {  3.0,  4.0 },
                                       {  5.0,  6.0 }
                                   }),
                                   MatrixUtils.createRealMatrix(new double[][] {
                                       {  7.0,  8.0 },
                                       {  9.0, 10.0 },
                                       { 11.0, 12.0 }
                                   }));
        final RealMatrix copy = cm.copy();
        Assertions.assertEquals(10.0, cm.getEntry(4, 1), 1.0e-12);
        Assertions.assertEquals(10.0, copy.getEntry(4, 1), 1.0e-12);
        copy.setEntry(4, 1, 3.5);
        Assertions.assertEquals(3.5, copy.getEntry(4, 1), 1.0e-12);
    }

}
