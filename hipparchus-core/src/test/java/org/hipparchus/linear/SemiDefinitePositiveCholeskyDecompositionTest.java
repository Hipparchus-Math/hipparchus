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

package org.hipparchus.linear;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.junit.Assert;
import org.junit.Test;

public class SemiDefinitePositiveCholeskyDecompositionTest {

    private double[][] testData = new double[][] {
        {  1,  2,   4,   7,  11 },
        {  2, 13,  23,  38,  58 },
        {  4, 23,  77, 122, 182 },
        {  7, 38, 122, 294, 430 },
        { 11, 58, 182, 430, 855 }
    };

    /** test dimensions */
    @Test
    public void testDimensions() {
        SemiDefinitePositiveCholeskyDecomposition llt =
            new SemiDefinitePositiveCholeskyDecomposition(MatrixUtils.createRealMatrix(testData));
        Assert.assertEquals(testData.length, llt.getL().getRowDimension());
        Assert.assertEquals(testData.length, llt.getL().getColumnDimension());
        Assert.assertEquals(testData.length, llt.getLT().getRowDimension());
        Assert.assertEquals(testData.length, llt.getLT().getColumnDimension());
    }

    /** test non-square matrix */
    @Test(expected = MathIllegalArgumentException.class)
    public void testNonSquare() {
        new SemiDefinitePositiveCholeskyDecomposition(MatrixUtils.createRealMatrix(new double[3][2]));
    }

    /** test negative definite matrix */
    @Test(expected = MathIllegalArgumentException.class)
    public void testNotPositiveDefinite() {
        new SemiDefinitePositiveCholeskyDecomposition(MatrixUtils.createRealMatrix(new double[][] {
                { -14, 11,  13, 15,  24 },
                {  11, 34,  13,  8,  25 },
                { -13, 13,  14, 15,  21 },
                {  15,  8, -15, 18,  23 },
                {  24, 25,  21, 23, -45 }
        }));
    }

    /** test A = LLT */
    @Test
    public void testAEqualLLT() {
        RealMatrix matrix = MatrixUtils.createRealMatrix(testData);
        SemiDefinitePositiveCholeskyDecomposition llt = new SemiDefinitePositiveCholeskyDecomposition(matrix);
        RealMatrix l  = llt.getL();
        RealMatrix lt = llt.getLT();
        double norm = l.multiply(lt).subtract(matrix).getNorm1();
        Assert.assertEquals(0, norm, 1.0e-15);
    }

    /** test that L is lower triangular */
    @Test
    public void testLLowerTriangular() {
        RealMatrix matrix = MatrixUtils.createRealMatrix(testData);
        RealMatrix l = new SemiDefinitePositiveCholeskyDecomposition(matrix).getL();
        for (int i = 0; i < l.getRowDimension(); i++) {
            for (int j = i + 1; j < l.getColumnDimension(); j++) {
                Assert.assertEquals(0.0, l.getEntry(i, j), 0.0);
            }
        }
    }

    /** test that LT is transpose of L */
    @Test
    public void testLTTransposed() {
        RealMatrix matrix = MatrixUtils.createRealMatrix(testData);
        SemiDefinitePositiveCholeskyDecomposition llt = new SemiDefinitePositiveCholeskyDecomposition(matrix);
        RealMatrix l  = llt.getL();
        RealMatrix lt = llt.getLT();
        double norm = l.subtract(lt.transpose()).getNorm1();
        Assert.assertEquals(0, norm, 1.0e-15);
    }

}
