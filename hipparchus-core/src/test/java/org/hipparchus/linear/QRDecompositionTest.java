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
import org.junit.jupiter.api.Test;

import java.util.Random;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;


class QRDecompositionTest {
    private double[][] testData3x3NonSingular = {
            { 12, -51, 4 },
            { 6, 167, -68 },
            { -4, 24, -41 }, };

    private double[][] testData3x3Singular = {
            { 1, 4, 7, },
            { 2, 5, 8, },
            { 3, 6, 9, }, };

    private double[][] testData3x4 = {
            { 12, -51, 4, 1 },
            { 6, 167, -68, 2 },
            { -4, 24, -41, 3 }, };

    private double[][] testData4x3 = {
            { 12, -51, 4, },
            { 6, 167, -68, },
            { -4, 24, -41, },
            { -5, 34, 7, }, };

    private static final double entryTolerance = 10e-16;

    private static final double normTolerance = 10e-14;

    /** test dimensions */
    @Test
    void testDimensions() {
        checkDimension(MatrixUtils.createRealMatrix(testData3x3NonSingular));

        checkDimension(MatrixUtils.createRealMatrix(testData4x3));

        checkDimension(MatrixUtils.createRealMatrix(testData3x4));

        Random r = new Random(643895747384642l);
        int    p = (5 * BlockRealMatrix.BLOCK_SIZE) / 4;
        int    q = (7 * BlockRealMatrix.BLOCK_SIZE) / 4;
        checkDimension(createTestMatrix(r, p, q));
        checkDimension(createTestMatrix(r, q, p));

    }

    private void checkDimension(RealMatrix m) {
        int rows = m.getRowDimension();
        int columns = m.getColumnDimension();
        QRDecomposition qr = new QRDecomposition(m);
        assertEquals(rows,    qr.getQ().getRowDimension());
        assertEquals(rows,    qr.getQ().getColumnDimension());
        assertEquals(rows,    qr.getR().getRowDimension());
        assertEquals(columns, qr.getR().getColumnDimension());
    }

    /** test A = QR */
    @Test
    void testAEqualQR() {
        checkAEqualQR(MatrixUtils.createRealMatrix(testData3x3NonSingular));

        checkAEqualQR(MatrixUtils.createRealMatrix(testData3x3Singular));

        checkAEqualQR(MatrixUtils.createRealMatrix(testData3x4));

        checkAEqualQR(MatrixUtils.createRealMatrix(testData4x3));

        Random r = new Random(643895747384642l);
        int    p = (5 * BlockRealMatrix.BLOCK_SIZE) / 4;
        int    q = (7 * BlockRealMatrix.BLOCK_SIZE) / 4;
        checkAEqualQR(createTestMatrix(r, p, q));

        checkAEqualQR(createTestMatrix(r, q, p));

    }

    private void checkAEqualQR(RealMatrix m) {
        QRDecomposition qr = new QRDecomposition(m);
        double norm = qr.getQ().multiply(qr.getR()).subtract(m).getNorm1();
        assertEquals(0, norm, normTolerance);
    }

    /** test the orthogonality of Q */
    @Test
    void testQOrthogonal() {
        checkQOrthogonal(MatrixUtils.createRealMatrix(testData3x3NonSingular));

        checkQOrthogonal(MatrixUtils.createRealMatrix(testData3x3Singular));

        checkQOrthogonal(MatrixUtils.createRealMatrix(testData3x4));

        checkQOrthogonal(MatrixUtils.createRealMatrix(testData4x3));

        Random r = new Random(643895747384642l);
        int    p = (5 * BlockRealMatrix.BLOCK_SIZE) / 4;
        int    q = (7 * BlockRealMatrix.BLOCK_SIZE) / 4;
        checkQOrthogonal(createTestMatrix(r, p, q));

        checkQOrthogonal(createTestMatrix(r, q, p));

    }

    private void checkQOrthogonal(RealMatrix m) {
        QRDecomposition qr = new QRDecomposition(m);
        RealMatrix eye = MatrixUtils.createRealIdentityMatrix(m.getRowDimension());
        double norm = qr.getQT().multiply(qr.getQ()).subtract(eye).getNorm1();
        assertEquals(0, norm, normTolerance);
    }

    /** test that R is upper triangular */
    @Test
    void testRUpperTriangular() {
        RealMatrix matrix = MatrixUtils.createRealMatrix(testData3x3NonSingular);
        checkUpperTriangular(new QRDecomposition(matrix).getR());

        matrix = MatrixUtils.createRealMatrix(testData3x3Singular);
        checkUpperTriangular(new QRDecomposition(matrix).getR());

        matrix = MatrixUtils.createRealMatrix(testData3x4);
        checkUpperTriangular(new QRDecomposition(matrix).getR());

        matrix = MatrixUtils.createRealMatrix(testData4x3);
        checkUpperTriangular(new QRDecomposition(matrix).getR());

        Random r = new Random(643895747384642l);
        int    p = (5 * BlockRealMatrix.BLOCK_SIZE) / 4;
        int    q = (7 * BlockRealMatrix.BLOCK_SIZE) / 4;
        matrix = createTestMatrix(r, p, q);
        checkUpperTriangular(new QRDecomposition(matrix).getR());

        matrix = createTestMatrix(r, p, q);
        checkUpperTriangular(new QRDecomposition(matrix).getR());

    }

    private void checkUpperTriangular(RealMatrix m) {
        m.walkInOptimizedOrder(new DefaultRealMatrixPreservingVisitor() {
            @Override
            public void visit(int row, int column, double value) {
                if (column < row) {
                    assertEquals(0.0, value, entryTolerance);
                }
            }
        });
    }

    /** test that H is trapezoidal */
    @Test
    void testHTrapezoidal() {
        RealMatrix matrix = MatrixUtils.createRealMatrix(testData3x3NonSingular);
        checkTrapezoidal(new QRDecomposition(matrix).getH());

        matrix = MatrixUtils.createRealMatrix(testData3x3Singular);
        checkTrapezoidal(new QRDecomposition(matrix).getH());

        matrix = MatrixUtils.createRealMatrix(testData3x4);
        checkTrapezoidal(new QRDecomposition(matrix).getH());

        matrix = MatrixUtils.createRealMatrix(testData4x3);
        checkTrapezoidal(new QRDecomposition(matrix).getH());

        Random r = new Random(643895747384642l);
        int    p = (5 * BlockRealMatrix.BLOCK_SIZE) / 4;
        int    q = (7 * BlockRealMatrix.BLOCK_SIZE) / 4;
        matrix = createTestMatrix(r, p, q);
        checkTrapezoidal(new QRDecomposition(matrix).getH());

        matrix = createTestMatrix(r, p, q);
        checkTrapezoidal(new QRDecomposition(matrix).getH());

    }

    private void checkTrapezoidal(RealMatrix m) {
        m.walkInOptimizedOrder(new DefaultRealMatrixPreservingVisitor() {
            @Override
            public void visit(int row, int column, double value) {
                if (column > row) {
                    assertEquals(0.0, value, entryTolerance);
                }
            }
        });
    }

    /** test matrices values */
    @Test
    void testMatricesValues() {
        QRDecomposition qr =
            new QRDecomposition(MatrixUtils.createRealMatrix(testData3x3NonSingular));
        RealMatrix qRef = MatrixUtils.createRealMatrix(new double[][] {
                { -12.0 / 14.0,   69.0 / 175.0,  -58.0 / 175.0 },
                {  -6.0 / 14.0, -158.0 / 175.0,    6.0 / 175.0 },
                {   4.0 / 14.0,  -30.0 / 175.0, -165.0 / 175.0 }
        });
        RealMatrix rRef = MatrixUtils.createRealMatrix(new double[][] {
                { -14.0,  -21.0, 14.0 },
                {   0.0, -175.0, 70.0 },
                {   0.0,    0.0, 35.0 }
        });
        RealMatrix hRef = MatrixUtils.createRealMatrix(new double[][] {
                { 26.0 / 14.0, 0.0, 0.0 },
                {  6.0 / 14.0, 648.0 / 325.0, 0.0 },
                { -4.0 / 14.0,  36.0 / 325.0, 2.0 }
        });

        // check values against known references
        RealMatrix q = qr.getQ();
        assertEquals(0, q.subtract(qRef).getNorm1(), 1.0e-13);
        RealMatrix qT = qr.getQT();
        assertEquals(0, qT.subtract(qRef.transpose()).getNorm1(), 1.0e-13);
        RealMatrix r = qr.getR();
        assertEquals(0, r.subtract(rRef).getNorm1(), 1.0e-13);
        RealMatrix h = qr.getH();
        assertEquals(0, h.subtract(hRef).getNorm1(), 1.0e-13);

        // check the same cached instance is returned the second time
        assertTrue(q == qr.getQ());
        assertTrue(r == qr.getR());
        assertTrue(h == qr.getH());

    }

    @Test
    void testNonInvertible() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            QRDecomposition qr =
                new QRDecomposition(MatrixUtils.createRealMatrix(testData3x3Singular));
            qr.getSolver().getInverse();
        });
    }

    @Test
    void testInvertTallSkinny() {
        RealMatrix a     = MatrixUtils.createRealMatrix(testData4x3);
        DecompositionSolver solver = new QRDecomposition(a).getSolver();
        RealMatrix pinv  = solver.getInverse();
        assertEquals(0, pinv.multiply(a).subtract(MatrixUtils.createRealIdentityMatrix(3)).getNorm1(), 1.0e-6);
        assertEquals(testData4x3.length,    solver.getRowDimension());
        assertEquals(testData4x3[0].length, solver.getColumnDimension());
    }

    @Test
    void testInvertShortWide() {
        RealMatrix a = MatrixUtils.createRealMatrix(testData3x4);
        DecompositionSolver solver = new QRDecomposition(a).getSolver();
        RealMatrix pinv  = solver.getInverse();
        assertEquals(0, a.multiply(pinv).subtract(MatrixUtils.createRealIdentityMatrix(3)).getNorm1(), 1.0e-6);
        assertEquals(0, pinv.multiply(a).getSubMatrix(0, 2, 0, 2).subtract(MatrixUtils.createRealIdentityMatrix(3)).getNorm1(), 1.0e-6);
        assertEquals(testData3x4.length,    solver.getRowDimension());
        assertEquals(testData3x4[0].length, solver.getColumnDimension());
    }

    private RealMatrix createTestMatrix(final Random r, final int rows, final int columns) {
        RealMatrix m = MatrixUtils.createRealMatrix(rows, columns);
        m.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor(){
            @Override
            public double visit(int row, int column, double value) {
                return 2.0 * r.nextDouble() - 1.0;
            }
        });
        return m;
    }

    @Test
    void testQRSingular() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            final RealMatrix a = MatrixUtils.createRealMatrix(new double[][]{
                {1, 6, 4}, {2, 4, -1}, {-1, 2, 5}
            });
            final RealVector b = new ArrayRealVector(new double[]{5, 6, 1});
            new QRDecomposer(1.0e-15).decompose(a).solve(b);
        });
    }

}
