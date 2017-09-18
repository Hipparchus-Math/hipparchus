/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
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

import java.util.Random;

import org.hipparchus.Field;
import org.hipparchus.RealFieldElement;
import org.hipparchus.analysis.differentiation.DSFactory;
import org.hipparchus.analysis.differentiation.DerivativeStructure;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.Decimal64Field;
import org.hipparchus.util.MathArrays;
import org.junit.Assert;
import org.junit.Test;


public class FieldQRDecompositionTest {
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

    DerivativeStructure zero = new DSFactory(1, 1).variable(0, 0);
    Field<DerivativeStructure> DSField = zero.getField();
    private static final double entryTolerance = 10e-16;

    private static final double normTolerance = 10e-14;

    /**Testing if the dimensions are correct.*/
    @Test
    public void testDimensions() {
        doTestDimensions(DSField);   
        doTestDimensions(Decimal64Field.getInstance());   
    }

    /**Testing if is impossible to solve QR.*/
    @Test(expected=MathIllegalArgumentException.class)
    public void testQRSingular(){
        QRSingular(DSField);
        QRSingular(Decimal64Field.getInstance());
    }

    /**Testing if Q is orthogonal*/
    @Test
    public void testQOrthogonal(){
        QOrthogonal(DSField);
        QOrthogonal(Decimal64Field.getInstance());
    }

    /**Testing if A = Q * R*/
    @Test
    public void testAEqualQR(){
        AEqualQR(DSField);
        AEqualQR(Decimal64Field.getInstance());
    }

    /**Testing if R is upper triangular.*/
    @Test
    public void testRUpperTriangular(){
        RUpperTriangular(DSField);
        RUpperTriangular(Decimal64Field.getInstance());
    }

    /**Testing if H is trapezoidal.*/
    @Test
    public void testHTrapezoidal(){
        HTrapezoidal(DSField);
        HTrapezoidal(Decimal64Field.getInstance());
    }
    /**Testing the values of the matrices.*/
    @Test
    public void testMatricesValues(){
        MatricesValues(DSField);
        MatricesValues(Decimal64Field.getInstance());
    }

    /**Testing if there is an error inverting a non invertible matrix.*/
    @Test(expected=MathIllegalArgumentException.class)
    public void testNonInvertible(){
        NonInvertible(DSField);
        NonInvertible(Decimal64Field.getInstance());
    }
    /**Testing to invert a tall and skinny matrix.*/
    @Test
    public void testInvertTallSkinny(){
        InvertTallSkinny(DSField);
        InvertTallSkinny(Decimal64Field.getInstance());
    }
    /**Testing to invert a short and wide matrix.*/
    @Test
    public void testInvertShortWide(){
        InvertShortWide(DSField);
        InvertShortWide(Decimal64Field.getInstance());
    }

    private <T extends RealFieldElement<T>> void doTestDimensions(Field<T> field){
        T[][] data3x3NS= convert(field, testData3x3NonSingular ); 
        T[][] data3x4= convert(field, testData3x4            );
        T[][] data4x3= convert(field, testData4x3            );
        checkDimension(MatrixUtils.createFieldMatrix( data3x3NS));

        checkDimension(MatrixUtils.createFieldMatrix( data4x3));

        checkDimension(MatrixUtils.createFieldMatrix( data3x4));

        Random r = new Random(643895747384642l);
        int    p = (5 * BlockFieldMatrix.BLOCK_SIZE) / 4;
        int    q = (7 * BlockFieldMatrix.BLOCK_SIZE) / 4;
        checkDimension(createTestMatrix(field,r, p, q));
        checkDimension(createTestMatrix(field, r, q, p));

    }

    private <T extends RealFieldElement<T>> void checkDimension(FieldMatrix<T> m) {
        int rows = m.getRowDimension();
        int columns = m.getColumnDimension();
        FieldQRDecomposition<T> qr = new FieldQRDecomposition<T>(m);
        Assert.assertEquals(rows,    qr.getQ().getRowDimension());
        Assert.assertEquals(rows,    qr.getQ().getColumnDimension());
        Assert.assertEquals(rows,    qr.getR().getRowDimension());
        Assert.assertEquals(columns, qr.getR().getColumnDimension());
    }

    private  <T extends RealFieldElement<T>> void AEqualQR(Field<T> field) {
        T[][] data3x3NS= convert(field, testData3x3NonSingular ); 
        T[][] data3x3S= convert(field, testData3x3Singular    );
        T[][] data3x4= convert(field, testData3x4            );
        T[][] data4x3= convert(field, testData4x3            );
        checkAEqualQR(MatrixUtils.createFieldMatrix( data3x3NS));

        checkAEqualQR(MatrixUtils.createFieldMatrix( data3x3S));

        checkAEqualQR(MatrixUtils.createFieldMatrix( data3x4));

        checkAEqualQR(MatrixUtils.createFieldMatrix( data4x3));

        Random r = new Random(643895747384642l);
        int    p = (5 * BlockFieldMatrix.BLOCK_SIZE) / 4;
        int    q = (7 * BlockFieldMatrix.BLOCK_SIZE) / 4;
        checkAEqualQR(createTestMatrix(field, r, p, q));

        checkAEqualQR(createTestMatrix(field, r, q, p));

    }

    private  <T extends RealFieldElement<T>> void checkAEqualQR(FieldMatrix<T> m) {
        FieldQRDecomposition<T> qr = new FieldQRDecomposition<>(m);
        T norm = norm(qr.getQ().multiply(qr.getR()).subtract(m));
        Assert.assertEquals(0, norm.getReal(), normTolerance);
    }

    private  <T extends RealFieldElement<T>> void QOrthogonal(Field<T> field) {
        T[][] data3x3NS= convert(field, testData3x3NonSingular ); 
        T[][] data3x3S= convert(field, testData3x3Singular    );
        T[][] data3x4= convert(field, testData3x4            );
        T[][] data4x3= convert(field, testData4x3            );
        checkQOrthogonal(MatrixUtils.createFieldMatrix( data3x3NS));

        checkQOrthogonal(MatrixUtils.createFieldMatrix( data3x3S));

        checkQOrthogonal(MatrixUtils.createFieldMatrix( data3x4));

        checkQOrthogonal(MatrixUtils.createFieldMatrix( data4x3));

        Random r = new Random(643895747384642l);
        int    p = (5 * BlockFieldMatrix.BLOCK_SIZE) / 4;
        int    q = (7 * BlockFieldMatrix.BLOCK_SIZE) / 4;
        checkQOrthogonal(createTestMatrix(field, r, p, q));

        checkQOrthogonal(createTestMatrix(field, r, q, p));

    }

    private  <T extends RealFieldElement<T>> void checkQOrthogonal(FieldMatrix<T> m) {
        FieldQRDecomposition<T> qr = new FieldQRDecomposition<T>(m);
        FieldMatrix<T> eye = MatrixUtils.createFieldIdentityMatrix(m.getField(),m.getRowDimension());
        T norm = norm(qr.getQT().multiply(qr.getQ()).subtract(eye));
        Assert.assertEquals(0, norm.getReal(), normTolerance);
    }

    private  <T extends RealFieldElement<T>> void RUpperTriangular(Field<T> field) {
        T[][] data3x3NS= convert(field, testData3x3NonSingular ); 
        T[][] data3x3S= convert(field, testData3x3Singular    );
        T[][] data3x4= convert(field, testData3x4            );
        T[][] data4x3= convert(field, testData4x3            );
        FieldMatrix<T> matrix = MatrixUtils.createFieldMatrix( data3x3NS);
        checkUpperTriangular(new FieldQRDecomposition<T>(matrix).getR());

        matrix = MatrixUtils.createFieldMatrix( data3x3S);
        checkUpperTriangular(new FieldQRDecomposition<T>(matrix).getR());

        matrix = MatrixUtils.createFieldMatrix( data3x4);
        checkUpperTriangular(new FieldQRDecomposition<T>(matrix).getR());

        matrix = MatrixUtils.createFieldMatrix( data4x3);
        checkUpperTriangular(new FieldQRDecomposition<T>(matrix).getR());

        Random r = new Random(643895747384642l);
        int    p = (5 * BlockFieldMatrix.BLOCK_SIZE) / 4;
        int    q = (7 * BlockFieldMatrix.BLOCK_SIZE) / 4;
        matrix = createTestMatrix(field, r, p, q);
        checkUpperTriangular(new FieldQRDecomposition<T>(matrix).getR());

        matrix = createTestMatrix(field, r, p, q);
        checkUpperTriangular(new FieldQRDecomposition<T>(matrix).getR());

    }

    private  <T extends RealFieldElement<T>> void checkUpperTriangular(FieldMatrix<T> m) {
        m.walkInOptimizedOrder(new DefaultFieldMatrixPreservingVisitor<T>(m.getField().getZero()) {
            @Override
            public void visit(int row, int column, T value) {
                if (column < row) {
                    Assert.assertEquals(0.0, value.getReal(), entryTolerance);
                }
            }
        });
    }

    private <T extends RealFieldElement<T>> void HTrapezoidal(Field<T> field) {
        T[][] data3x3NS= convert(field, testData3x3NonSingular ); 
        T[][] data3x3S= convert(field, testData3x3Singular    );
        T[][] data3x4= convert(field, testData3x4            );
        T[][] data4x3= convert(field, testData4x3            );
        FieldMatrix<T> matrix = MatrixUtils.createFieldMatrix( data3x3NS);
        checkTrapezoidal(new FieldQRDecomposition<T>(matrix).getH());

        matrix = MatrixUtils.createFieldMatrix( data3x3S);
        checkTrapezoidal(new FieldQRDecomposition<T>(matrix).getH());

        matrix = MatrixUtils.createFieldMatrix( data3x4);
        checkTrapezoidal(new FieldQRDecomposition<T>(matrix).getH());

        matrix = MatrixUtils.createFieldMatrix( data4x3);
        checkTrapezoidal(new FieldQRDecomposition<T>(matrix).getH());

        Random r = new Random(643895747384642l);
        int    p = (5 * BlockFieldMatrix.BLOCK_SIZE) / 4;
        int    q = (7 * BlockFieldMatrix.BLOCK_SIZE) / 4;
        matrix = createTestMatrix(field, r, p, q);
        checkTrapezoidal(new FieldQRDecomposition<T>(matrix).getH());

        matrix = createTestMatrix(field, r, p, q);
        checkTrapezoidal(new FieldQRDecomposition<T>(matrix).getH());

    }

    private  <T extends RealFieldElement<T>> void checkTrapezoidal(FieldMatrix<T> m) {
        m.walkInOptimizedOrder(new DefaultFieldMatrixPreservingVisitor<T>(m.getField().getZero()) {
            @Override
            public void visit(int row, int column, T value) {
                if (column > row) {
                    Assert.assertEquals(0.0, value.getReal(), entryTolerance);
                }
            }
        });
    }
    
    private  <T extends RealFieldElement<T>> void MatricesValues(Field<T> field) {
        T[][] data3x3NS= convert(field, testData3x3NonSingular ); 
        FieldQRDecomposition<T> qr =
            new FieldQRDecomposition<T>(MatrixUtils.createFieldMatrix( data3x3NS));
        FieldMatrix<T> qRef = MatrixUtils.createFieldMatrix( convert(field,new double[][] {
                { -12.0 / 14.0,   69.0 / 175.0,  -58.0 / 175.0 },
                {  -6.0 / 14.0, -158.0 / 175.0,    6.0 / 175.0 },
                {   4.0 / 14.0,  -30.0 / 175.0, -165.0 / 175.0 }
        }));
        FieldMatrix<T> rRef = MatrixUtils.createFieldMatrix( convert(field, new double[][] {
                { -14.0,  -21.0, 14.0 },
                {   0.0, -175.0, 70.0 },
                {   0.0,    0.0, 35.0 }
        }));
        FieldMatrix<T> hRef = MatrixUtils.createFieldMatrix( convert(field,new double[][] {
                { 26.0 / 14.0, 0.0, 0.0 },
                {  6.0 / 14.0, 648.0 / 325.0, 0.0 },
                { -4.0 / 14.0,  36.0 / 325.0, 2.0 }
        }));

        // check values against known references
        FieldMatrix<T> q = qr.getQ();
        Assert.assertEquals(0, norm(q.subtract(qRef)).getReal(), 1.0e-13);
        FieldMatrix<T> qT = qr.getQT();
        Assert.assertEquals(0, norm(qT.subtract(qRef.transpose())).getReal(), 1.0e-13);
        FieldMatrix<T> r = qr.getR();
        Assert.assertEquals(0, norm(r.subtract(rRef)).getReal(), 1.0e-13);
        FieldMatrix<T> h = qr.getH();
        Assert.assertEquals(0, norm(h.subtract(hRef)).getReal(), 1.0e-13);

        // check the same cached instance is returned the second time
        Assert.assertTrue(q == qr.getQ());
        Assert.assertTrue(r == qr.getR());
        Assert.assertTrue(h == qr.getH());

    }

    private  <T extends RealFieldElement<T>> void NonInvertible(Field<T> field) {
        T[][] data3x3S= convert(field, testData3x3Singular    );
        FieldQRDecomposition<T> qr =
            new FieldQRDecomposition<T>(MatrixUtils.createFieldMatrix( data3x3S));
        qr.getSolver().getInverse();
    }

    private  <T extends RealFieldElement<T>> void InvertTallSkinny(Field<T> field) {
        T[][] data4x3= convert(field, testData4x3            );
        FieldMatrix<T> a     = MatrixUtils.createFieldMatrix(data4x3);
        FieldMatrix<T> pinv  = new FieldQRDecomposition<T>(a).getSolver().getInverse();
        Assert.assertEquals(0, norm(pinv.multiply(a).subtract(MatrixUtils.createFieldIdentityMatrix(field, 3))).getReal(), 1.0e-6);
    }

    private  <T extends RealFieldElement<T>> void InvertShortWide(Field<T> field) {
        T[][] data3x4= convert(field, testData3x4            );
        FieldMatrix<T> a = MatrixUtils.createFieldMatrix( data3x4);
        FieldMatrix<T> pinv  = new FieldQRDecomposition<T>(a).getSolver().getInverse();
        Assert.assertEquals(0,norm( a.multiply(pinv).subtract(MatrixUtils.createFieldIdentityMatrix(field, 3))).getReal(), 1.0e-6);
        Assert.assertEquals(0,norm( pinv.multiply(a).getSubMatrix(0, 2, 0, 2).subtract(MatrixUtils.createFieldIdentityMatrix(field, 3))).getReal(), 1.0e-6);
    }

    private  <T extends RealFieldElement<T>> FieldMatrix<T> createTestMatrix(Field<T> field, final Random r, final int rows, final int columns) {
        FieldMatrix<T> m = MatrixUtils.createFieldMatrix(field, rows, columns);
        m.walkInOptimizedOrder(new DefaultFieldMatrixChangingVisitor<T>(field.getOne()){
            @Override
            public T visit(int row, int column,T value) {
                return field.getZero().add(2.0 * r.nextDouble() - 1.0);
            }
        });
        return m;
    }

    private  <T extends RealFieldElement<T>> void QRSingular(Field<T> field) {
        final FieldMatrix<T> a = MatrixUtils.createFieldMatrix(convert( field, new double[][] {
            { 1, 6, 4 }, { 2, 4, -1 }, { -1, 2, 5 }
        }));
        T[] vv = MathArrays.buildArray(field,3);
        vv[0] = field.getZero().add(5);
        vv[1] = field.getZero().add(6);
        vv[2] = field.getZero().add(1);
        
        final FieldVector<T> b = new ArrayFieldVector<T>(field, vv);
        new FieldQRDecomposition<T>(a, field.getZero().add(1.0e-15)).getSolver().solve(b);
    }
    
    private <T extends RealFieldElement<T>> T norm(FieldMatrix<T> FM ){
        return walkInColumnOrder(FM, new FieldMatrixPreservingVisitor<T>() {

            /** Last row index. */
            private double endRow;

            /** Sum of absolute values on one column. */
            private T columnSum;

            /** Maximal sum across all columns. */
            private T maxColSum;

            /** {@inheritDoc} */
            @Override
            public void start(final int rows, final int columns,
                              final int startRow, final int endRow,
                              final int startColumn, final int endColumn) {
                this.endRow = endRow;
                columnSum   = FM.getField().getZero();
                maxColSum   = FM.getField().getZero();
            }

            /** {@inheritDoc} */
            @Override
            public void visit(final int row, final int column, final T value) {
                columnSum = columnSum.add(value).abs();
                if (row == endRow) {
                    maxColSum = (maxColSum.getReal() > columnSum.getReal()) ? maxColSum : columnSum ;
                    columnSum = FM.getField().getZero();
                }
            }

            /** {@inheritDoc} */
            @Override
            public T end() {
                return maxColSum;
            }

        });
    }
    
    private <T extends RealFieldElement<T>> T walkInColumnOrder(FieldMatrix<T> FM, FieldMatrixPreservingVisitor<T> visitor) {
        final int rows    = FM.getRowDimension();
        final int columns = FM.getColumnDimension();
        visitor.start(rows, columns, 0, rows - 1, 0, columns - 1);
        for (int column = 0; column < columns; ++column) {
            for (int row = 0; row < rows; ++row) {
                visitor.visit(row, column, FM.getEntry(row, column));
            }
        }
        return visitor.end();
    }
    
    private <T extends RealFieldElement<T>> T[][] convert(Field<T> field, double[][] value){
        T[][] res = MathArrays.buildArray(field, value.length, value[0].length);
        for (int ii = 0; ii < (value.length); ii++){
            for (int jj = 0; jj < (value[0].length); jj++){
                res[ii][jj] = field.getZero().add(value[ii][jj]);
            }
        }
        return res;
        
    }


}
