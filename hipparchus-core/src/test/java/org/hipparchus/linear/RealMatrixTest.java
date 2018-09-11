/*
 * Licensed to the Hipparchus project under one or more
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
package org.hipparchus.linear;

import org.junit.Assert;
import org.junit.Test;

public class RealMatrixTest {

    @Test
    public void testDefaultMultiplyTransposed() {
        RealMatrix a = MatrixUtils.createRealMatrix(new double[][] { {1d,2d,3d}, {2d,5d,3d}, {1d,0d,8d} });
        RealMatrix b = MatrixUtils.createRealMatrix(new double[][] { {4d, -5d, 6d} });
        RealMatrix abTRef = a.multiplyTransposed(b);
        DefaultMatrix dma = new DefaultMatrix(a);
        DefaultMatrix dmb = new DefaultMatrix(b);
        RealMatrix abT = dma.multiplyTransposed(dmb);
        Assert.assertEquals(0.0, abT.subtract(abTRef).getNorm(), 1.0e-10);
    }

    @Test
    public void testDefaultTransposeMultiply() {
        RealMatrix a = MatrixUtils.createRealMatrix(new double[][] { {1d,2d,3d}, {2d,5d,3d}, {1d,0d,8d} });
        RealMatrix b = MatrixUtils.createRealMatrix(new double[][] { {4d}, {-5d}, {6d} });
        RealMatrix aTbRef = a.transposeMultiply(b);
        DefaultMatrix dma = new DefaultMatrix(a);
        DefaultMatrix dmb = new DefaultMatrix(b);
        RealMatrix aTb = dma.transposeMultiply(dmb);
        Assert.assertEquals(0.0, aTb.subtract(aTbRef).getNorm(), 1.0e-10);
    }

    // local class that does NOT override multiplyTransposed nor transposeMultiply
    // so the default methods are called
    private class DefaultMatrix extends AbstractRealMatrix {

        RealMatrix m;
        public DefaultMatrix(RealMatrix m) {
            this.m = m;
        }

        @Override
        public int getRowDimension() {
            return m.getRowDimension();
        }

        @Override
        public int getColumnDimension() {
            return m.getColumnDimension();
        }

        @Override
        public RealMatrix createMatrix(int rowDimension, int columnDimension) {
            return m.createMatrix(rowDimension, columnDimension);
        }

        @Override
        public RealMatrix copy() {
            return m.copy();
        }

        @Override
        public double getEntry(int row, int column) {
            return m.getEntry(row, column);
        }

        @Override
        public void setEntry(int row, int column, double value) {
            m.setEntry(row, column, value);
        }

    }

}
