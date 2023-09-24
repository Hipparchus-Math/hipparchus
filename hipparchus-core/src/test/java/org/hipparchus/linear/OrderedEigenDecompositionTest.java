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

import static org.junit.Assert.assertNotNull;

import org.hipparchus.UnitTestUtils;
import org.hipparchus.complex.Complex;
import org.hipparchus.complex.ComplexField;
import org.hipparchus.util.FastMath;
import org.junit.Assert;
import org.junit.Test;

public class OrderedEigenDecompositionTest {

    /**
     *
     */
    @Test
    public void testReal() {
        // AA = [1 2;1 -3];

        RealMatrix A = MatrixUtils.createRealMatrix(new double[][] {
            { 1, 2 }, { 1, -3 }
        });

        OrderedEigenDecomposition ordEig = new OrderedEigenDecomposition(A);

        assertNotNull(ordEig.getD());
        assertNotNull(ordEig.getV());

        final double s2 = FastMath.sqrt(2);
        final double s3 = FastMath.sqrt(3);
        RealMatrix D_expected = MatrixUtils.createRealMatrix(new double[][] {
            { -1 - s2 * s3, 0 }, { 0, s2 * s3 - 1.0 }
        });

        RealMatrix V_expected = MatrixUtils.createRealMatrix(new double[][] {
            { (s3 - 3 * s2)     /  6, (s2 + 2 * s3) / 5 },
            { (4 * s3 + 3 * s2) / 12, (2 * s2 - s3) / 5 }
        });

        UnitTestUtils.assertEquals("D", D_expected, ordEig.getD(), 1.0e-15);
        UnitTestUtils.assertEquals("V", V_expected, ordEig.getV(), 1.0e-15);

        // checking definition of the decomposition A = V*D*inv(V)
        UnitTestUtils.assertEquals("A", A,
                                   ordEig.getV().multiply(ordEig.getD()).multiply(MatrixUtils.inverse(ordEig.getV())),
                                   2.0e-15);

    }

    /**
     *
     */
    @Test
    public void testImaginary() {
        // AA = [3 -2;4 -1];

        RealMatrix A = MatrixUtils.createRealMatrix(new double[][] {
            { 3, -2 }, { 4, -1 }
        });

        OrderedEigenDecomposition ordEig = new OrderedEigenDecomposition(A);

        assertNotNull(ordEig.getD());
        assertNotNull(ordEig.getV());

        RealMatrix D_expected = MatrixUtils.createRealMatrix(new double[][] {
            { 1, 2 }, { -2, 1 }
        });

        RealMatrix V_expected = MatrixUtils.createRealMatrix(new double[][] {
            { -0.5, 0.5 }, { 0, 1 }
        });

        UnitTestUtils.assertEquals("D", D_expected, ordEig.getD(), 1.0e-15);
        UnitTestUtils.assertEquals("V", V_expected, ordEig.getV(), 1.0e-15);

        // checking definition of the decomposition A = V*D*inv(V)
        UnitTestUtils.assertEquals("A", A, ordEig.getV().multiply(ordEig.getD()).multiply(MatrixUtils.inverse(ordEig.getV())), 1.0e-15);

        // checking A*V = D*V
        FieldMatrix<Complex> ac = new Array2DRowFieldMatrix<>(ComplexField.getInstance(), 2, 2);
        for (int i = 0; i < ac.getRowDimension(); ++i) {
            for (int j = 0; j < ac.getColumnDimension(); ++j) {
                ac.setEntry(i, j, new Complex(A.getEntry(i, j)));
            }
        }
        for (int i = 0; i < ordEig.getEigenvalues().length; ++i) {
        	final Complex              li  = ordEig.getEigenvalue(i);
        	final FieldVector<Complex> ei  = ordEig.getEigenvector(i);
        	final FieldVector<Complex> aei = ac.operate(ei);
        	for (int j = 0; j < ei.getDimension(); ++j) {
        		Assert.assertEquals(aei.getEntry(j).getRealPart(),      li.multiply(ei.getEntry(j)).getRealPart(),      1.0e-10);
        		Assert.assertEquals(aei.getEntry(j).getImaginaryPart(), li.multiply(ei.getEntry(j)).getImaginaryPart(), 1.0e-10);
        	}
        }

    }

    /**
    *
    */
   @Test
   public void testImaginary33() {
       // AA = [3 -2 0;4 -1 0;1 1 1];

       RealMatrix A = MatrixUtils.createRealMatrix(new double[][] {
           { 3, -2, 0 }, { 4, -1, 0 }, { 1, 1, 1 }
       });

       OrderedEigenDecomposition ordEig = new OrderedEigenDecomposition(A);

       assertNotNull(ordEig.getD());
       assertNotNull(ordEig.getV());

       RealMatrix D_expected = MatrixUtils.createRealMatrix(new double[][] {
           { 1, 2, 0 }, { -2, 1, 0}, {0, 0, 1 }
       });

       final double a = FastMath.sqrt(8.0 / 17.0);
       final double b = FastMath.sqrt(17.0) / 2.0;
       RealMatrix V_expected = MatrixUtils.createRealMatrix(new double[][] {
           { 0,        a, 0 },
           { a,        a, 0 },
           { a, -0.5 * a, b }
       });
       UnitTestUtils.assertEquals("D", D_expected, ordEig.getD(), 1.0e-15);
       UnitTestUtils.assertEquals("V", V_expected, ordEig.getV(), 1.0e-15);

       // checking definition of the decomposition A = V*D*inv(V)
       UnitTestUtils.assertEquals("A", A,
                                  ordEig.getV().multiply(ordEig.getD()).multiply(MatrixUtils.inverse(ordEig.getV())),
                                  8.0e-15);

       // checking A*V = D*V
       FieldMatrix<Complex> ac = new Array2DRowFieldMatrix<>(ComplexField.getInstance(), 3, 3);
       for (int i = 0; i < ac.getRowDimension(); ++i) {
           for (int j = 0; j < ac.getColumnDimension(); ++j) {
               ac.setEntry(i, j, new Complex(A.getEntry(i, j)));
           }
       }
       for (int i = 0; i < ordEig.getEigenvalues().length; ++i) {
       	final Complex              li  = ordEig.getEigenvalue(i);
       	final FieldVector<Complex> ei  = ordEig.getEigenvector(i);
       	final FieldVector<Complex> aei = ac.operate(ei);
       	for (int j = 0; j < ei.getDimension(); ++j) {
       		Assert.assertEquals(aei.getEntry(j).getRealPart(),      li.multiply(ei.getEntry(j)).getRealPart(),      1.0e-10);
       		Assert.assertEquals(aei.getEntry(j).getImaginaryPart(), li.multiply(ei.getEntry(j)).getImaginaryPart(), 1.0e-10);
       	}
       }

   }

   /**
   *
   */
  @Test
  public void testImaginaryNullEigenvalue() {
      // AA = [3 -2 0;4 -1 0;3 -2 0];

      RealMatrix A = MatrixUtils.createRealMatrix(new double[][] {
          { 3, -2, 0 }, { 4, -1, 0 }, { 3, -2, 0 }
      });

      OrderedEigenDecomposition ordEig = new OrderedEigenDecomposition(A);

      assertNotNull(ordEig.getD());
      assertNotNull(ordEig.getV());

      RealMatrix D_expected = MatrixUtils.createRealMatrix(new double[][] {
          { 0, 0, 0 }, { 0, 1, 2 }, { 0, -2, 1 }
      });

      final double a  = FastMath.sqrt(11.0 / 50.0);
      final double s2 = FastMath.sqrt(2.0);
      RealMatrix V_expected = MatrixUtils.createRealMatrix(new double[][] {
          { 0.0,     -a, -2 * a },
          { 0.0, -3 * a,     -a },
          { s2,      -a, -2 * a }
      });

      UnitTestUtils.assertEquals("D", D_expected, ordEig.getD(), 1.0e-15);
      UnitTestUtils.assertEquals("V", V_expected, ordEig.getV(), 2.0e-15);

      // checking definition of the decomposition A = V*D*inv(V)
      UnitTestUtils.assertEquals("A", A,
                                 ordEig.getV().multiply(ordEig.getD()).multiply(MatrixUtils.inverse(ordEig.getV())),
                                 7.0e-15);

      // checking A*V = D*V
      FieldMatrix<Complex> ac = new Array2DRowFieldMatrix<>(ComplexField.getInstance(), 3, 3);
      for (int i = 0; i < ac.getRowDimension(); ++i) {
          for (int j = 0; j < ac.getColumnDimension(); ++j) {
              ac.setEntry(i, j, new Complex(A.getEntry(i, j)));
          }
      }

      for (int i = 0; i < ordEig.getEigenvalues().length; ++i) {
      	final Complex              li  = ordEig.getEigenvalue(i);
      	final FieldVector<Complex> ei  = ordEig.getEigenvector(i);
      	final FieldVector<Complex> aei = ac.operate(ei);
      	for (int j = 0; j < ei.getDimension(); ++j) {
      		Assert.assertEquals(aei.getEntry(j).getRealPart(),      li.multiply(ei.getEntry(j)).getRealPart(),      1.0e-10);
      		Assert.assertEquals(aei.getEntry(j).getImaginaryPart(), li.multiply(ei.getEntry(j)).getImaginaryPart(), 1.0e-10);
      	}
      }

  }

}
