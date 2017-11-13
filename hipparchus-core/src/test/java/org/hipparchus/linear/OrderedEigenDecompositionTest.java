/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.text.DecimalFormat;

import org.junit.Test;

public class OrderedEigenDecompositionTest {

    static private final RealMatrixFormat f =
                    new RealMatrixFormat("", "", "\n", "", "", "\t\t",
                                         new DecimalFormat(" ##############0.0000;-##############0.0000"));

    /**
     *
     */
    @Test
    public void testOrdEigenDecomposition_Real() {
        // AA = [1 2;1 -3];

        RealMatrix A = MatrixUtils.createRealMatrix(new double[][] {
            { 1, 2 }, { 1, -3 }
        });

        OrderedEigenDecomposition ordEig = new OrderedEigenDecomposition(A);

        assertNotNull(ordEig.getD());
        assertNotNull(ordEig.getV());

        RealMatrix D_expected = MatrixUtils.createRealMatrix(new double[][] {
            { -3.4495, 0 }, { 0, 1.4495 }
        });

        RealMatrix V_expected = MatrixUtils.createRealMatrix(new double[][] {
            { -0.4184, 0.9757 }, { 0.9309, 0.2193 }
        });

        assertEquals(f.format(D_expected), f.format(ordEig.getD()));
        assertEquals(f.format(V_expected), f.format(ordEig.getV()));

        // checking definition of the decomposition A = V*D*inv(V)
        assertEquals(f.format(A),
                     f.format(ordEig.getV().multiply(ordEig.getD()).multiply(MatrixUtils.inverse(ordEig.getV())).scalarAdd(0)));
    }

    /**
     *
     */
    @Test
    public void testOrdEigenDecomposition_Imaginary() {
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

        assertEquals(f.format(D_expected), f.format(ordEig.getD()));
        assertEquals(f.format(V_expected), f.format(ordEig.getV()));

        // checking definition of the decomposition A = V*D*inv(V)
        assertEquals(f.format(A),
                     f.format(ordEig.getV().multiply(ordEig.getD()).multiply(MatrixUtils.inverse(ordEig.getV())).scalarAdd(0)));
    }

    /**
     *
     */
    @Test
    public void testOrdEigenDecomposition_Imaginary_33() {
        // AA = [3 -2 0;4 -1 0;1 1 1];

        RealMatrix A = MatrixUtils.createRealMatrix(new double[][] {
            { 3, -2, 0 }, { 4, -1, 0 }, { 3, -2, 0 }
        });

        OrderedEigenDecomposition ordEig = new OrderedEigenDecomposition(A);

        assertNotNull(ordEig.getD());
        assertNotNull(ordEig.getV());

        RealMatrix D_expected = MatrixUtils.createRealMatrix(new double[][] {
            { 0, 0, 0 }, { 0, 1, 2 }, { 0, -2, 1 }
        });

        RealMatrix V_expected = MatrixUtils.createRealMatrix(new double[][] {
            { -0.00001, -0.4690, -0.9381 }, { -0.0, -1.4071, -0.4690 },
            { 1.4142, -0.4690, -0.9381 }
        });

        assertEquals(f.format(D_expected), f.format(ordEig.getD()));
        assertEquals(f.format(V_expected), f.format(ordEig.getV()));

        // checking definition of the decomposition A = V*D*inv(V)
        assertEquals(A.subtract(ordEig.getV().multiply(ordEig.getD()).multiply(MatrixUtils.inverse(ordEig.getV()))).getFrobeniusNorm(),
                     0,
                     1E-10);
    }

}
