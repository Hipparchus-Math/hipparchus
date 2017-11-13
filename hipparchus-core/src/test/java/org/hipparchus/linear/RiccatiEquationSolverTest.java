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

import java.text.DecimalFormat;

import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealMatrixFormat;
import org.hipparchus.linear.RiccatiEquationSolverImpl;
import org.junit.Test;

public class RiccatiEquationSolverTest {

    static private final RealMatrixFormat f =
                    new RealMatrixFormat("", "", "\n", "", "", "\t\t",
                                         new DecimalFormat(" ##############0.0000;-##############0.0000"));

    @Test
    public void test_real_2_2() {
        // AA = [-3 2;1 1];
        // BB = [0;1];

        RealMatrix A = MatrixUtils.createRealMatrix(new double[][] {
            { -3, 2 }, { 1, 1 }
        });

        RealMatrix B = MatrixUtils.createRealMatrix(new double[][] {
            { 0 }, { 1 }
        });

        RealMatrix R = MatrixUtils.createRealIdentityMatrix(1);
        RealMatrix Q = MatrixUtils.createRealIdentityMatrix(2);

        RiccatiEquationSolver a = new RiccatiEquationSolverImpl(A, B, Q, R);

        RealMatrix P_expected = MatrixUtils.createRealMatrix(new double[][] {
            { 0.3221, 0.7407 }, { 0.7407, 3.2277 }
        });

        assertEquals(f.format(P_expected), f.format(a.getP()));

        // P
        // 0.3221 0.7407
        // 0.7407 3.2277
        // K
        // 0.7407 3.2277
    }

    @Test
    public void test_imaginary_2_2() {
        // AA = [3 -2;4 -1];
        // BB = [0;1];

        RealMatrix A = MatrixUtils.createRealMatrix(new double[][] {
            { 3, -2 }, { 4, -1 }
        });

        RealMatrix B = MatrixUtils.createRealMatrix(new double[][] {
            { 0 }, { 1 }
        });

        RealMatrix R = MatrixUtils.createRealIdentityMatrix(1);
        RealMatrix Q = MatrixUtils.createRealIdentityMatrix(2);

        RiccatiEquationSolver a = new RiccatiEquationSolverImpl(A, B, Q, R);

        RealMatrix P_expected = MatrixUtils.createRealMatrix(new double[][] {
            { 19.7598, -7.6430 }, { -7.6430, 4.7072 }
        });

        assertEquals(f.format(P_expected), f.format(a.getP()));

        // P
        // 19.7598 -7.6430
        // -7.6430 4.7072
        //
        // K
        // -7.6430 4.7072

    }

    @Test
    public void test_imaginary_6_6() {

        RealMatrix A = MatrixUtils.createRealMatrix(new double[][] {
            { 1, 0, 0, 1, 0, 0 }, { 1, 0, 0, 0, 1, 0 },
            { 1, 0, 0, 0, 0, 1 }, { 0, 0, 0, 0, 0, 0 },
            { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 }
        });

        RealMatrix B = MatrixUtils.createRealMatrix(new double[][] {
            { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { -0.0032, 0, 0 },
            { 0, -0.0028, 0 }, { 0, 0, -0.0019 }
        });

        RealMatrix R = MatrixUtils.createRealIdentityMatrix(3);
        RealMatrix Q = MatrixUtils.createRealIdentityMatrix(6);

        RiccatiEquationSolver a = new RiccatiEquationSolverImpl(A, B, Q, R);

        RealMatrix P_expected = MatrixUtils.createRealMatrix(new double[][] {
            { 2.2791, 0.0036, 0.0045, 2.1121, 0.0628, 0.0982 },
            { 0.0036, 0.0002, -0.0000, 0.0017, 0.0029, -0.0011 },
            { 0.0045, -0.0000, 0.0003, 0.0022, -0.0011, 0.0034 },
            { 2.1121, 0.0017, 0.0022, 2.0307, 0.0305, 0.0479 },
            { 0.0628, 0.0029, -0.0011, 0.0305, 0.0746, -0.0387 },
            { 0.0982, -0.0011, 0.0034, 0.0479, -0.0387, 0.0967 } });

        assertEquals(f.format(P_expected),
                     f.format(a.getP().scalarMultiply(1e-05)));

        // 1.0e+05 *
        // 2.2791 0.0036 0.0045 2.1121 0.0628 0.0982
        // 0.0036 0.0002 -0.0000 0.0017 0.0029 -0.0011
        // 0.0045 -0.0000 0.0003 0.0022 -0.0011 0.0034
        // 2.1121 0.0017 0.0022 2.0307 0.0305 0.0479
        // 0.0628 0.0029 -0.0011 0.0305 0.0746 -0.0387
        // 0.0982 -0.0011 0.0034 0.0479 -0.0387 0.0967
    }

    @Test
    public void test_imaginary_6_6_ill_conditioned() {
        // A = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; 0 0 0 0 0 0; 0 0 0 0 0 0;
        // 0 0 0 0 0 0];
        // B = [ 0 0 0; 0 0 0; 0 0 0; -0.0032 0 0; 0 -0.0028 0; 0 0 -0.0019];
        // R= [1 0 0; 0 1 0; 0 0 1];
        // Q = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0;
        // 0 0 0 0 0 1];

        RealMatrix A = MatrixUtils.createRealMatrix(new double[][] {
            { 0, 0, 0, 1, 0, 0 }, { 0, 0, 0, 0, 1, 0 },
            { 0, 0, 0, 0, 0, 1 }, { 0, 0, 0, 0, 0, 0 },
            { 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0 }
        });

        RealMatrix B = MatrixUtils.createRealMatrix(new double[][] {
            { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 }, { -0.0032, 0, 0 },
            { 0, -0.0028, 0 }, { 0, 0, -0.0019 }
        });

        RealMatrix R = MatrixUtils.createRealIdentityMatrix(3);
        RealMatrix Q = MatrixUtils.createRealIdentityMatrix(6);
        RiccatiEquationSolver a = new RiccatiEquationSolverImpl(A, B, Q, R);

        RealMatrix P_expected = MatrixUtils.createRealMatrix(new double[][] {
            { 25.02, -1E-12, 0, 312.5, -1E-12, 0 },
            { -1E-12, 26.7448, 0, 0, 357.1429, 0 },
            { 0, 0, 32.4597, -1E-12, -1E-12, 526.3158 },
            { 312.5, -1E-12, -1E-12, 7818.7475, 0, 0 },
            { -1E-12, 357.1429, -1E-12, 0, 9551.7235, -1E-12 },
            { 0, 0, 526.3158, 0, -1E-12, 17084.0482 }
        });

        assertEquals(f.format(P_expected), f.format(a.getP()));
    }

}
