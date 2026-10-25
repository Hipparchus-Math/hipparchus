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

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

public class QRUpdaterTest {

    /** Tolerance for quantities of order 1 after a few Givens rotations. */
    private static final double TOL = 1.0e-14;

    /** Two independent constraint columns, |A1| = 3. */
    private static final RealVector A1 = new ArrayRealVector(new double[] { 1.0,  2.0,  2.0 });
    private static final RealVector A2 = new ArrayRealVector(new double[] { 0.0,  1.0, -1.0 });

    @Test
    public void testGetRDependsOnActiveCount() {
        final QRUpdater updater = new QRUpdater(MatrixUtils.createRealIdentityMatrix(3));
        Assertions.assertNull(updater.getR());
        Assertions.assertEquals(0, updater.solveRT(new ArrayRealVector(3)).getDimension());

        Assertions.assertTrue(add(updater, A1));
        Assertions.assertEquals(1, updater.getR().getRowDimension());
        // |a1| = 3, so R = [±3]
        Assertions.assertEquals(3.0, FastMath.abs(updater.getR().getEntry(0, 0)), TOL);

        Assertions.assertTrue(add(updater, A2));
        Assertions.assertEquals(2, updater.getR().getRowDimension());

        Assertions.assertTrue(add(updater, 1.0, 0.0, 0.0));
        Assertions.assertEquals(3, updater.getIq());
        Assertions.assertEquals(3, updater.getR().getRowDimension());
    }

    @Test
    public void testFactorizationReproducesActiveColumns() {
        final QRUpdater updater = new QRUpdater(MatrixUtils.createRealIdentityMatrix(3));
        Assertions.assertTrue(add(updater, A1));
        Assertions.assertTrue(add(updater, A2));
        checkFactorization(updater, A1, A2);
    }

    @Test
    public void testDependentConstraintWithoutRotation() {
        // after adding e1, the dependent vector 2 e1 needs no rotation at all:
        // the rollback loop only meets (cc, ss) = (0, 0) markers
        final QRUpdater updater = new QRUpdater(MatrixUtils.createRealIdentityMatrix(3));
        Assertions.assertTrue(add(updater, 1.0, 0.0, 0.0));
        Assertions.assertFalse(add(updater, 2.0, 0.0, 0.0));
        Assertions.assertEquals(1, updater.getIq());
        checkIdentity(updater.getJ());
    }

    @Test
    public void testDependentConstraintRollbackWithZeroCosine() {
        // the last component is non-zero but too small: a rotation with cc = 0 is
        // applied, the constraint is rejected and the rotation must be undone
        final QRUpdater updater = new QRUpdater(MatrixUtils.createRealIdentityMatrix(3));
        Assertions.assertTrue(add(updater, 1.0, 0.0, 0.0));
        Assertions.assertFalse(add(updater, 2.0, 0.0, 1.0e-17));
        Assertions.assertEquals(1, updater.getIq());
        checkIdentity(updater.getJ());
    }

    @Test
    public void testDependentConstraintRollbackWithNonZeroCosine() {
        final QRUpdater updater = new QRUpdater(MatrixUtils.createRealIdentityMatrix(3));
        Assertions.assertTrue(add(updater, 1.0, 0.0, 0.0));
        Assertions.assertFalse(add(updater, 2.0, 1.0e-17, 1.0e-17));
        Assertions.assertEquals(1, updater.getIq());
        checkIdentity(updater.getJ());
    }

    @Test
    public void testDeleteOutOfRangeIsIgnored() {
        final QRUpdater updater = new QRUpdater(MatrixUtils.createRealIdentityMatrix(3));
        Assertions.assertTrue(add(updater, A1));
        final RealMatrix before = updater.getR().copy();
        updater.deleteConstraint(-1);
        updater.deleteConstraint(1);
        Assertions.assertEquals(1, updater.getIq());
        Assertions.assertEquals(0.0, updater.getR().subtract(before).getNorm1(), 0.0);
    }

    @Test
    public void testDeleteKeepsFactorization() {
        final RealVector a3 = new ArrayRealVector(new double[] { 3.0, -1.0,  0.5 });
        final QRUpdater updater = new QRUpdater(MatrixUtils.createRealIdentityMatrix(3));
        Assertions.assertTrue(add(updater, A1));
        Assertions.assertTrue(add(updater, A2));
        Assertions.assertTrue(add(updater, a3));
        updater.deleteConstraint(0);
        checkFactorization(updater, A2, a3);
        updater.deleteConstraint(1);
        checkFactorization(updater, A2);
        updater.deleteConstraint(0);
        Assertions.assertEquals(0, updater.getIq());
        Assertions.assertNull(updater.getR());
    }

    private static boolean add(final QRUpdater updater, final double... a) {
        return add(updater, new ArrayRealVector(a));
    }

    private static boolean add(final QRUpdater updater, final RealVector a) {
        return updater.addConstraint(updater.computeD(a));
    }

    /** Check that J_active * R reproduces the active constraint columns (J orthogonal since L = I). */
    private static void checkFactorization(final QRUpdater updater, final RealVector... active) {
        final int iq = updater.getIq();
        Assertions.assertEquals(active.length, iq);
        final RealMatrix jActive = updater.getJ().getSubMatrix(0, updater.getJ().getRowDimension() - 1, 0, iq - 1);
        final RealMatrix product = jActive.multiply(updater.getR());
        for (int k = 0; k < iq; ++k) {
            Assertions.assertEquals(0.0, product.getColumnVector(k).subtract(active[k]).getNorm(), TOL);
        }
    }

    private static void checkIdentity(final RealMatrix m) {
        Assertions.assertEquals(0.0,
                                m.subtract(MatrixUtils.createRealIdentityMatrix(m.getRowDimension())).getNorm1(),
                                TOL);
    }

}
