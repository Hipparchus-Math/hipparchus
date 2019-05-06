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

import java.util.TreeSet;

import org.hipparchus.complex.Complex;
import org.hipparchus.complex.ComplexComparator;

/**
 * Given a matrix A, it computes an eigen decomposition A = VDV^{T}.
 *
 * It also ensures that eigen values in the diagonal of D are in ascending
 * order.
 *
 */
public class OrderedEigenDecomposition extends EigenDecomposition {

    /**
     * Constructor using the EigenDecomposition as starting point for ordering.
     *
     * @param matrix matrix to decompose
     */
    public OrderedEigenDecomposition(final RealMatrix matrix) {
        super(matrix);

        final RealMatrix D = this.getD();
        final RealMatrix V = this.getV();

        // getting eigen values
        TreeSet<Complex> eigenValues = new TreeSet<>(new ComplexComparator());
        for (int ij = 0; ij < matrix.getRowDimension(); ij++) {
            eigenValues.add(new Complex(getRealEigenvalue(ij),
                                        getImagEigenvalue(ij)));
        }

        // ordering
        for (int ij = 0; ij < matrix.getRowDimension() - 1; ij++) {
            final Complex eigValue = eigenValues.pollFirst();
            int currentIndex = -1;
            // searching the current index
            for (currentIndex = ij; currentIndex < matrix.getRowDimension(); currentIndex++) {
                Complex compCurrent = null;
                if (currentIndex == 0) {
                    compCurrent = new Complex(D.getEntry(currentIndex,
                                                         currentIndex), D.getEntry(currentIndex + 1,
                                                                                   currentIndex));
                } else if (currentIndex + 1 == matrix.getRowDimension()) {
                    compCurrent = new Complex(D.getEntry(currentIndex,
                                                         currentIndex), D.getEntry(currentIndex - 1,
                                                                                   currentIndex));
                } else {
                    if (D.getEntry(currentIndex - 1, currentIndex) != 0) {
                        compCurrent = new Complex(D.getEntry(currentIndex,
                                                             currentIndex), D.getEntry(currentIndex - 1,
                                                                                       currentIndex));
                    } else {
                        compCurrent = new Complex(D.getEntry(currentIndex,
                                                             currentIndex), D.getEntry(currentIndex + 1,
                                                                                       currentIndex));

                    }

                }

                if (eigValue.equals(compCurrent)) {
                    break;
                }
            }

            if (ij == currentIndex) {
                continue;
            }

            // exchanging D
            Complex previousValue = null;
            if (ij == 0) {
                previousValue = new Complex(D.getEntry(ij, ij), D.getEntry(ij + 1, ij));
            } else if (ij + 1 == matrix.getRowDimension()) {
                previousValue = new Complex(D.getEntry(ij, ij), D.getEntry(ij - 1, ij));
            } else {
                if (D.getEntry(ij - 1, ij) != 0) {
                    previousValue = new Complex(D.getEntry(ij, ij), D.getEntry(ij - 1, ij));
                } else {
                    previousValue = new Complex(D.getEntry(ij, ij), D.getEntry(ij + 1, ij));

                }
            }
            // moved eigenvalue
            D.setEntry(ij, ij, eigValue.getReal());
            if (ij == 0) {
                D.setEntry(ij + 1, ij, eigValue.getImaginary());
            } else if ((ij + 1) == matrix.getRowDimension()) {
                D.setEntry(ij - 1, ij, eigValue.getImaginary());
            } else {
                if (eigValue.getImaginary() > 0) {
                    D.setEntry(ij - 1, ij, eigValue.getImaginary());
                    D.setEntry(ij + 1, ij, 0);
                } else {
                    D.setEntry(ij + 1, ij, eigValue.getImaginary());
                    D.setEntry(ij - 1, ij, 0);
                }
            }
            // previous eigen value
            D.setEntry(currentIndex, currentIndex, previousValue.getReal());
            if (currentIndex == 0) {
                D.setEntry(currentIndex + 1, currentIndex,
                           previousValue.getImaginary());
            } else if ((currentIndex + 1) == matrix.getRowDimension()) {
                D.setEntry(currentIndex - 1, currentIndex,
                           previousValue.getImaginary());
            } else {
                if (previousValue.getImaginary() > 0) {
                    D.setEntry(currentIndex - 1, currentIndex,
                               previousValue.getImaginary());
                    D.setEntry(currentIndex + 1, currentIndex, 0);
                } else {
                    D.setEntry(currentIndex + 1, currentIndex,
                               previousValue.getImaginary());
                    D.setEntry(currentIndex - 1, currentIndex, 0);
                }
            }

            // exchanging V
            final double[] previousColumnV = V.getColumn(ij);
            V.setColumn(ij, V.getColumn(currentIndex));
            V.setColumn(currentIndex, previousColumnV);
        }

    }

    /** {@inheritDoc} */
    @Override
    public RealMatrix getVT() {
        return getV().transpose();
    }
}
