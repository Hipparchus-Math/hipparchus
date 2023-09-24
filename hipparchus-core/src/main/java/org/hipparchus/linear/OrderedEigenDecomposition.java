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

import java.util.TreeSet;
import java.util.Arrays;
import java.util.Comparator;

import org.hipparchus.complex.Complex;
import org.hipparchus.complex.ComplexComparator;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.MathRuntimeException;

/**
 * Given a matrix A, it computes an eigen decomposition A = VDV^{T}.
 *
 * It also ensures that eigen values in the diagonal of D are in ascending
 * order.
 *
 */
public class OrderedEigenDecomposition extends EigenDecompositionNonSymmetric {

    /**
     * Constructor using the EigenDecomposition as starting point for ordering.
     *
     * @param matrix matrix to decompose
     */
    public OrderedEigenDecomposition(final RealMatrix matrix) {
        this(matrix, DEFAULT_EPSILON, new ComplexComparator());
    }

    /**
     * Calculates the eigen decomposition of the given real matrix.
     * <p>
     * Supports decomposition of a general matrix since 3.1.
     *
     * @param matrix Matrix to decompose.
     * @param epsilon Epsilon used for internal tests (e.g. is singular, eigenvalue ratio, etc.)
     * @throws MathIllegalStateException if the algorithm fails to converge.
     * @param eigenValuesComparator comparator for sorting eigen values
     * @throws MathRuntimeException if the decomposition of a general matrix
     * results in a matrix with zero norm
     * @since 3.0
     */
    public OrderedEigenDecomposition(final RealMatrix matrix, final double epsilon,
                                     final Comparator<Complex> eigenValuesComparator) {
        super(matrix, epsilon);

        final RealMatrix D = this.getD();
        final RealMatrix V = this.getV();
        final RealMatrix newD = new Array2DRowRealMatrix(D.getData(), true);
        final RealMatrix newV = new Array2DRowRealMatrix(V.getData(), true);

        // getting eigen values
        TreeSet<Complex> eigenValues = new TreeSet<>(eigenValuesComparator);
        IndexedEigenvalue[] newEigenValues = new IndexedEigenvalue[newD.getRowDimension()];
        for (int ij = 0; ij < matrix.getRowDimension(); ij++) {
            newEigenValues[ij] = new IndexedEigenvalue(ij, getEigenvalue(ij));
            eigenValues.add(getEigenvalue(ij));
        }

        // ordering
        Arrays.sort(newEigenValues, (v1, v2) -> eigenValuesComparator.compare(v1.getEigenvalue(), v2.getEigenvalue()));
        for (int ij = 0; ij < matrix.getRowDimension() - 1; ij++) {
            final IndexedEigenvalue eij = newEigenValues[ij];
            final Complex eigValue = eigenValues.pollFirst();
            int currentIndex;
            // searching the current index
            for (currentIndex = ij; currentIndex < matrix.getRowDimension(); currentIndex++) {
                Complex compCurrent;
                if (currentIndex == 0) {
                    compCurrent = new Complex(D.getEntry(currentIndex,     currentIndex),
                                              D.getEntry(currentIndex + 1, currentIndex));
                } else if (currentIndex + 1 == matrix.getRowDimension()) {
                    compCurrent = new Complex(D.getEntry(currentIndex,     currentIndex),
                                              D.getEntry(currentIndex - 1, currentIndex));
                } else {
                    if (D.getEntry(currentIndex - 1, currentIndex) != 0) {
                        compCurrent = new Complex(D.getEntry(currentIndex,     currentIndex),
                                                  D.getEntry(currentIndex - 1, currentIndex));
                    } else {
                        compCurrent = new Complex(D.getEntry(currentIndex,     currentIndex),
                                                  D.getEntry(currentIndex + 1, currentIndex));

                    }

                }

                if (eigValue.equals(compCurrent)) {
                    break;
                }
            }

            if (ij == eij.getIndex()) {
                if (ij == currentIndex) {
                    continue;
                }
            }

            // exchanging D
            Complex previousValue;
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
            Complex newPreviousValue;
            if (ij == 0) {
                newPreviousValue = new Complex(newD.getEntry(ij, ij), newD.getEntry(ij + 1, ij));
            } else if (ij + 1 == matrix.getRowDimension()) {
                newPreviousValue = new Complex(newD.getEntry(ij, ij), newD.getEntry(ij - 1, ij));
            } else {
                if (newD.getEntry(ij - 1, ij) != 0) {
                    newPreviousValue = new Complex(newD.getEntry(ij, ij), newD.getEntry(ij - 1, ij));
                } else {
                    newPreviousValue = new Complex(newD.getEntry(ij, ij), newD.getEntry(ij + 1, ij));

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
            newD.setEntry(ij, ij, eij.getEigenvalue().getReal());
            if (ij == 0) {
                newD.setEntry(ij + 1, ij, eij.getEigenvalue().getImaginary());
            } else if ((ij + 1) == matrix.getRowDimension()) {
                newD.setEntry(ij - 1, ij, eij.getEigenvalue().getImaginary());
            } else {
                if (eij.getEigenvalue().getImaginary() > 0) {
                    newD.setEntry(ij - 1, ij, eij.getEigenvalue().getImaginary());
                    newD.setEntry(ij + 1, ij, 0);
                } else {
                    newD.setEntry(ij + 1, ij, eij.getEigenvalue().getImaginary());
                    newD.setEntry(ij - 1, ij, 0);
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
            newD.setEntry(eij.getIndex(), eij.getIndex(), newPreviousValue.getReal());
            if (eij.getIndex() == 0) {
                newD.setEntry(eij.getIndex() + 1, eij.getIndex(),
                              newPreviousValue.getImaginary());
            } else if ((eij.getIndex() + 1) == matrix.getRowDimension()) {
                newD.setEntry(eij.getIndex() - 1, eij.getIndex(),
                              newPreviousValue.getImaginary());
            } else {
                if (newPreviousValue.getImaginary() > 0) {
                    newD.setEntry(eij.getIndex() - 1, eij.getIndex(),
                                  newPreviousValue.getImaginary());
                    newD.setEntry(eij.getIndex() + 1, eij.getIndex(), 0);
                } else {
                    newD.setEntry(eij.getIndex() + 1, eij.getIndex(),
                                  newPreviousValue.getImaginary());
                    newD.setEntry(eij.getIndex() - 1, eij.getIndex(), 0);
                }
            }

            // exchanging V
            final double[] previousColumnV = V.getColumn(ij);
            V.setColumn(ij, V.getColumn(currentIndex));
            V.setColumn(currentIndex, previousColumnV);
            final double[] newPreviousColumnV = newV.getColumn(ij);
            newV.setColumn(ij, newV.getColumn(eij.getIndex()));
            newV.setColumn(eij.getIndex(), newPreviousColumnV);

//            // exchanging eigenvalue
//            for (int k = ij + 1; k < matrix.getRowDimension(); ++k) {
//                if (eigenValues[k].getIndex() == ij) {
//                    eigenValues[k].setIndex(eij.getIndex());
//                    break;
//                }
//            }

        }

    }

}
