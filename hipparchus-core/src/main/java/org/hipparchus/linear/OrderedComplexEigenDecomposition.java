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

import java.util.Arrays;

import org.hipparchus.complex.Complex;

/**
 * Given a matrix A, it computes a complex eigen decomposition A = VDV^{T}.
 *
 * It ensures that eigen values in the diagonal of D are in ascending order.
 *
 */
public class OrderedComplexEigenDecomposition extends ComplexEigenDecomposition {

    /**
     * Constructor for the decomposition.
     *
     * @param matrix real matrix.
     */
    public OrderedComplexEigenDecomposition(final RealMatrix matrix) {
        super(matrix);

        final FieldMatrix<Complex> D = this.getD();
        final FieldMatrix<Complex> V = this.getV();

        // getting eigen values
        IndexedEigenvalue[] eigenValues = new IndexedEigenvalue[D.getRowDimension()];
        for (int ij = 0; ij < matrix.getRowDimension(); ij++) {
            eigenValues[ij] = new IndexedEigenvalue(ij, D.getEntry(ij, ij));
        }

        // ordering
        Arrays.sort(eigenValues);
        for (int ij = 0; ij < matrix.getRowDimension() - 1; ij++) {
            final IndexedEigenvalue eij = eigenValues[ij];

            if (ij == eij.index) {
                continue;
            }

            // exchanging D
            final Complex previousValue = D.getEntry(ij, ij);
            D.setEntry(ij, ij, eij.eigenValue);
            D.setEntry(eij.index, eij.index, previousValue);

            // exchanging V
            for (int k = 0; k  < matrix.getRowDimension(); ++k) {
                final Complex previous = V.getEntry(k, ij);
                V.setEntry(k, ij, V.getEntry(k, eij.index));
                V.setEntry(k, eij.index, previous);
            }

            // exchanging eigenvalue
            for (int k = ij + 1; k < matrix.getRowDimension(); ++k) {
                if (eigenValues[k].index == ij) {
                    eigenValues[k].index = eij.index;
                    break;
                }
            }
        }

        checkDefinition(matrix);
    }

    /** {@inheritDoc} */
    @Override
    public FieldMatrix<Complex> getVT() {
        return getV().transpose();
    }

    /** Container for index and eigenvalue pair. */
    private static class IndexedEigenvalue implements Comparable<IndexedEigenvalue> {

        /** Index in the diagonal matrix. */
        private int index;

        /** Eigenvalue. */
        private final Complex eigenValue;

        /** Build the container from its fields.
         * @param index index in the diagonal matrix
         * @param eigenalue eigenvalue
         */
        IndexedEigenvalue(final int index, final Complex eigenvalue) {
            this.index      = index;
            this.eigenValue = eigenvalue;
        }

        /** {@inheritDoc}
         * <p>
         * Ordering uses real ordering as the primary sort order and
         * imaginary ordering as the secondary sort order.
         * </p>
         */
        @Override
        public int compareTo(final IndexedEigenvalue other) {
            final int cR = Double.compare(eigenValue.getReal(), other.eigenValue.getReal());
            if (cR == 0) {
                return Double.compare(eigenValue.getImaginary(),other.eigenValue.getImaginary());
            } else {
                return cR;
            }
        }

    }

}
