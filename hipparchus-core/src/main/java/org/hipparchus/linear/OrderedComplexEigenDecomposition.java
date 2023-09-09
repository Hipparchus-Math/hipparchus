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

import java.util.Arrays;
import java.util.Comparator;

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
        this(matrix,
             ComplexEigenDecomposition.DEFAULT_EIGENVECTORS_EQUALITY,
             ComplexEigenDecomposition.DEFAULT_EPSILON,
             ComplexEigenDecomposition.DEFAULT_EPSILON_AV_VD_CHECK);
    }

    /**
     * Constructor for decomposition.
     * <p>
     * The {@code eigenVectorsEquality} threshold is used to ensure the L∞-normalized
     * eigenvectors found using inverse iteration are different from each other.
     * if \(min(|e_i-e_j|,|e_i+e_j|)\) is smaller than this threshold, the algorithm
     * considers it has found again an already known vector, so it drops it and attempts
     * a new inverse iteration with a different start vector. This value should be
     * much larger than {@code epsilon} which is used for convergence
     * </p>
     * <p>
     * This constructor calls {@link #OrderedComplexEigenDecomposition(RealMatrix, double,
     * double, double, Comparator)} with a comparator using real ordering as the primary
     * sort order and imaginary ordering as the secondary sort order..
     * </p>
     * @param matrix real matrix.
     * @param eigenVectorsEquality threshold below which eigenvectors are considered equal
     * @param epsilon Epsilon used for internal tests (e.g. is singular, eigenvalue ratio, etc.)
     * @param epsilonAVVDCheck Epsilon criteria for final AV=VD check
     * @since 1.9
     */
    public OrderedComplexEigenDecomposition(final RealMatrix matrix, final double eigenVectorsEquality,
                                            final double epsilon, final double epsilonAVVDCheck) {
        this(matrix, eigenVectorsEquality, epsilon, epsilonAVVDCheck,
             (c1, c2) -> {
                 final int cR = Double.compare(c1.getReal(), c2.getReal());
                 if (cR == 0) {
                     return Double.compare(c1.getImaginary(), c2.getImaginary());
                 } else {
                     return cR;
                 }
             });
    }

    /**
     * Constructor for decomposition.
     * <p>
     * The {@code eigenVectorsEquality} threshold is used to ensure the L∞-normalized
     * eigenvectors found using inverse iteration are different from each other.
     * if \(min(|e_i-e_j|,|e_i+e_j|)\) is smaller than this threshold, the algorithm
     * considers it has found again an already known vector, so it drops it and attempts
     * a new inverse iteration with a different start vector. This value should be
     * much larger than {@code epsilon} which is used for convergence
     * </p>
     * @param matrix real matrix.
     * @param eigenVectorsEquality threshold below which eigenvectors are considered equal
     * @param epsilon Epsilon used for internal tests (e.g. is singular, eigenvalue ratio, etc.)
     * @param epsilonAVVDCheck Epsilon criteria for final AV=VD check
     * @param eigenValuesComparator comparator for sorting eigen values
     * @since 3.0
     */
    public OrderedComplexEigenDecomposition(final RealMatrix matrix, final double eigenVectorsEquality,
                                            final double epsilon, final double epsilonAVVDCheck,
                                            final Comparator<Complex> eigenValuesComparator) {
        super(matrix, eigenVectorsEquality, epsilon, epsilonAVVDCheck);
        final FieldMatrix<Complex> D = this.getD();
        final FieldMatrix<Complex> V = this.getV();

        // getting eigen values
        IndexedEigenvalue[] eigenValues = new IndexedEigenvalue[D.getRowDimension()];
        for (int ij = 0; ij < matrix.getRowDimension(); ij++) {
            eigenValues[ij] = new IndexedEigenvalue(ij, D.getEntry(ij, ij));
        }

        // ordering
        Arrays.sort(eigenValues, (v1, v2) -> eigenValuesComparator.compare(v1.eigenValue, v2.eigenValue));
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

        // reorder the eigenvalues and eigenvector s array in base class
        matricesToEigenArrays();

        checkDefinition(matrix);

    }

    /** {@inheritDoc} */
    @Override
    public FieldMatrix<Complex> getVT() {
        return getV().transpose();
    }

    /** Container for index and eigenvalue pair. */
    private static class IndexedEigenvalue {

        /** Index in the diagonal matrix. */
        private int index;

        /** Eigenvalue. */
        private final Complex eigenValue;

        /** Build the container from its fields.
         * @param index index in the diagonal matrix
         * @param eigenvalue eigenvalue
         */
        IndexedEigenvalue(final int index, final Complex eigenvalue) {
            this.index      = index;
            this.eigenValue = eigenvalue;
        }

        /** {@inheritDoc} */
        @Override
        public boolean equals(final Object other) {

            if (this == other) {
                return true;
            }

            if (other instanceof IndexedEigenvalue) {
                final IndexedEigenvalue rhs = (IndexedEigenvalue) other;
                return eigenValue.equals(rhs.eigenValue);
            }

            return false;

        }

        /**
         * Get a hashCode for the pair.
         * @return a hash code value for this object
         */
        @Override
        public int hashCode() {
            return 4563 + index + eigenValue.hashCode();
        }

    }

}
