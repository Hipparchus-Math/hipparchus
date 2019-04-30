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
        TreeSet<Complex> eigenValues = new TreeSet<>(new ComplexComparator());
        for (int ij = 0; ij < matrix.getRowDimension(); ij++) {
            eigenValues.add(D.getEntry(ij, ij));
        }

        // ordering
        for (int ij = 0; ij < matrix.getRowDimension() - 1; ij++) {
            final Complex eigValue = eigenValues.pollFirst();
            int currentIndex = -1;
            // searching the current index
            for (currentIndex = ij; currentIndex < matrix.getRowDimension(); currentIndex++) {
                Complex compCurrent = D.getEntry(currentIndex, currentIndex);
                if (eigValue.equals(compCurrent)) {
                    break;
                }
            }

            if (ij == currentIndex) {
                continue;
            }

            // exchanging D
            Complex previousValue = D.getEntry(ij, ij);
            D.setEntry(ij, ij, eigValue);
            D.setEntry(currentIndex, currentIndex, previousValue);

            // exchanging V
            final Complex[] previousColumnV = V.getColumn(ij);
            V.setColumn(ij, V.getColumn(currentIndex));
            V.setColumn(currentIndex, previousColumnV);
        }

        checkDefinition(matrix);
    }

    /** {@inheritDoc} */
    @Override
    public FieldMatrix<Complex> getVT() {
        return getV().transpose();
    }
}
