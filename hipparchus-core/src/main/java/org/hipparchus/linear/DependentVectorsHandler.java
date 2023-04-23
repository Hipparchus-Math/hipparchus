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

import java.util.List;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;

/** Enumerate to specify how dependent vectors should be handled in
 * {@link MatrixUtils#orthonormalize(List, double, DependentVectorsHandler)} and
 * {@link MatrixUtils#orthonormalize(Field, List, CalculusFieldElement, DependentVectorsHandler)}.
 * @since 2.1
 */
public enum DependentVectorsHandler {

    /** Generate a {@link MathIllegalArgumentException} if dependent vectors are found. */
    GENERATE_EXCEPTION {

        /** {@inheritDoc} */
        @Override
        public int manageDependent(final int index, final List<RealVector> basis) {
            // generate exception, dependent vectors are forbidden with this settings
            throw new MathIllegalArgumentException(LocalizedCoreFormats.ZERO_NORM);
        }


        /** {@inheritDoc} */
        @Override
        public <T extends CalculusFieldElement<T>> int manageDependent(final Field<T> field,
                                                                       final int index,
                                                                       final List<FieldVector<T>> basis) {
            // generate exception, dependent vectors are forbidden with this settings
            throw new MathIllegalArgumentException(LocalizedCoreFormats.ZERO_NORM);
        }

    },

    /** Replace dependent vectors by vectors with norm 0.
     * <p>
     * This behavior matches the Wolfram language API. It keeps the
     * number of output vectors equal to the number of input vectors.
     * The only two norms output vectors can have are 0 and 1.
     * </p>
     */
    ADD_ZERO_VECTOR {

        /** {@inheritDoc} */
        @Override
        public int manageDependent(final int index, final List<RealVector> basis) {
            // add a zero vector, preserving output vector size (and dropping its normalization property)
            basis.set(index, MatrixUtils.createRealVector(basis.get(index).getDimension()));
            return index + 1;
        }


        /** {@inheritDoc} */
        @Override
        public <T extends CalculusFieldElement<T>> int manageDependent(final Field<T> field,
                                                                       final int index,
                                                                       final List<FieldVector<T>> basis) {
            // add a zero vector, preserving output vector size (and dropping its normalization property)
            basis.set(index, MatrixUtils.createFieldVector(field, basis.get(index).getDimension()));
            return index + 1;
        }

    },

    /** Ignore dependent vectors.
     * <p>
     * This behavior ensures the output vectors form an orthonormal
     * basis, i.e. all vectors are independent and they all have norm 1.
     * The number of output vectors may be smaller than the number of
     * input vectors, this number corresponds to the dimension of the
     * span of the input vectors.
     * </p>
     */
    REDUCE_BASE_TO_SPAN {

        /** {@inheritDoc} */
        @Override
        public int manageDependent(final int index, final List<RealVector> basis) {
            // remove dependent vector
            basis.remove(index);
            return index;
        }


        /** {@inheritDoc} */
        @Override
        public <T extends CalculusFieldElement<T>> int manageDependent(final Field<T> field,
                                                                       final int index,
                                                                       final List<FieldVector<T>> basis) {
            // remove dependent vector
            basis.remove(index);
            return index;
        }

    };

    /** Manage a dependent vector.
     * @param index of the vector in the basis
     * @param basis placeholder for basis vectors
     * @return next index to manage
     */
    public abstract int manageDependent(int index, List<RealVector> basis);

    /** Manage a dependent vector.
     * @param <T> type of the vectors components
     * @param field field to which the vectors belong
     * @param index of the vector in the basis
     * @param basis placeholder for basis vectors
     * @return next index to manage
     */
    public abstract <T extends CalculusFieldElement<T>> int manageDependent(Field<T> field, int index, List<FieldVector<T>> basis);

}
