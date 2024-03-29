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
package org.hipparchus.analysis.differentiation;

import java.io.Serializable;

import org.hipparchus.Field;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;

/** Factory for {@link DerivativeStructure}.
 * <p>This class is a factory for {@link DerivativeStructure} instances.</p>
 * <p>Instances of this class are guaranteed to be immutable.</p>
 * @see DerivativeStructure
 * @since 1.1
 */
public class DSFactory implements Serializable {

    /** Serializable UID. */
    private static final long serialVersionUID = 20161222L;

    /** Compiler for the current dimensions. */
    private final transient DSCompiler compiler;

    /** Field the {@link DerivativeStructure} instances belong to. */
    private final transient DSField derivativeField;

    /** Simple constructor.
     * @param parameters number of free parameters
     * @param order derivation order
     */
    public DSFactory(final int parameters, final int order) {
        this.compiler        = DSCompiler.getCompiler(parameters, order);
        this.derivativeField = new DSField(constant(0.0), constant(1.0), constant(FastMath.PI));
    }

    /** Get the {@link Field} the {@link DerivativeStructure} instances belong to.
     * @return {@link Field} the {@link DerivativeStructure} instances belong to
     */
    public DSField getDerivativeField() {
        return derivativeField;
    }

    /** Build a {@link DerivativeStructure} representing a constant value.
     * @param value value of the constant
     * @return a {@link DerivativeStructure} representing a constant value
     */
    public DerivativeStructure constant(double value) {
        final DerivativeStructure ds = new DerivativeStructure(this);
        ds.setDerivativeComponent(0, value);
        return ds;
    }

    /** Build a {@link DerivativeStructure} representing a variable.
     * <p>Instances built using this method are considered
     * to be the free variables with respect to which differentials
     * are computed. As such, their differential with respect to
     * themselves is +1.</p>
     * @param index index of the variable (from 0 to
     * {@link #getCompiler()}.{@link DSCompiler#getFreeParameters() getFreeParameters()} - 1)
     * @param value value of the variable
     * @exception MathIllegalArgumentException if index if greater or
     * equal to {@link #getCompiler()}.{@link DSCompiler#getFreeParameters() getFreeParameters()}.
     * @return a {@link DerivativeStructure} representing a variable
     */
    public DerivativeStructure variable(final int index, final double value)
        throws MathIllegalArgumentException {

        if (index >= getCompiler().getFreeParameters()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_LARGE_BOUND_EXCLUDED,
                                                   index, getCompiler().getFreeParameters());
        }

        final DerivativeStructure ds = new DerivativeStructure(this);
        ds.setDerivativeComponent(0, value);

        if (getCompiler().getOrder() > 0) {
            // the derivative of the variable with respect to itself is 1.
            ds.setDerivativeComponent(DSCompiler.getCompiler(index, getCompiler().getOrder()).getSize(), 1.0);
        }

        return ds;

    }

    /** Build a {@link DerivativeStructure} from all its derivatives.
     * @param derivatives derivatives sorted according to
     * {@link DSCompiler#getPartialDerivativeIndex(int...)}
     * @return a {@link DerivativeStructure} with specified derivatives
     * @exception MathIllegalArgumentException if derivatives array does not match the
     * {@link DSCompiler#getSize() size} expected by the compiler
     * @exception MathIllegalArgumentException if order is too large
     * @see DerivativeStructure#getAllDerivatives()
     */
    @SafeVarargs
    public final DerivativeStructure build(final double ... derivatives)
        throws MathIllegalArgumentException {

        if (derivatives.length != compiler.getSize()) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   derivatives.length, compiler.getSize());
        }

        return new DerivativeStructure(this, derivatives);

    }

    /** Build a {@link DerivativeStructure} with an uninitialized array.
     * <p>This method is intended only for DerivativeStructure internal use.</p>
     * @return a {@link DerivativeStructure} with an uninitialized array
     */
    DerivativeStructure build() {
        return new DerivativeStructure(this);
    }

    /** Get the compiler for the current dimensions.
     * @return compiler for the current dimensions
     */
    public DSCompiler getCompiler() {
        return compiler;
    }

    /** Check rules set compatibility.
     * @param factory other factory field to check against instance
     * @exception MathIllegalArgumentException if number of free parameters or orders are inconsistent
     */
    void checkCompatibility(final DSFactory factory) throws MathIllegalArgumentException {
        compiler.checkCompatibility(factory.compiler);
    }

    /**
     * Replace the instance with a data transfer object for serialization.
     * @return data transfer object that will be serialized
     */
    private Object writeReplace() {
        return new DataTransferObject(compiler.getFreeParameters(), compiler.getOrder());
    }

    /** Internal class used only for serialization. */
    private static class DataTransferObject implements Serializable {

        /** Serializable UID. */
        private static final long serialVersionUID = 20161222L;

        /** Number of variables.
         * @serial
         */
        private final int variables;

        /** Derivation order.
         * @serial
         */
        private final int order;

        /** Simple constructor.
         * @param variables number of variables
         * @param order derivation order
         */
        DataTransferObject(final int variables, final int order) {
            this.variables = variables;
            this.order     = order;
        }

        /** Replace the deserialized data transfer object with a {@link DSFactory}.
         * @return replacement {@link DSFactory}
         */
        private Object readResolve() {
            return new DSFactory(variables, order);
        }

    }

    /** Field for {link DerivativeStructure} instances.
     */
    public static class DSField implements Field<DerivativeStructure> {

        /** Constant function evaluating to 0.0. */
        private final DerivativeStructure zero;

        /** Constant function evaluating to 1.0. */
        private final DerivativeStructure one;

        /** Constant function evaluating to π. */
        private final DerivativeStructure pi;

        /** Simple constructor.
         * @param zero constant function evaluating to 0.0
         * @param one constant function evaluating to 1.0
         * @param pi constant function evaluating to π
         */
        DSField(final DerivativeStructure zero, final DerivativeStructure one, final DerivativeStructure pi) {
            this.zero = zero;
            this.one  = one;
            this.pi   = pi;
        }

        /** {@inheritDoc} */
        @Override
        public DerivativeStructure getZero() {
            return zero;
        }

        /** {@inheritDoc} */
        @Override
        public DerivativeStructure getOne() {
            return one;
        }

        /** Get the Archimedes constant π.
         * <p>
         * Archimedes constant is the ratio of a circle's circumference to its diameter.
         * </p>
         * @return Archimedes constant π
         * @since 2.0
         */
        public DerivativeStructure getPi() {
            return pi;
        }

        /** {@inheritDoc} */
        @Override
        public Class<DerivativeStructure> getRuntimeClass() {
            return DerivativeStructure.class;
        }

        /** {@inheritDoc} */
        @Override
        public boolean equals(final Object other) {
            if (this == other) {
                return true;
            } else if (other instanceof DSField) {
                DSFactory lhsFactory = zero.getFactory();
                DSFactory rhsFactory = ((DSField) other).zero.getFactory();
                return lhsFactory.compiler == rhsFactory.compiler;
            } else {
                return false;
            }
        }

        /** {@inheritDoc} */
        @Override
        public int hashCode() {
            final DSCompiler compiler = zero.getFactory().getCompiler();
            return 0x9943b886 ^ (compiler.getFreeParameters() << 16 & compiler.getOrder());
        }

    }

}
