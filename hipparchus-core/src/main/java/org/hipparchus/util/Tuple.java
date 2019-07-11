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
package org.hipparchus.util;

import java.util.Arrays;

import org.hipparchus.Field;
import org.hipparchus.FieldElement;
import org.hipparchus.RealFieldElement;
import org.hipparchus.exception.MathIllegalArgumentException;

/**
 * This class allows to perform the same computation of all components of a Tuple at once.
 * @since 1.2
 */
public class Tuple implements RealFieldElement<Tuple> {

    /** Components of the tuple. */
    private final double[] values;

    /** Field the instance belongs to. */
    private final transient TupleField field;

    /** Creates a new instance from its components.
     * @param x components of the tuple
     */
    public Tuple(final double... x) {
        this(new TupleField(x.length), x.clone());
    }

    /** Creates a new instance from its components.
     * @param field field the instance belongs to
     * @param x components of the tuple (beware, it is <em>not</em> copied, it is shared with caller)
     */
    private Tuple(final TupleField field, final double[] x) {// NOPMD - storing user-supplied array is intentional and documented here
        this.values = x;
        this.field  = field;
    }

    /** Get the dimension of the tuple.
     * @return dimension of the tuple
     */
    public int getDimension() {
        return values.length;
    }

    /** Get one component of the tuple.
     * @param index index of the component, between 0 and {@link #getDimension() getDimension()} - 1
     * @return value of the component
     */
    public double getComponent(final int index) {
        return values[index];
    }

    /** Get all components of the tuple.
     * @return all components
     */
    public double[] getComponents() {
        return values.clone();
    }

    /** {@inheritDoc} */
    @Override
    public Field<Tuple> getField() {
        return field;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple add(final Tuple a) {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i] + a.values[i];
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple subtract(final Tuple a) {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i] - a.values[i];
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple negate() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = -values[i];
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple multiply(final Tuple a) {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i] * a.values[i];
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple multiply(final int n) {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i] * n;
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple divide(final Tuple a) {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i] / a.values[i];
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple reciprocal() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = 1.0 / values[i];
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public boolean equals(final Object obj) {
        if (obj instanceof Tuple) {
            final Tuple that = (Tuple) obj;
            if (getDimension() == that.getDimension()) {
                boolean equals = true;
                for (int i = 0; i < values.length; ++i) {
                    equals &= Double.doubleToRawLongBits(values[i]) == Double.doubleToRawLongBits(that.values[i]);
                }
                return equals;
            }
        }
        return false;
    }

    /** {@inheritDoc} */
    @Override
    public int hashCode() {
        return  0x34b1a444 + Arrays.hashCode(values);
    }

    /** {@inheritDoc} */
    @Override
    public double getReal() {
        return values[0];
    }

    /** {@inheritDoc} */
    @Override
    public Tuple add(final double a) {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i] + a;
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple subtract(final double a) {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i] - a;
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple multiply(final double a) {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i] * a;
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple divide(final double a) {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i] / a;
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple remainder(final double a) {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.IEEEremainder(values[i], a);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple remainder(final Tuple a) {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.IEEEremainder(values[i], a.values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple abs() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.abs(values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple ceil() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.ceil(values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple floor() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.floor(values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple rint() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.rint(values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public long round() {
        return FastMath.round(values[0]);
    }

    /** {@inheritDoc} */
    @Override
    public Tuple signum() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.signum(values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple copySign(final Tuple sign) {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.copySign(values[i], sign.values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple copySign(final double sign) {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.copySign(values[i], sign);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple scalb(final int n) {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.scalb(values[i], n);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple hypot(final Tuple y) {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.hypot(values[i], y.values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple sqrt() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.sqrt(values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple cbrt() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.cbrt(values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple rootN(final int n) {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            if (values[i] < 0) {
                result.values[i] = -FastMath.pow(-values[i], 1.0 / n);
            } else {
                result.values[i] = FastMath.pow(values[i], 1.0 / n);
            }
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple pow(final double p) {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.pow(values[i], p);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple pow(final int n) {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.pow(values[i], n);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple pow(final Tuple e) {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.pow(values[i], e.values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple exp() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.exp(values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple expm1() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.expm1(values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple log() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.log(values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple log1p() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.log1p(values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple log10() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.log10(values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple cos() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.cos(values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple sin() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.sin(values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldSinCos<Tuple> sinCos() {
        final Tuple sin = new Tuple(field, new double[values.length]);
        final Tuple cos = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            final SinCos sc = FastMath.sinCos(values[i]);
            sin.values[i] = sc.sin();
            cos.values[i] = sc.cos();
        }
        return new FieldSinCos<>(sin, cos);
    }

    /** {@inheritDoc} */
    @Override
    public Tuple tan() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.tan(values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple acos() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.acos(values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple asin() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.asin(values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple atan() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.atan(values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple atan2(final Tuple x) {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.atan2(values[i], x.values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple cosh() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.cosh(values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple sinh() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.sinh(values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple tanh() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.tanh(values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple acosh() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.acosh(values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple asinh() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.asinh(values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple atanh() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.atanh(values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple toDegrees() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.toDegrees(values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple toRadians() {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = FastMath.toRadians(values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple linearCombination(final Tuple[] a, final Tuple[] b)
        throws MathIllegalArgumentException {
        final Tuple result = new Tuple(field, new double[values.length]);
        MathUtils.checkDimension(a.length, b.length);
        final double[] aDouble = new double[a.length];
        final double[] bDouble = new double[b.length];
        for (int i = 0; i < values.length; ++i) {
            for (int j = 0; j < a.length; ++j) {
                aDouble[j] = a[j].values[i];
                bDouble[j] = b[j].values[i];
            }
            result.values[i] = MathArrays.linearCombination(aDouble, bDouble);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple linearCombination(final double[] a, final Tuple[] b)
        throws MathIllegalArgumentException {
        final Tuple result = new Tuple(field, new double[values.length]);
        MathUtils.checkDimension(a.length, b.length);
        final double[] bDouble = new double[b.length];
        for (int i = 0; i < values.length; ++i) {
            for (int j = 0; j < a.length; ++j) {
                bDouble[j] = b[j].values[i];
            }
            result.values[i] = MathArrays.linearCombination(a, bDouble);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple linearCombination(final Tuple a1, final Tuple b1,
                                   final Tuple a2, final Tuple b2) {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = MathArrays.linearCombination(a1.values[i], b1.values[i],
                                                            a2.values[i], b2.values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple linearCombination(final double a1, final Tuple b1,
                                   final double a2, final Tuple b2) {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = MathArrays.linearCombination(a1, b1.values[i],
                                                            a2, b2.values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple linearCombination(final Tuple a1, final Tuple b1,
                                   final Tuple a2, final Tuple b2,
                                   final Tuple a3, final Tuple b3) {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = MathArrays.linearCombination(a1.values[i], b1.values[i],
                                                            a2.values[i], b2.values[i],
                                                            a3.values[i], b3.values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple linearCombination(final double a1, final Tuple b1,
                                   final double a2, final Tuple b2,
                                   final double a3, final Tuple b3) {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = MathArrays.linearCombination(a1, b1.values[i],
                                                            a2, b2.values[i],
                                                            a3, b3.values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple linearCombination(final Tuple a1, final Tuple b1,
                                   final Tuple a2, final Tuple b2,
                                   final Tuple a3, final Tuple b3,
                                   final Tuple a4, final Tuple b4) {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = MathArrays.linearCombination(a1.values[i], b1.values[i],
                                                            a2.values[i], b2.values[i],
                                                            a3.values[i], b3.values[i],
                                                            a4.values[i], b4.values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public Tuple linearCombination(final double a1, final Tuple b1,
                                   final double a2, final Tuple b2,
                                   final double a3, final Tuple b3,
                                   final double a4, final Tuple b4) {
        final Tuple result = new Tuple(field, new double[values.length]);
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = MathArrays.linearCombination(a1, b1.values[i],
                                                            a2, b2.values[i],
                                                            a3, b3.values[i],
                                                            a4, b4.values[i]);
        }
        return result;
    }

    /** Field for {link Tuple} instances.
     */
    private static class TupleField implements Field<Tuple> {

        /** Constant function evaluating to 0.0. */
        private final Tuple zero;

        /** Constant function evaluating to 1.0. */
        private final Tuple one;

        /** Simple constructor.
         * @param dimension dimension of the tuple
         */
        TupleField(final int dimension) {
            final double[] zeroData = new double[dimension];
            final double[] oneData  = new double[dimension];
            Arrays.fill(oneData, 1.0);
            this.zero = new Tuple(this, zeroData);
            this.one  = new Tuple(this, oneData);
        }

        /** {@inheritDoc} */
        @Override
        public Tuple getZero() {
            return zero;
        }

        /** {@inheritDoc} */
        @Override
        public Tuple getOne() {
            return one;
        }

        /** {@inheritDoc} */
        @Override
        public Class<? extends FieldElement<Tuple>> getRuntimeClass() {
            return Tuple.class;
        }

        /** {@inheritDoc} */
        @Override
        public boolean equals(final Object other) {
            if (other instanceof TupleField) {
                return zero.getDimension() == ((TupleField) other).zero.getDimension();
            } else {
                return false;
            }
        }

        /** {@inheritDoc} */
        @Override
        public int hashCode() {
            return 0x6672493d ^ zero.getDimension();
        }

    }

}
