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
 * @param <T> the type of the field elements
 * @since 1.2
 */
public class FieldTuple<T extends RealFieldElement<T>> implements RealFieldElement<FieldTuple<T>> {

    /** Components of the tuple. */
    private final T[] values;

    /** Field the instance belongs to. */
    private final transient FieldTupleField<T> field;

    /** Creates a new instance from its components.
     * @param x components of the tuple
     */
    @SafeVarargs
    public FieldTuple(final T... x) {
        this(new FieldTupleField<>(x[0].getField(), x.length), x.clone());
    }

    /** Creates a new instance from its components.
     * @param field field the instance belongs to
     * @param x components of the tuple (beware, it is <em>not</em> copied, it is shared with caller)
     */
    private FieldTuple(final FieldTupleField<T> field, final T[] x) { // NOPMD - storing user-supplied array is intentional and documented here
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
    public T getComponent(final int index) {
        return values[index];
    }

    /** Get all components of the tuple.
     * @return all components
     */
    public T[] getComponents() {
        return values.clone();
    }

    /** {@inheritDoc} */
    @Override
    public Field<FieldTuple<T>> getField() {
        return field;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> add(final FieldTuple<T> a) {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].add(a.values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> subtract(final FieldTuple<T> a) {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].subtract(a.values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> negate() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].negate();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> multiply(final FieldTuple<T> a) {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].multiply(a.values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> multiply(final int n) {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].multiply(n);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> divide(final FieldTuple<T> a) {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].divide(a.values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> reciprocal() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].reciprocal();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public boolean equals(final Object obj) {
        if (obj instanceof FieldTuple<?> ) {
            @SuppressWarnings("unchecked")
            final FieldTuple<T> that = (FieldTuple<T>) obj;
            if (getDimension() == that.getDimension()) {
                boolean equals = true;
                for (int i = 0; i < values.length; ++i) {
                    equals &= values[i].equals(that.values[i]);
                }
                return equals;
            }
        }
        return false;
    }

    /** {@inheritDoc} */
    @Override
    public int hashCode() {
        return  0x58f61de5 + Arrays.hashCode(values);
    }

    /** {@inheritDoc} */
    @Override
    public double getReal() {
        return values[0].getReal();
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> add(final double a) {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].add(a);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> subtract(final double a) {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].subtract(a);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> multiply(final double a) {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].multiply(a);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> divide(final double a) {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].divide(a);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> remainder(final double a) {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].remainder(a);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> remainder(final FieldTuple<T> a) {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].remainder(a.values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> abs() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].abs();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> ceil() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].ceil();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> floor() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].floor();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> rint() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].rint();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public long round() {
        return values[0].round();
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> signum() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].signum();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> copySign(final FieldTuple<T> sign) {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].copySign(sign.values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> copySign(final double sign) {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].copySign(sign);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> scalb(final int n) {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].scalb(n);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> hypot(final FieldTuple<T> y) {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].hypot(y.values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> sqrt() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].sqrt();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> cbrt() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].cbrt();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> rootN(final int n) {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].rootN(n);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> pow(final double p) {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].pow(p);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> pow(final int n) {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].pow(n);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> pow(final FieldTuple<T> e) {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].pow(e.values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> exp() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].exp();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> expm1() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].expm1();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> log() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].log();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> log1p() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].log1p();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> log10() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].log10();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> cos() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].cos();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> sin() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].sin();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldSinCos<FieldTuple<T>> sinCos() {
        final FieldTuple<T> sin = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        final FieldTuple<T> cos = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            final FieldSinCos<T> sc = FastMath.sinCos(values[i]);
            sin.values[i] = sc.sin();
            cos.values[i] = sc.cos();
        }
        return new FieldSinCos<>(sin, cos);
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> tan() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].tan();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> acos() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].acos();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> asin() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].asin();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> atan() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].atan();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> atan2(final FieldTuple<T> x) {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].atan2(x.values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> cosh() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].cosh();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> sinh() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].sinh();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> tanh() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].tanh();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> acosh() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].acosh();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> asinh() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].asinh();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> atanh() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].atanh();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> toDegrees() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].toDegrees();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> toRadians() {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = values[i].toRadians();
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> linearCombination(final FieldTuple<T>[] a, final FieldTuple<T>[] b)
        throws MathIllegalArgumentException {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        MathUtils.checkDimension(a.length, b.length);
        final T[] aT = MathArrays.buildArray(values[0].getField(), a.length);
        final T[] bT = MathArrays.buildArray(values[0].getField(), b.length);
        for (int i = 0; i < values.length; ++i) {
            for (int j = 0; j < a.length; ++j) {
                aT[j] = a[j].values[i];
                bT[j] = b[j].values[i];
            }
            result.values[i] = aT[0].linearCombination(aT, bT);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> linearCombination(final double[] a, final FieldTuple<T>[] b)
        throws MathIllegalArgumentException {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        MathUtils.checkDimension(a.length, b.length);
        final T[] bT = MathArrays.buildArray(values[0].getField(), b.length);
        for (int i = 0; i < values.length; ++i) {
            for (int j = 0; j < a.length; ++j) {
                bT[j] = b[j].values[i];
            }
            result.values[i] = bT[0].linearCombination(a, bT);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> linearCombination(final FieldTuple<T> a1, final FieldTuple<T> b1,
                                           final FieldTuple<T> a2, final FieldTuple<T> b2) {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = a1.values[0].linearCombination(a1.values[i], b1.values[i],
                                                              a2.values[i], b2.values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> linearCombination(final double a1, final FieldTuple<T> b1,
                                           final double a2, final FieldTuple<T> b2) {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = b1.values[0].linearCombination(a1, b1.values[i],
                                                              a2, b2.values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> linearCombination(final FieldTuple<T> a1, final FieldTuple<T> b1,
                                           final FieldTuple<T> a2, final FieldTuple<T> b2,
                                           final FieldTuple<T> a3, final FieldTuple<T> b3) {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = a1.values[0].linearCombination(a1.values[i], b1.values[i],
                                                              a2.values[i], b2.values[i],
                                                              a3.values[i], b3.values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> linearCombination(final double a1, final FieldTuple<T> b1,
                                           final double a2, final FieldTuple<T> b2,
                                           final double a3, final FieldTuple<T> b3) {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = b1.values[0].linearCombination(a1, b1.values[i],
                                                              a2, b2.values[i],
                                                              a3, b3.values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> linearCombination(final FieldTuple<T> a1, final FieldTuple<T> b1,
                                           final FieldTuple<T> a2, final FieldTuple<T> b2,
                                           final FieldTuple<T> a3, final FieldTuple<T> b3,
                                           final FieldTuple<T> a4, final FieldTuple<T> b4) {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = a1.values[0].linearCombination(a1.values[i], b1.values[i],
                                                              a2.values[i], b2.values[i],
                                                              a3.values[i], b3.values[i],
                                                              a4.values[i], b4.values[i]);
        }
        return result;
    }

    /** {@inheritDoc} */
    @Override
    public FieldTuple<T> linearCombination(final double a1, final FieldTuple<T> b1,
                                           final double a2, final FieldTuple<T> b2,
                                           final double a3, final FieldTuple<T> b3,
                                           final double a4, final FieldTuple<T> b4) {
        final FieldTuple<T> result = new FieldTuple<>(field, MathArrays.buildArray(values[0].getField(), values.length));
        for (int i = 0; i < values.length; ++i) {
            result.values[i] = b1.values[0].linearCombination(a1, b1.values[i],
                                                              a2, b2.values[i],
                                                              a3, b3.values[i],
                                                              a4, b4.values[i]);
        }
        return result;
    }

    /** Field for {link FieldTuple} instances.
     * @param <T> the type of the field elements
     */
    private static class FieldTupleField<T extends RealFieldElement<T>> implements Field<FieldTuple<T>> {

        /** Constant function evaluating to 0.0. */
        private final FieldTuple<T> zero;

        /** Constant function evaluating to 1.0. */
        private final FieldTuple<T> one;

        /** Simple constructor.
         * @param field field to which the elements belong
         * @param dimension dimension of the tuple
         */
        FieldTupleField(final Field<T> field, final int dimension) {
            final T[] zeroData = MathArrays.buildArray(field, dimension);
            Arrays.fill(zeroData, field.getZero());
            final T[] oneData  = MathArrays.buildArray(field, dimension);
            Arrays.fill(oneData, field.getOne());
            this.zero = new FieldTuple<>(this, zeroData);
            this.one  = new FieldTuple<>(this, oneData);
        }

        /** {@inheritDoc} */
        @Override
        public FieldTuple<T> getZero() {
            return zero;
        }

        /** {@inheritDoc} */
        @Override
        public FieldTuple<T> getOne() {
            return one;
        }

        /** {@inheritDoc} */
        @SuppressWarnings("unchecked")
        @Override
        public Class<? extends FieldElement<FieldTuple<T>>> getRuntimeClass() {
            return (Class<? extends FieldElement<FieldTuple<T>>>) zero.getClass();
        }

        /** {@inheritDoc} */
        @Override
        public boolean equals(final Object other) {
            if (other instanceof FieldTupleField) {
                @SuppressWarnings("unchecked")
                final FieldTupleField<T> that = (FieldTupleField<T>) other;
                return zero.getDimension() == that.zero.getDimension();
            } else {
                return false;
            }
        }

        /** {@inheritDoc} */
        @Override
        public int hashCode() {
            return 0xb4a533e1 ^ zero.getDimension();
        }

    }

}
