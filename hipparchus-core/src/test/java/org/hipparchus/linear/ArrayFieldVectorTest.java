/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */
package org.hipparchus.linear;

import org.hipparchus.Field;
import org.hipparchus.FieldElement;
import org.hipparchus.UnitTestUtils;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.fraction.Fraction;
import org.hipparchus.fraction.FractionField;
import org.junit.jupiter.api.Test;

import java.io.Serializable;
import java.lang.reflect.Array;
import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotSame;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test cases for the {@link ArrayFieldVector} class.
 *
 */
class ArrayFieldVectorTest {

    //
    protected Fraction[][] ma1 = {
            {new Fraction(1), new Fraction(2), new Fraction(3)},
            {new Fraction(4), new Fraction(5), new Fraction(6)},
            {new Fraction(7), new Fraction(8), new Fraction(9)}
    };
    protected Fraction[] vec1 = {new Fraction(1), new Fraction(2), new Fraction(3)};
    protected Fraction[] vec2 = {new Fraction(4), new Fraction(5), new Fraction(6)};
    protected Fraction[] vec3 = {new Fraction(7), new Fraction(8), new Fraction(9)};
    protected Fraction[] vec4 = { new Fraction(1), new Fraction(2), new Fraction(3),
                                  new Fraction(4), new Fraction(5), new Fraction(6),
                                  new Fraction(7), new Fraction(8), new Fraction(9)};
    protected Fraction[] vec_null = {Fraction.ZERO, Fraction.ZERO, Fraction.ZERO};
    protected Fraction[] dvec1 = {new Fraction(1), new Fraction(2), new Fraction(3),
                                  new Fraction(4), new Fraction(5), new Fraction(6),
                                  new Fraction(7), new Fraction(8), new Fraction(9)};
    protected Fraction[][] mat1 = {
            {new Fraction(1), new Fraction(2), new Fraction(3)},
            {new Fraction(4), new Fraction(5), new Fraction(6)},
            {new Fraction(7), new Fraction(8), new Fraction(9)}
    };

    // Testclass to test the FieldVector<Fraction> interface
    // only with enough content to support the test
    public static class FieldVectorTestImpl<T extends FieldElement<T>>
        implements FieldVector<T>, Serializable {

        private static final long serialVersionUID = 3970959016014158539L;

        private final Field<T> field;

        /** Entries of the vector. */
        protected T[] data;

        /** Build an array of elements.
         * @param length size of the array to build
         * @return a new array
         */
        @SuppressWarnings("unchecked") // field is of type T
        private T[] buildArray(final int length) {
            return (T[]) Array.newInstance(field.getRuntimeClass(), length);
        }

        public FieldVectorTestImpl(T[] d) {
            field = d[0].getField();
            data = d.clone();
        }

        public Field<T> getField() {
            return field;
        }

        private UnsupportedOperationException unsupported() {
            return new UnsupportedOperationException("Not supported, unneeded for test purposes");
        }

        public FieldVector<T> copy() {
            throw unsupported();
        }

        public FieldVector<T> add(FieldVector<T> v) {
            throw unsupported();
        }

        public FieldVector<T> add(T[] v) {
            throw unsupported();
        }

        public FieldVector<T> subtract(FieldVector<T> v) {
            throw unsupported();
        }

        public FieldVector<T> subtract(T[] v) {
            throw unsupported();
        }

        public FieldVector<T> mapAdd(T d) {
            throw unsupported();
        }

        public FieldVector<T> mapAddToSelf(T d) {
            throw unsupported();
        }

        public FieldVector<T> mapSubtract(T d) {
            throw unsupported();
        }

        public FieldVector<T> mapSubtractToSelf(T d) {
            throw unsupported();
        }

        public FieldVector<T> mapMultiply(T d) {
            T[] out = buildArray(data.length);
            for (int i = 0; i < data.length; i++) {
                out[i] = data[i].multiply(d);
            }
            return new FieldVectorTestImpl<T>(out);
        }

        public FieldVector<T> mapMultiplyToSelf(T d) {
            throw unsupported();
        }

        public FieldVector<T> mapDivide(T d) {
            throw unsupported();
        }

        public FieldVector<T> mapDivideToSelf(T d) {
            throw unsupported();
        }

        public FieldVector<T> mapInv() {
            throw unsupported();
        }

        public FieldVector<T> mapInvToSelf() {
            throw unsupported();
        }

        public FieldVector<T> ebeMultiply(FieldVector<T> v) {
            throw unsupported();
        }

        public FieldVector<T> ebeMultiply(T[] v) {
            throw unsupported();
        }

        public FieldVector<T> ebeDivide(FieldVector<T> v) {
            throw unsupported();
        }

        public FieldVector<T> ebeDivide(T[] v) {
            throw unsupported();
        }

        public T[] getData() {
            return data.clone();
        }

        public T dotProduct(FieldVector<T> v) {
            T dot = field.getZero();
            for (int i = 0; i < data.length; i++) {
                dot = dot.add(data[i].multiply(v.getEntry(i)));
            }
            return dot;
        }

        public T dotProduct(T[] v) {
            T dot = field.getZero();
            for (int i = 0; i < data.length; i++) {
                dot = dot.add(data[i].multiply(v[i]));
            }
            return dot;
        }

        public FieldVector<T> projection(FieldVector<T> v) {
            throw unsupported();
        }

        public FieldVector<T> projection(T[] v) {
            throw unsupported();
        }

        public FieldMatrix<T> outerProduct(FieldVector<T> v) {
            throw unsupported();
        }

        public FieldMatrix<T> outerProduct(T[] v) {
            throw unsupported();
        }

        public T getEntry(int index) {
            return data[index];
        }

        public int getDimension() {
            return data.length;
        }

        public FieldVector<T> append(FieldVector<T> v) {
            throw unsupported();
        }

        public FieldVector<T> append(T d) {
            throw unsupported();
        }

        public FieldVector<T> append(T[] a) {
            throw unsupported();
        }

        public FieldVector<T> getSubVector(int index, int n) {
            throw unsupported();
        }

        public void setEntry(int index, T value) {
            throw unsupported();
        }

        public void setSubVector(int index, FieldVector<T> v) {
            throw unsupported();
        }

        public void setSubVector(int index, T[] v) {
            throw unsupported();
        }

        public void set(T value) {
            throw unsupported();
        }

        public T[] toArray() {
            return data.clone();
        }

    }

    @Test
    void testConstructors() {

        ArrayFieldVector<Fraction> v0 = new ArrayFieldVector<Fraction>(FractionField.getInstance());
        assertEquals(0, v0.getDimension());

        ArrayFieldVector<Fraction> v1 = new ArrayFieldVector<Fraction>(FractionField.getInstance(), 7);
        assertEquals(7, v1.getDimension());
        assertEquals(Fraction.ZERO, v1.getEntry(6));

        ArrayFieldVector<Fraction> v2 = new ArrayFieldVector<Fraction>(5, new Fraction(123, 100));
        assertEquals(5, v2.getDimension());
        assertEquals(new Fraction(123, 100), v2.getEntry(4));

        ArrayFieldVector<Fraction> v3 = new ArrayFieldVector<Fraction>(FractionField.getInstance(), vec1);
        assertEquals(3, v3.getDimension());
        assertEquals(new Fraction(2), v3.getEntry(1));

        ArrayFieldVector<Fraction> v4 = new ArrayFieldVector<Fraction>(FractionField.getInstance(), vec4, 3, 2);
        assertEquals(2, v4.getDimension());
        assertEquals(new Fraction(4), v4.getEntry(0));
        try {
            new ArrayFieldVector<Fraction>(vec4, 8, 3);
            fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected behavior
        }

        FieldVector<Fraction> v5_i = new ArrayFieldVector<Fraction>(dvec1);
        assertEquals(9, v5_i.getDimension());
        assertEquals(new Fraction(9), v5_i.getEntry(8));

        ArrayFieldVector<Fraction> v5 = new ArrayFieldVector<Fraction>(dvec1);
        assertEquals(9, v5.getDimension());
        assertEquals(new Fraction(9), v5.getEntry(8));

        ArrayFieldVector<Fraction> v6 = new ArrayFieldVector<Fraction>(dvec1, 3, 2);
        assertEquals(2, v6.getDimension());
        assertEquals(new Fraction(4), v6.getEntry(0));
        try {
            new ArrayFieldVector<Fraction>(dvec1, 8, 3);
            fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected behavior
        }

        ArrayFieldVector<Fraction> v7 = new ArrayFieldVector<Fraction>(v1);
        assertEquals(7, v7.getDimension());
        assertEquals(Fraction.ZERO, v7.getEntry(6));

        FieldVectorTestImpl<Fraction> v7_i = new FieldVectorTestImpl<Fraction>(vec1);

        ArrayFieldVector<Fraction> v7_2 = new ArrayFieldVector<Fraction>(v7_i);
        assertEquals(3, v7_2.getDimension());
        assertEquals(new Fraction(2), v7_2.getEntry(1));

        ArrayFieldVector<Fraction> v8 = new ArrayFieldVector<Fraction>(v1, true);
        assertEquals(7, v8.getDimension());
        assertEquals(Fraction.ZERO, v8.getEntry(6));
        assertNotSame(v1.getDataRef(), v8.getDataRef(), "testData not same object ");

        ArrayFieldVector<Fraction> v8_2 = new ArrayFieldVector<Fraction>(v1, false);
        assertEquals(7, v8_2.getDimension());
        assertEquals(Fraction.ZERO, v8_2.getEntry(6));
        assertArrayEquals(v1.getDataRef(), v8_2.getDataRef());

        ArrayFieldVector<Fraction> v9 = new ArrayFieldVector<Fraction>((FieldVector<Fraction>) v1, (FieldVector<Fraction>) v3);
        assertEquals(10, v9.getDimension());
        assertEquals(new Fraction(1), v9.getEntry(7));

    }

    @Test
    void testDataInOut() {

        ArrayFieldVector<Fraction> v1 = new ArrayFieldVector<Fraction>(vec1);
        ArrayFieldVector<Fraction> v2 = new ArrayFieldVector<Fraction>(vec2);
        ArrayFieldVector<Fraction> v4 = new ArrayFieldVector<Fraction>(vec4);
        FieldVectorTestImpl<Fraction> v2_t = new FieldVectorTestImpl<Fraction>(vec2);

        FieldVector<Fraction> v_append_1 = v1.append(v2);
        assertEquals(6, v_append_1.getDimension());
        assertEquals(new Fraction(4), v_append_1.getEntry(3));

        FieldVector<Fraction> v_append_2 = v1.append(new Fraction(2));
        assertEquals(4, v_append_2.getDimension());
        assertEquals(new Fraction(2), v_append_2.getEntry(3));

        FieldVector<Fraction> v_append_4 = v1.append(v2_t);
        assertEquals(6, v_append_4.getDimension());
        assertEquals(new Fraction(4), v_append_4.getEntry(3));

        FieldVector<Fraction> v_copy = v1.copy();
        assertEquals(3, v_copy.getDimension());
        assertNotSame(v1.getDataRef(), v_copy.toArray(), "testData not same object ");

        Fraction[] a_frac = v1.toArray();
        assertEquals(3, a_frac.length);
        assertNotSame(v1.getDataRef(), a_frac, "testData not same object ");


//      ArrayFieldVector<Fraction> vout4 = (ArrayFieldVector<Fraction>) v1.clone();
//      Assertions.assertEquals(3, vout4.getDimension());
//      Assertions.assertEquals(v1.getDataRef(), vout4.getDataRef());


        FieldVector<Fraction> vout5 = v4.getSubVector(3, 3);
        assertEquals(3, vout5.getDimension());
        assertEquals(new Fraction(5), vout5.getEntry(1));
        try {
            v4.getSubVector(3, 7);
            fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected behavior
        }

        ArrayFieldVector<Fraction> v_set1 = (ArrayFieldVector<Fraction>) v1.copy();
        v_set1.setEntry(1, new Fraction(11));
        assertEquals(new Fraction(11), v_set1.getEntry(1));
        try {
            v_set1.setEntry(3, new Fraction(11));
            fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected behavior
        }

        ArrayFieldVector<Fraction> v_set2 = (ArrayFieldVector<Fraction>) v4.copy();
        v_set2.set(3, v1);
        assertEquals(new Fraction(1), v_set2.getEntry(3));
        assertEquals(new Fraction(7), v_set2.getEntry(6));
        try {
            v_set2.set(7, v1);
            fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected behavior
        }

        ArrayFieldVector<Fraction> v_set3 = (ArrayFieldVector<Fraction>) v1.copy();
        v_set3.set(new Fraction(13));
        assertEquals(new Fraction(13), v_set3.getEntry(2));

        try {
            v_set3.getEntry(23);
            fail("ArrayIndexOutOfBoundsException expected");
        } catch (ArrayIndexOutOfBoundsException ex) {
            // expected behavior
        }

        ArrayFieldVector<Fraction> v_set4 = (ArrayFieldVector<Fraction>) v4.copy();
        v_set4.setSubVector(3, v2_t);
        assertEquals(new Fraction(4), v_set4.getEntry(3));
        assertEquals(new Fraction(7), v_set4.getEntry(6));
        try {
            v_set4.setSubVector(7, v2_t);
            fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected behavior
        }


        ArrayFieldVector<Fraction> vout10 = (ArrayFieldVector<Fraction>) v1.copy();
        ArrayFieldVector<Fraction> vout10_2 = (ArrayFieldVector<Fraction>) v1.copy();
        assertEquals(vout10, vout10_2);
        vout10_2.setEntry(0, new Fraction(11, 10));
        assertNotSame(vout10, vout10_2);

    }

    @Test
    void testMapFunctions() {
        ArrayFieldVector<Fraction> v1 = new ArrayFieldVector<Fraction>(vec1);

        //octave =  v1 .+ 2.0
        FieldVector<Fraction> v_mapAdd = v1.mapAdd(new Fraction(2));
        Fraction[] result_mapAdd = {new Fraction(3), new Fraction(4), new Fraction(5)};
        checkArray("compare vectors" ,result_mapAdd,v_mapAdd.toArray());

        //octave =  v1 .+ 2.0
        FieldVector<Fraction> v_mapAddToSelf = v1.copy();
        v_mapAddToSelf.mapAddToSelf(new Fraction(2));
        Fraction[] result_mapAddToSelf = {new Fraction(3), new Fraction(4), new Fraction(5)};
        checkArray("compare vectors" ,result_mapAddToSelf,v_mapAddToSelf.toArray());

        //octave =  v1 .- 2.0
        FieldVector<Fraction> v_mapSubtract = v1.mapSubtract(new Fraction(2));
        Fraction[] result_mapSubtract = {new Fraction(-1), Fraction.ZERO, new Fraction(1)};
        checkArray("compare vectors" ,result_mapSubtract,v_mapSubtract.toArray());

        //octave =  v1 .- 2.0
        FieldVector<Fraction> v_mapSubtractToSelf = v1.copy();
        v_mapSubtractToSelf.mapSubtractToSelf(new Fraction(2));
        Fraction[] result_mapSubtractToSelf = {new Fraction(-1), Fraction.ZERO, new Fraction(1)};
        checkArray("compare vectors" ,result_mapSubtractToSelf,v_mapSubtractToSelf.toArray());

        //octave =  v1 .* 2.0
        FieldVector<Fraction> v_mapMultiply = v1.mapMultiply(new Fraction(2));
        Fraction[] result_mapMultiply = {new Fraction(2), new Fraction(4), new Fraction(6)};
        checkArray("compare vectors" ,result_mapMultiply,v_mapMultiply.toArray());

        //octave =  v1 .* 2.0
        FieldVector<Fraction> v_mapMultiplyToSelf = v1.copy();
        v_mapMultiplyToSelf.mapMultiplyToSelf(new Fraction(2));
        Fraction[] result_mapMultiplyToSelf = {new Fraction(2), new Fraction(4), new Fraction(6)};
        checkArray("compare vectors" ,result_mapMultiplyToSelf,v_mapMultiplyToSelf.toArray());

        //octave =  v1 ./ 2.0
        FieldVector<Fraction> v_mapDivide = v1.mapDivide(new Fraction(2));
        Fraction[] result_mapDivide = {new Fraction(1, 2), new Fraction(1), new Fraction(3, 2)};
        checkArray("compare vectors" ,result_mapDivide,v_mapDivide.toArray());

        //octave =  v1 ./ 2.0
        FieldVector<Fraction> v_mapDivideToSelf = v1.copy();
        v_mapDivideToSelf.mapDivideToSelf(new Fraction(2));
        Fraction[] result_mapDivideToSelf = {new Fraction(1, 2), new Fraction(1), new Fraction(3, 2)};
        checkArray("compare vectors" ,result_mapDivideToSelf,v_mapDivideToSelf.toArray());

        //octave =  v1 .^-1
        FieldVector<Fraction> v_mapInv = v1.mapInv();
        Fraction[] result_mapInv = {new Fraction(1),new Fraction(1, 2),new Fraction(1, 3)};
        checkArray("compare vectors" ,result_mapInv,v_mapInv.toArray());

        //octave =  v1 .^-1
        FieldVector<Fraction> v_mapInvToSelf = v1.copy();
        v_mapInvToSelf.mapInvToSelf();
        Fraction[] result_mapInvToSelf = {new Fraction(1),new Fraction(1, 2),new Fraction(1, 3)};
        checkArray("compare vectors" ,result_mapInvToSelf,v_mapInvToSelf.toArray());

    }

    @Test
    void testBasicFunctions() {
        ArrayFieldVector<Fraction> v1 = new ArrayFieldVector<Fraction>(vec1);
        ArrayFieldVector<Fraction> v2 = new ArrayFieldVector<Fraction>(vec2);
        new ArrayFieldVector<Fraction>(vec_null);

        FieldVectorTestImpl<Fraction> v2_t = new FieldVectorTestImpl<Fraction>(vec2);

        //octave =  v1 + v2
        ArrayFieldVector<Fraction> v_add = v1.add(v2);
        Fraction[] result_add = {new Fraction(5), new Fraction(7), new Fraction(9)};
        checkArray("compare vect" ,v_add.toArray(),result_add);

        FieldVectorTestImpl<Fraction> vt2 = new FieldVectorTestImpl<Fraction>(vec2);
        FieldVector<Fraction> v_add_i = v1.add(vt2);
        Fraction[] result_add_i = {new Fraction(5), new Fraction(7), new Fraction(9)};
        checkArray("compare vect" ,v_add_i.toArray(),result_add_i);

        //octave =  v1 - v2
        ArrayFieldVector<Fraction> v_subtract = v1.subtract(v2);
        Fraction[] result_subtract = {new Fraction(-3), new Fraction(-3), new Fraction(-3)};
        checkArray("compare vect" ,v_subtract.toArray(),result_subtract);

        FieldVector<Fraction> v_subtract_i = v1.subtract(vt2);
        Fraction[] result_subtract_i = {new Fraction(-3), new Fraction(-3), new Fraction(-3)};
        checkArray("compare vect" ,v_subtract_i.toArray(),result_subtract_i);

        // octave v1 .* v2
        ArrayFieldVector<Fraction>  v_ebeMultiply = v1.ebeMultiply(v2);
        Fraction[] result_ebeMultiply = {new Fraction(4), new Fraction(10), new Fraction(18)};
        checkArray("compare vect" ,v_ebeMultiply.toArray(),result_ebeMultiply);

        FieldVector<Fraction>  v_ebeMultiply_2 = v1.ebeMultiply(v2_t);
        Fraction[] result_ebeMultiply_2 = {new Fraction(4), new Fraction(10), new Fraction(18)};
        checkArray("compare vect" ,v_ebeMultiply_2.toArray(),result_ebeMultiply_2);

        // octave v1 ./ v2
        ArrayFieldVector<Fraction>  v_ebeDivide = v1.ebeDivide(v2);
        Fraction[] result_ebeDivide = {new Fraction(1, 4), new Fraction(2, 5), new Fraction(1, 2)};
        checkArray("compare vect" ,v_ebeDivide.toArray(),result_ebeDivide);

        FieldVector<Fraction>  v_ebeDivide_2 = v1.ebeDivide(v2_t);
        Fraction[] result_ebeDivide_2 = {new Fraction(1, 4), new Fraction(2, 5), new Fraction(1, 2)};
        checkArray("compare vect" ,v_ebeDivide_2.toArray(),result_ebeDivide_2);

        // octave  dot(v1,v2)
        Fraction dot =  v1.dotProduct(v2);
        assertEquals(new Fraction(32), dot, "compare val ");

        // octave  dot(v1,v2_t)
        Fraction dot_2 =  v1.dotProduct(v2_t);
        assertEquals(new Fraction(32), dot_2, "compare val ");

        FieldMatrix<Fraction> m_outerProduct = v1.outerProduct(v2);
        assertEquals(new Fraction(4), m_outerProduct.getEntry(0,0), "compare val ");

        FieldMatrix<Fraction> m_outerProduct_2 = v1.outerProduct(v2_t);
        assertEquals(new Fraction(4), m_outerProduct_2.getEntry(0,0), "compare val ");

        ArrayFieldVector<Fraction> v_projection = v1.projection(v2);
        Fraction[] result_projection = {new Fraction(128, 77), new Fraction(160, 77), new Fraction(192, 77)};
        checkArray("compare vect", v_projection.toArray(), result_projection);

        FieldVector<Fraction> v_projection_2 = v1.projection(v2_t);
        Fraction[] result_projection_2 = {new Fraction(128, 77), new Fraction(160, 77), new Fraction(192, 77)};
        checkArray("compare vect", v_projection_2.toArray(), result_projection_2);

     }

    @Test
    void testMisc() {
        ArrayFieldVector<Fraction> v1 = new ArrayFieldVector<Fraction>(vec1);
        ArrayFieldVector<Fraction> v4 = new ArrayFieldVector<Fraction>(vec4);
        FieldVector<Fraction> v4_2 = new ArrayFieldVector<Fraction>(vec4);

        assertEquals("{1; 2; 3}",  v1.toString());
        assertEquals("{1; 2; 3; 4; 5; 6; 7; 8; 9}",  v4.toString());

        /*
         Fraction[] dout1 = v1.copyOut();
        Assertions.assertEquals(3, dout1.length);
        assertNotSame("testData not same object ", v1.getDataRef(), dout1);
         */
        try {
            v1.checkVectorDimensions(2);
            fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected behavior
        }

       try {
            v1.checkVectorDimensions(v4);
            fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected behavior
        }

        try {
            v1.checkVectorDimensions(v4_2);
            fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected behavior
        }

    }

    @Test
    void testSerial()  {
        ArrayFieldVector<Fraction> v = new ArrayFieldVector<Fraction>(vec1);
        assertEquals(v,UnitTestUtils.serializeAndRecover(v));
    }

    @Test
    void testZeroVectors() {

        // when the field is not specified, array cannot be empty
        try {
            new ArrayFieldVector<Fraction>(new Fraction[0]);
            fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected behavior
        }
        try {
            new ArrayFieldVector<Fraction>(new Fraction[0], true);
            fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected behavior
        }
        try {
            new ArrayFieldVector<Fraction>(new Fraction[0], false);
            fail("MathIllegalArgumentException expected");
        } catch (MathIllegalArgumentException ex) {
            // expected behavior
        }

        // when the field is specified, array can be empty
        assertEquals(0, new ArrayFieldVector<Fraction>(FractionField.getInstance(), new Fraction[0]).getDimension());
        assertEquals(0, new ArrayFieldVector<Fraction>(FractionField.getInstance(), new Fraction[0], true).getDimension());
        assertEquals(0, new ArrayFieldVector<Fraction>(FractionField.getInstance(), new Fraction[0], false).getDimension());

    }

    @Test
    void testOuterProduct() {
        final ArrayFieldVector<Fraction> u
            = new ArrayFieldVector<Fraction>(FractionField.getInstance(),
                                             new Fraction[] {new Fraction(1),
                                                             new Fraction(2),
                                                             new Fraction(-3)});
        final ArrayFieldVector<Fraction> v
            = new ArrayFieldVector<Fraction>(FractionField.getInstance(),
                                             new Fraction[] {new Fraction(4),
                                                             new Fraction(-2)});

        final FieldMatrix<Fraction> uv = u.outerProduct(v);

        final double tol = Math.ulp(1d);
        assertEquals(new Fraction(4).doubleValue(), uv.getEntry(0, 0).doubleValue(), tol);
        assertEquals(new Fraction(-2).doubleValue(), uv.getEntry(0, 1).doubleValue(), tol);
        assertEquals(new Fraction(8).doubleValue(), uv.getEntry(1, 0).doubleValue(), tol);
        assertEquals(new Fraction(-4).doubleValue(), uv.getEntry(1, 1).doubleValue(), tol);
        assertEquals(new Fraction(-12).doubleValue(), uv.getEntry(2, 0).doubleValue(), tol);
        assertEquals(new Fraction(6).doubleValue(), uv.getEntry(2, 1).doubleValue(), tol);
    }

    /** verifies that two vectors are equals */
    protected void checkArray(String msg, Fraction[] m, Fraction[] n) {
        if (m.length != n.length) {
            fail("vectors have different lengths");
        }
        for (int i = 0; i < m.length; i++) {
            assertEquals(m[i],n[i],msg + " " +  i + " elements differ");
        }
    }

    /*
     * TESTS OF THE VISITOR PATTERN
     */

    /** The whole vector is visited. */
    @Test
    void testWalkInDefaultOrderPreservingVisitor1() {
        final Fraction[] data = new Fraction[] {
            Fraction.ZERO, Fraction.ONE, Fraction.ZERO,
            Fraction.ZERO, Fraction.TWO, Fraction.ZERO,
            Fraction.ZERO, Fraction.ZERO, new Fraction(3)
        };
        final ArrayFieldVector<Fraction> v = new ArrayFieldVector<Fraction>(data);
        final FieldVectorPreservingVisitor<Fraction> visitor;
        visitor = new FieldVectorPreservingVisitor<Fraction>() {

            private int expectedIndex;

            public void visit(final int actualIndex, final Fraction actualValue) {
                assertEquals(expectedIndex, actualIndex);
                assertEquals(data[actualIndex], actualValue, Integer.toString(actualIndex));
                ++expectedIndex;
            }

            public void start(final int actualSize, final int actualStart,
                              final int actualEnd) {
                assertEquals(data.length, actualSize);
                assertEquals(0, actualStart);
                assertEquals(data.length - 1, actualEnd);
                expectedIndex = 0;
            }

            public Fraction end() {
                return Fraction.ZERO;
            }
        };
        v.walkInDefaultOrder(visitor);
    }

    /** Visiting an invalid subvector. */
    @Test
    void testWalkInDefaultOrderPreservingVisitor2() {
        final ArrayFieldVector<Fraction> v = create(5);
        final FieldVectorPreservingVisitor<Fraction> visitor;
        visitor = new FieldVectorPreservingVisitor<Fraction>() {

            public void visit(int index, Fraction value) {
                // Do nothing
            }

            public void start(int dimension, int start, int end) {
                // Do nothing
            }

            public Fraction end() {
                return Fraction.ZERO;
            }
        };
        try {
            v.walkInDefaultOrder(visitor, -1, 4);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 5, 4);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 0, -1);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 0, 5);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 4, 0);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
    }

    /** Visiting a valid subvector. */
    @Test
    void testWalkInDefaultOrderPreservingVisitor3() {
        final Fraction[] data = new Fraction[] {
            Fraction.ZERO, Fraction.ONE, Fraction.ZERO,
            Fraction.ZERO, Fraction.TWO, Fraction.ZERO,
            Fraction.ZERO, Fraction.ZERO, new Fraction(3)
        };
        final ArrayFieldVector<Fraction> v = new ArrayFieldVector<Fraction>(data);
        final int expectedStart = 2;
        final int expectedEnd = 7;
        final FieldVectorPreservingVisitor<Fraction> visitor;
        visitor = new FieldVectorPreservingVisitor<Fraction>() {

            private int expectedIndex;

            public void visit(final int actualIndex, final Fraction actualValue) {
                assertEquals(expectedIndex, actualIndex);
                assertEquals(data[actualIndex], actualValue, Integer.toString(actualIndex));
                ++expectedIndex;
            }

            public void start(final int actualSize, final int actualStart,
                              final int actualEnd) {
                assertEquals(data.length, actualSize);
                assertEquals(expectedStart, actualStart);
                assertEquals(expectedEnd, actualEnd);
                expectedIndex = expectedStart;
            }

            public Fraction end() {
                return Fraction.ZERO;
            }
        };
        v.walkInDefaultOrder(visitor, expectedStart, expectedEnd);
    }

    /** The whole vector is visited. */
    @Test
    void testWalkInOptimizedOrderPreservingVisitor1() {
        final Fraction[] data = new Fraction[] {
            Fraction.ZERO, Fraction.ONE, Fraction.ZERO,
            Fraction.ZERO, Fraction.TWO, Fraction.ZERO,
            Fraction.ZERO, Fraction.ZERO, new Fraction(3)
        };
        final ArrayFieldVector<Fraction> v = new ArrayFieldVector<Fraction>(data);
        final FieldVectorPreservingVisitor<Fraction> visitor;
        visitor = new FieldVectorPreservingVisitor<Fraction>() {
            private final boolean[] visited = new boolean[data.length];

            public void visit(final int actualIndex, final Fraction actualValue) {
                visited[actualIndex] = true;
                assertEquals(data[actualIndex], actualValue, Integer.toString(actualIndex));
            }

            public void start(final int actualSize, final int actualStart,
                              final int actualEnd) {
                assertEquals(data.length, actualSize);
                assertEquals(0, actualStart);
                assertEquals(data.length - 1, actualEnd);
                Arrays.fill(visited, false);
            }

            public Fraction end() {
                for (int i = 0; i < data.length; i++) {
                    assertTrue(visited[i],
                                      "entry " + i + "has not been visited");
                }
                return Fraction.ZERO;
            }
        };
        v.walkInOptimizedOrder(visitor);
    }

    /** Visiting an invalid subvector. */
    @Test
    void testWalkInOptimizedOrderPreservingVisitor2() {
        final ArrayFieldVector<Fraction> v = create(5);
        final FieldVectorPreservingVisitor<Fraction> visitor;
        visitor = new FieldVectorPreservingVisitor<Fraction>() {

            public void visit(int index, Fraction value) {
                // Do nothing
            }

            public void start(int dimension, int start, int end) {
                // Do nothing
            }

            public Fraction end() {
                return Fraction.ZERO;
            }
        };
        try {
            v.walkInOptimizedOrder(visitor, -1, 4);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 5, 4);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 0, -1);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 0, 5);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 4, 0);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
    }

    /** Visiting a valid subvector. */
    @Test
    void testWalkInOptimizedOrderPreservingVisitor3() {
        final Fraction[] data = new Fraction[] {
            Fraction.ZERO, Fraction.ONE, Fraction.ZERO,
            Fraction.ZERO, Fraction.TWO, Fraction.ZERO,
            Fraction.ZERO, Fraction.ZERO, new Fraction(3)
        };
        final ArrayFieldVector<Fraction> v = new ArrayFieldVector<Fraction>(data);
        final int expectedStart = 2;
        final int expectedEnd = 7;
        final FieldVectorPreservingVisitor<Fraction> visitor;
        visitor = new FieldVectorPreservingVisitor<Fraction>() {
            private final boolean[] visited = new boolean[data.length];

            public void visit(final int actualIndex, final Fraction actualValue) {
                assertEquals(data[actualIndex], actualValue, Integer.toString(actualIndex));
                visited[actualIndex] = true;
            }

            public void start(final int actualSize, final int actualStart,
                              final int actualEnd) {
                assertEquals(data.length, actualSize);
                assertEquals(expectedStart, actualStart);
                assertEquals(expectedEnd, actualEnd);
                Arrays.fill(visited, true);
            }

            public Fraction end() {
                for (int i = expectedStart; i <= expectedEnd; i++) {
                    assertTrue(visited[i],
                                      "entry " + i + "has not been visited");
                }
                return Fraction.ZERO;
            }
        };
        v.walkInOptimizedOrder(visitor, expectedStart, expectedEnd);
    }

    /** The whole vector is visited. */
    @Test
    void testWalkInDefaultOrderChangingVisitor1() {
        final Fraction[] data = new Fraction[] {
            Fraction.ZERO, Fraction.ONE, Fraction.ZERO,
            Fraction.ZERO, Fraction.TWO, Fraction.ZERO,
            Fraction.ZERO, Fraction.ZERO, new Fraction(3)
        };
        final ArrayFieldVector<Fraction> v = new ArrayFieldVector<Fraction>(data);
        final FieldVectorChangingVisitor<Fraction> visitor;
        visitor = new FieldVectorChangingVisitor<Fraction>() {

            private int expectedIndex;

            public Fraction visit(final int actualIndex, final Fraction actualValue) {
                assertEquals(expectedIndex, actualIndex);
                assertEquals(data[actualIndex], actualValue, Integer.toString(actualIndex));
                ++expectedIndex;
                return actualValue.add(actualIndex);
            }

            public void start(final int actualSize, final int actualStart,
                              final int actualEnd) {
                assertEquals(data.length, actualSize);
                assertEquals(0, actualStart);
                assertEquals(data.length - 1, actualEnd);
                expectedIndex = 0;
            }

            public Fraction end() {
                return Fraction.ZERO;
            }
        };
        v.walkInDefaultOrder(visitor);
        for (int i = 0; i < data.length; i++) {
            assertEquals(data[i].add(i), v.getEntry(i), "entry " + i);
        }
    }

    /** Visiting an invalid subvector. */
    @Test
    void testWalkInDefaultOrderChangingVisitor2() {
        final ArrayFieldVector<Fraction> v = create(5);
        final FieldVectorChangingVisitor<Fraction> visitor;
        visitor = new FieldVectorChangingVisitor<Fraction>() {

            public Fraction visit(int index, Fraction value) {
                return Fraction.ZERO;
            }

            public void start(int dimension, int start, int end) {
                // Do nothing
            }

            public Fraction end() {
                return Fraction.ZERO;
            }
        };
        try {
            v.walkInDefaultOrder(visitor, -1, 4);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 5, 4);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 0, -1);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 0, 5);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInDefaultOrder(visitor, 4, 0);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
    }

    /** Visiting a valid subvector. */
    @Test
    void testWalkInDefaultOrderChangingVisitor3() {
        final Fraction[] data = new Fraction[] {
            Fraction.ZERO, Fraction.ONE, Fraction.ZERO,
            Fraction.ZERO, Fraction.TWO, Fraction.ZERO,
            Fraction.ZERO, Fraction.ZERO, new Fraction(3)
        };
        final ArrayFieldVector<Fraction> v = new ArrayFieldVector<Fraction>(data);
        final int expectedStart = 2;
        final int expectedEnd = 7;
        final FieldVectorChangingVisitor<Fraction> visitor;
        visitor = new FieldVectorChangingVisitor<Fraction>() {

            private int expectedIndex;

            public Fraction visit(final int actualIndex, final Fraction actualValue) {
                assertEquals(expectedIndex, actualIndex);
                assertEquals(data[actualIndex], actualValue, Integer.toString(actualIndex));
                ++expectedIndex;
                return actualValue.add(actualIndex);
            }

            public void start(final int actualSize, final int actualStart,
                              final int actualEnd) {
                assertEquals(data.length, actualSize);
                assertEquals(expectedStart, actualStart);
                assertEquals(expectedEnd, actualEnd);
                expectedIndex = expectedStart;
            }

            public Fraction end() {
                return Fraction.ZERO;
            }
        };
        v.walkInDefaultOrder(visitor, expectedStart, expectedEnd);
        for (int i = expectedStart; i <= expectedEnd; i++) {
            assertEquals(data[i].add(i), v.getEntry(i), "entry " + i);
        }
    }

    /** The whole vector is visited. */
    @Test
    void testWalkInOptimizedOrderChangingVisitor1() {
        final Fraction[] data = new Fraction[] {
            Fraction.ZERO, Fraction.ONE, Fraction.ZERO,
            Fraction.ZERO, Fraction.TWO, Fraction.ZERO,
            Fraction.ZERO, Fraction.ZERO, new Fraction(3)
        };
        final ArrayFieldVector<Fraction> v = new ArrayFieldVector<Fraction>(data);
        final FieldVectorChangingVisitor<Fraction> visitor;
        visitor = new FieldVectorChangingVisitor<Fraction>() {
            private final boolean[] visited = new boolean[data.length];

            public Fraction visit(final int actualIndex, final Fraction actualValue) {
                visited[actualIndex] = true;
                assertEquals(data[actualIndex], actualValue, Integer.toString(actualIndex));
                return actualValue.add(actualIndex);
            }

            public void start(final int actualSize, final int actualStart,
                              final int actualEnd) {
                assertEquals(data.length, actualSize);
                assertEquals(0, actualStart);
                assertEquals(data.length - 1, actualEnd);
                Arrays.fill(visited, false);
            }

            public Fraction end() {
                for (int i = 0; i < data.length; i++) {
                    assertTrue(visited[i],
                                      "entry " + i + "has not been visited");
                }
                return Fraction.ZERO;
            }
        };
        v.walkInOptimizedOrder(visitor);
        for (int i = 0; i < data.length; i++) {
            assertEquals(data[i].add(i), v.getEntry(i), "entry " + i);
        }
    }

    /** Visiting an invalid subvector. */
    @Test
    void testWalkInOptimizedOrderChangingVisitor2() {
        final ArrayFieldVector<Fraction> v = create(5);
        final FieldVectorChangingVisitor<Fraction> visitor;
        visitor = new FieldVectorChangingVisitor<Fraction>() {

            public Fraction visit(int index, Fraction value) {
                return Fraction.ZERO;
            }

            public void start(int dimension, int start, int end) {
                // Do nothing
            }

            public Fraction end() {
                return Fraction.ZERO;
            }
        };
        try {
            v.walkInOptimizedOrder(visitor, -1, 4);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 5, 4);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 0, -1);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 0, 5);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
        try {
            v.walkInOptimizedOrder(visitor, 4, 0);
            fail();
        } catch (MathIllegalArgumentException e) {
            // Expected behavior
        }
    }

    /** Visiting a valid subvector. */
    @Test
    void testWalkInOptimizedOrderChangingVisitor3() {
        final Fraction[] data = new Fraction[] {
            Fraction.ZERO, Fraction.ONE, Fraction.ZERO,
            Fraction.ZERO, Fraction.TWO, Fraction.ZERO,
            Fraction.ZERO, Fraction.ZERO, new Fraction(3)
        };
        final ArrayFieldVector<Fraction> v = new ArrayFieldVector<Fraction>(data);
        final int expectedStart = 2;
        final int expectedEnd = 7;
        final FieldVectorChangingVisitor<Fraction> visitor;
        visitor = new FieldVectorChangingVisitor<Fraction>() {
            private final boolean[] visited = new boolean[data.length];

            public Fraction visit(final int actualIndex, final Fraction actualValue) {
                assertEquals(data[actualIndex], actualValue, Integer.toString(actualIndex));
                visited[actualIndex] = true;
                return actualValue.add(actualIndex);
            }

            public void start(final int actualSize, final int actualStart,
                              final int actualEnd) {
                assertEquals(data.length, actualSize);
                assertEquals(expectedStart, actualStart);
                assertEquals(expectedEnd, actualEnd);
                Arrays.fill(visited, true);
            }

            public Fraction end() {
                for (int i = expectedStart; i <= expectedEnd; i++) {
                    assertTrue(visited[i],
                                      "entry " + i + "has not been visited");
                }
                return Fraction.ZERO;
            }
        };
        v.walkInOptimizedOrder(visitor, expectedStart, expectedEnd);
        for (int i = expectedStart; i <= expectedEnd; i++) {
            assertEquals(data[i].add(i), v.getEntry(i), "entry " + i);
        }
    }

    private ArrayFieldVector<Fraction> create(int n) {
        Fraction[] t = new Fraction[n];
        for (int i = 0; i < n; ++i) {
            t[i] = Fraction.ZERO;
        }
        return new ArrayFieldVector<Fraction>(t);
    }
}
