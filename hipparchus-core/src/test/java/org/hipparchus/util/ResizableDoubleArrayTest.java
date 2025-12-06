/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements. See the NOTICE file distributed with this
 * work for additional information regarding copyright ownership. The ASF
 * licenses this file to You under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * https://www.apache.org/licenses/LICENSE-2.0 Unless required by applicable law
 * or agreed to in writing, software distributed under the License is
 * distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 * KIND, either express or implied. See the License for the specific language
 * governing permissions and limitations under the License.
 */
package org.hipparchus.util;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.random.RandomDataGenerator;
import org.hipparchus.util.ResizableDoubleArray.ExpansionMode;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertDoesNotThrow;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * This class contains test cases for the ResizableDoubleArray.
 */
class ResizableDoubleArrayTest {

    protected ResizableDoubleArray da = null;

    // Array used to test rolling
    protected ResizableDoubleArray ra = null;

    @AfterEach
    void tearDown()
        throws Exception {
        da = null;
        ra = null;
    }

    @BeforeEach
    void setUp()
        throws Exception {
        da = new ResizableDoubleArray();
        ra = new ResizableDoubleArray();
    }

    @Test
    void testConstructors() {
        float defaultExpansionFactor = 2.0f;
        double defaultContractionCriteria = 2.5;
        ExpansionMode defaultMode = ResizableDoubleArray.ExpansionMode.MULTIPLICATIVE;

        ResizableDoubleArray testDa = new ResizableDoubleArray(2);
        assertEquals(0, testDa.getNumElements());
        assertEquals(2, testDa.getCapacity());
        assertEquals(defaultExpansionFactor, testDa.getExpansionFactor(), 0);
        assertEquals(defaultContractionCriteria, testDa.getContractionCriterion(), 0);
        assertEquals(defaultMode, testDa.getExpansionMode());
        try {
            da = new ResizableDoubleArray(-1);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        testDa = new ResizableDoubleArray((double[]) null);
        assertEquals(0, testDa.getNumElements());

        double[] initialArray = new double[] { 0, 1, 2 };

        testDa = new ResizableDoubleArray(initialArray);
        assertEquals(3, testDa.getNumElements());

        testDa = new ResizableDoubleArray(2, 2.0);
        assertEquals(0, testDa.getNumElements());
        assertEquals(2, testDa.getCapacity());
        assertEquals(defaultExpansionFactor, testDa.getExpansionFactor(), 0);
        assertEquals(defaultContractionCriteria, testDa.getContractionCriterion(), 0);
        assertEquals(defaultMode, testDa.getExpansionMode());

        try {
            da = new ResizableDoubleArray(2, 0.5);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        testDa = new ResizableDoubleArray(2, 3.0);
        assertEquals(3.0f, testDa.getExpansionFactor(), 0);
        assertEquals(3.5f, testDa.getContractionCriterion(), 0);

        testDa = new ResizableDoubleArray(2, 2.0, 3.0);
        assertEquals(0, testDa.getNumElements());
        assertEquals(2, testDa.getCapacity());
        assertEquals(defaultExpansionFactor, testDa.getExpansionFactor(), 0);
        assertEquals(3.0f, testDa.getContractionCriterion(), 0);
        assertEquals(defaultMode, testDa.getExpansionMode());

        try {
            da = new ResizableDoubleArray(2, 2.0, 1.5);
            fail("Expecting MathIllegalArgumentException");
        } catch (MathIllegalArgumentException ex) {
            // expected
        }

        testDa = new ResizableDoubleArray(2, 2.0, 3.0, ResizableDoubleArray.ExpansionMode.ADDITIVE);
        assertEquals(0, testDa.getNumElements());
        assertEquals(2, testDa.getCapacity());
        assertEquals(defaultExpansionFactor, testDa.getExpansionFactor(), 0);
        assertEquals(3.0f, testDa.getContractionCriterion(), 0);
        assertEquals(ResizableDoubleArray.ExpansionMode.ADDITIVE, testDa.getExpansionMode());

        try {
            da = new ResizableDoubleArray(2, 2.0d, 2.5d, null);
            fail("Expecting NullArgumentException");
        } catch (NullArgumentException ex) {
            // expected
        }

        // Copy constructor
        testDa = new ResizableDoubleArray(2, 2.0, 3.0, ResizableDoubleArray.ExpansionMode.ADDITIVE);
        testDa.addElement(2.0);
        testDa.addElement(3.2);
        ResizableDoubleArray copyDa = new ResizableDoubleArray(testDa);
        assertEquals(copyDa, testDa);
        assertEquals(testDa, copyDa);

        // JIRA: MATH-1252
        final double[] values = { 1 };
        testDa = new ResizableDoubleArray(values);
        assertArrayEquals(values, testDa.getElements(), 0);
        assertEquals(1, testDa.getNumElements());
        assertEquals(1, testDa.getElement(0), 0);
    }

    @Test
    void testGetValues() {
        double[] controlArray = {
            2.0, 4.0, 6.0
        };

        da.addElement(2.0);
        da.addElement(4.0);
        da.addElement(6.0);

        double[] testArray = da.getElements();

        for (int i = 0; i < da.getNumElements(); i++) {
            assertEquals(testArray[i],
                         controlArray[i], Double.MIN_VALUE, "The testArray values should equal the controlArray values, " +
                         "index i: " + i + " does not match");
        }
    }

    @Test
    void testMinMax() {
        da.addElement(2.0);
        da.addElement(22.0);
        da.addElement(-2.0);
        da.addElement(21.0);
        da.addElement(22.0);
        da.addElement(42.0);
        da.addElement(62.0);
        da.addElement(22.0);
        da.addElement(122.0);
        da.addElement(1212.0);

        assertEquals(-2.0,
                     Arrays.stream(da.getElements()).min().getAsDouble(),
                     Double.MIN_VALUE,
                     "Min should be -2.0");

        assertEquals(1212.0,
                     Arrays.stream(da.getElements()).max().getAsDouble(),
                     Double.MIN_VALUE,
                     "Max should be 1212.0");
    }

    @Test
    void testSetElementArbitraryExpansion1() {

        // MULTIPLICATIVE_MODE
        da.addElement(2.0);
        da.addElement(4.0);
        da.addElement(6.0);
        da.setElement(1, 3.0);

        // Expand the array arbitrarily to 1000 items
        da.setElement(1000, 3.4);

        assertEquals(1001, da.getNumElements(), "The number of elements should now be 1001, it isn't");

        assertEquals(0.0, da.getElement(760), Double.MIN_VALUE, "Uninitialized Elements are default value of 0.0, index 766 wasn't");

        assertEquals(3.4,
                     da.getElement(1000), Double.MIN_VALUE, "The 1000th index should be 3.4, it isn't");
        assertEquals(2.0,
                     da.getElement(0), Double.MIN_VALUE, "The 0th index should be 2.0, it isn't");
    }

    @Test
    void testSetElementArbitraryExpansion2() {
        // Make sure numElements and expansion work correctly for expansion
        // boundary cases
        da.addElement(2.0);
        da.addElement(4.0);
        da.addElement(6.0);
        assertEquals(16, da.getCapacity());
        assertEquals(3, da.getNumElements());
        da.setElement(3, 7.0);
        assertEquals(16, da.getCapacity());
        assertEquals(4, da.getNumElements());
        da.setElement(10, 10.0);
        assertEquals(16, da.getCapacity());
        assertEquals(11, da.getNumElements());
        da.setElement(9, 10.0);
        assertEquals(16, da.getCapacity());
        assertEquals(11, da.getNumElements());

        try {
            da.setElement(-2, 3);
            fail("Expecting ArrayIndexOutOfBoundsException for negative index");
        } catch (ArrayIndexOutOfBoundsException ex) {
            // expected
        }

        // ADDITIVE_MODE

        ResizableDoubleArray testDa =
                new ResizableDoubleArray(2, 2.0, 3.0, ResizableDoubleArray.ExpansionMode.ADDITIVE);
        assertEquals(2, testDa.getCapacity());
        testDa.addElement(1d);
        testDa.addElement(1d);
        assertEquals(2, testDa.getCapacity());
        testDa.addElement(1d);
        assertEquals(4, testDa.getCapacity());
    }

    @Test
    void testAdd1000() {
        for (int i = 0; i < 1000; i++) {
            da.addElement(i);
        }

        assertEquals(1000, da.getNumElements(), "Number of elements should be equal to 1000 after adding 1000 values");

        assertEquals(56.0,
                     da.getElement(56), Double.MIN_VALUE, "The element at the 56th index should be 56");

        assertEquals(1024,
                     da.getCapacity(),
                     "Internal Storage length should be 1024 if we started out with initial capacity of " +
                     "16 and an expansion factor of 2.0");
    }

    @Test
    void testAddElements() {
        ResizableDoubleArray testDa = new ResizableDoubleArray();

        // MULTIPLICATIVE_MODE
        testDa.addElements(new double[] { 4, 5, 6 });
        assertEquals(3, testDa.getNumElements(), 0);
        assertEquals(4, testDa.getElement(0), 0);
        assertEquals(5, testDa.getElement(1), 0);
        assertEquals(6, testDa.getElement(2), 0);

        testDa.addElements(new double[] { 4, 5, 6 });
        assertEquals(6, testDa.getNumElements());

        // ADDITIVE_MODE (x's are occupied storage locations, 0's are open)
        testDa = new ResizableDoubleArray(2, 2.0, 2.5, ResizableDoubleArray.ExpansionMode.ADDITIVE);
        assertEquals(2, testDa.getCapacity());
        testDa.addElements(new double[] { 1d }); // x,0
        testDa.addElements(new double[] { 2d }); // x,x
        testDa.addElements(new double[] { 3d }); // x,x,x,0 -- expanded
        assertEquals(1d, testDa.getElement(0), 0);
        assertEquals(2d, testDa.getElement(1), 0);
        assertEquals(3d, testDa.getElement(2), 0);
        assertEquals(4, testDa.getCapacity()); // x,x,x,0
        assertEquals(3, testDa.getNumElements());
    }

    @Test
    void testAddElementRolling() {
        ra.addElement(0.5);
        ra.addElement(1.0);
        ra.addElement(1.0);
        ra.addElement(1.0);
        ra.addElement(1.0);
        ra.addElement(1.0);
        ra.addElementRolling(2.0);

        assertEquals(6,
                     ra.getNumElements(),
                     "There should be 6 elements in the eda");
        assertEquals(2.0,
                     Arrays.stream(ra.getElements()).max().getAsDouble(),
                     Double.MIN_VALUE,
                     "The max element should be 2.0");
        assertEquals(1.0,
                     Arrays.stream(ra.getElements()).min().getAsDouble(),
                     Double.MIN_VALUE,
                     "The min element should be 1.0");

        for (int i = 0; i < 1024; i++) {
            ra.addElementRolling(i);
        }

        assertEquals(6, ra.getNumElements(), "We just inserted 1024 rolling elements, num elements should still be 6");

        // MULTIPLICATIVE_MODE
        da.clear();
        da.addElement(1);
        da.addElement(2);
        da.addElementRolling(3);
        assertEquals(3, da.getElement(1), 0);
        da.addElementRolling(4);
        assertEquals(3, da.getElement(0), 0);
        assertEquals(4, da.getElement(1), 0);
        da.addElement(5);
        assertEquals(5, da.getElement(2), 0);
        da.addElementRolling(6);
        assertEquals(4, da.getElement(0), 0);
        assertEquals(5, da.getElement(1), 0);
        assertEquals(6, da.getElement(2), 0);

        // ADDITIVE_MODE (x's are occupied storage locations, 0's are open)
        ResizableDoubleArray testDa =
                new ResizableDoubleArray(2, 2.0, 2.5, ResizableDoubleArray.ExpansionMode.ADDITIVE);
        assertEquals(2, testDa.getCapacity());
        testDa.addElement(1d); // x,0
        testDa.addElement(2d); // x,x
        testDa.addElement(3d); // x,x,x,0 -- expanded
        assertEquals(1d, testDa.getElement(0), 0);
        assertEquals(2d, testDa.getElement(1), 0);
        assertEquals(3d, testDa.getElement(2), 0);
        assertEquals(4, testDa.getCapacity()); // x,x,x,0
        assertEquals(3, testDa.getNumElements());
        testDa.addElementRolling(4d);
        assertEquals(2d, testDa.getElement(0), 0);
        assertEquals(3d, testDa.getElement(1), 0);
        assertEquals(4d, testDa.getElement(2), 0);
        assertEquals(4, testDa.getCapacity()); // 0,x,x,x
        assertEquals(3, testDa.getNumElements());
        testDa.addElementRolling(5d); // 0,0,x,x,x,0 -- time to contract
        assertEquals(3d, testDa.getElement(0), 0);
        assertEquals(4d, testDa.getElement(1), 0);
        assertEquals(5d, testDa.getElement(2), 0);
        assertEquals(4, testDa.getCapacity()); // contracted -- x,x,x,0
        assertEquals(3, testDa.getNumElements());
        try {
            testDa.getElement(4);
            fail("Expecting ArrayIndexOutOfBoundsException");
        } catch (ArrayIndexOutOfBoundsException ex) {
            // expected
        }
        try {
            testDa.getElement(-1);
            fail("Expecting ArrayIndexOutOfBoundsException");
        } catch (ArrayIndexOutOfBoundsException ex) {
            // expected
        }
    }

    @Test
    void testSetNumberOfElements() {
        da.addElement(1.0);
        da.addElement(1.0);
        da.addElement(1.0);
        da.addElement(1.0);
        da.addElement(1.0);
        da.addElement(1.0);
        assertEquals(6, da.getNumElements(), "Number of elements should equal 6");

        da.setNumElements(3);
        assertEquals(3, da.getNumElements(), "Number of elements should equal 3");

        try {
            da.setNumElements(-3);
            fail("Setting number of elements to negative should've thrown an exception");
        } catch (MathIllegalArgumentException iae) {
        }

        da.setNumElements(1024);
        assertEquals(1024, da.getNumElements(), "Number of elements should now be 1024");
        assertEquals(0.0, da.getElement(453), Double.MIN_VALUE, "Element 453 should be a default double");
    }

    @Test
    void testWithInitialCapacity() {

        ResizableDoubleArray eDA2 = new ResizableDoubleArray(2);
        assertEquals(0, eDA2.getNumElements(), "Initial number of elements should be 0");

        final RandomDataGenerator gen = new RandomDataGenerator(1000);
        final int iterations = gen.nextInt(100, 1000);

        for (int i = 0; i < iterations; i++) {
            eDA2.addElement(i);
        }

        assertEquals(iterations, eDA2.getNumElements(), "Number of elements should be equal to " + iterations);

        eDA2.addElement(2.0);

        assertEquals(iterations + 1, eDA2.getNumElements(), "Number of elements should be equals to " +
                     (iterations + 1));
    }

    @Test
    void testWithInitialCapacityAndExpansionFactor() {

        ResizableDoubleArray eDA3 = new ResizableDoubleArray(3, 3.0, 3.5);
        assertEquals(0, eDA3.getNumElements(), "Initial number of elements should be 0");

        final RandomDataGenerator gen = new RandomDataGenerator(1000);
        final int iterations = gen.nextInt(100, 1000);

        for (int i = 0; i < iterations; i++) {
            eDA3.addElement(i);
        }

        assertEquals(iterations, eDA3.getNumElements(), "Number of elements should be equal to " + iterations);

        eDA3.addElement(2.0);

        assertEquals(iterations + 1, eDA3.getNumElements(), "Number of elements should be equals to " +
                     (iterations + 1));

        assertEquals(3.0f,
                     eDA3.getExpansionFactor(), Double.MIN_VALUE, "Expansion factor should equal 3.0");
    }

    @Test
    void testDiscard() {
        da.addElement(2.0);
        da.addElement(2.0);
        da.addElement(2.0);
        da.addElement(2.0);
        da.addElement(2.0);
        da.addElement(2.0);
        da.addElement(2.0);
        da.addElement(2.0);
        da.addElement(2.0);
        da.addElement(2.0);
        da.addElement(2.0);
        assertEquals(11, da.getNumElements(), "Number of elements should be 11");

        da.discardFrontElements(5);
        assertEquals(6, da.getNumElements(), "Number of elements should be 6");

        da.addElement(2.0);
        da.addElement(2.0);
        da.addElement(2.0);
        da.addElement(2.0);
        assertEquals(10, da.getNumElements(), "Number of elements should be 10");

        da.discardMostRecentElements(2);
        assertEquals(8, da.getNumElements(), "Number of elements should be 8");

        try {
            da.discardFrontElements(-1);
            fail("Trying to discard a negative number of element is not allowed");
        } catch (Exception e) {
        }

        try {
            da.discardMostRecentElements(-1);
            fail("Trying to discard a negative number of element is not allowed");
        } catch (Exception e) {
        }

        try {
            da.discardFrontElements(10000);
            fail("You can't discard more elements than the array contains");
        } catch (Exception e) {
        }

        try {
            da.discardMostRecentElements(10000);
            fail("You can't discard more elements than the array contains");
        } catch (Exception e) {
        }

    }

    @Test
    void testSubstitute() {

        da.addElement(2.0);
        da.addElement(2.0);
        da.addElement(2.0);
        da.addElement(2.0);
        da.addElement(2.0);
        da.addElement(2.0);
        da.addElement(2.0);
        da.addElement(2.0);
        da.addElement(2.0);
        da.addElement(2.0);
        da.addElement(2.0);
        assertEquals(11, da.getNumElements(), "Number of elements should be 11");

        da.substituteMostRecentElement(24);

        assertEquals(11, da.getNumElements(), "Number of elements should be 11");

        assertDoesNotThrow(() -> da.discardMostRecentElements(10), "Trying to discard a negative number of element is not allowed");

        da.substituteMostRecentElement(24);

        assertEquals(1, da.getNumElements(), "Number of elements should be 1");

    }

    @SuppressWarnings("unlikely-arg-type")
    @Test
    void testEqualsAndHashCode()
        throws Exception {

        // Wrong type
        ResizableDoubleArray first = new ResizableDoubleArray();
        Double other = Double.valueOf(2);
        assertNotEquals(first, other);

        // Null
        other = null;
        assertNotEquals(first, other);

        // Reflexive
        assertEquals(first, first);

        // Argumentless constructor
        ResizableDoubleArray second = new ResizableDoubleArray();
        verifyEquality(first, second);

        // Equals iff same data, same properties
        ResizableDoubleArray third = new ResizableDoubleArray(3, 2.0, 2.0);
        verifyInequality(third, first);
        ResizableDoubleArray fourth = new ResizableDoubleArray(3, 2.0, 2.0);
        ResizableDoubleArray fifth = new ResizableDoubleArray(2, 2.0, 2.0);
        verifyEquality(third, fourth);
        verifyInequality(third, fifth);
        third.addElement(4.1);
        third.addElement(4.2);
        third.addElement(4.3);
        fourth.addElement(4.1);
        fourth.addElement(4.2);
        fourth.addElement(4.3);
        verifyEquality(third, fourth);

        // expand
        fourth.addElement(4.4);
        verifyInequality(third, fourth);
        third.addElement(4.4);
        verifyEquality(third, fourth);
        fourth.addElement(4.4);
        verifyInequality(third, fourth);
        third.addElement(4.4);
        verifyEquality(third, fourth);
        fourth.addElementRolling(4.5);
        third.addElementRolling(4.5);
        verifyEquality(third, fourth);

        // discard
        third.discardFrontElements(1);
        verifyInequality(third, fourth);
        fourth.discardFrontElements(1);
        verifyEquality(third, fourth);

        // discard recent
        third.discardMostRecentElements(2);
        fourth.discardMostRecentElements(2);
        verifyEquality(third, fourth);

        // wrong order
        third.addElement(18);
        fourth.addElement(17);
        third.addElement(17);
        fourth.addElement(18);
        verifyInequality(third, fourth);

        // Copy constructor
        verifyEquality(fourth, new ResizableDoubleArray(fourth));

        // Instance copy
        verifyEquality(fourth, fourth.copy());
    }

    @Test
    void testGetArrayRef() {
        final ResizableDoubleArray a = new ResizableDoubleArray();

        // Modify "a" through the public API.
        final int index = 20;
        final double v1 = 1.2;
        a.setElement(index, v1);

        // Modify the internal storage through the protected API.
        final double v2 = v1 + 3.4;
        final double[] aInternalArray = a.getArrayRef();
        aInternalArray[a.getStartIndex() + index] = v2;

        assertEquals(v2, a.getElement(index), 0d);
    }

    @Test
    void testCompute() {
        final ResizableDoubleArray a = new ResizableDoubleArray();
        final int max = 20;
        for (int i = 1; i <= max; i++) {
            a.setElement(i, i);
        }

        final MathArrays.Function add = new MathArrays.Function() {

            @Override
            public double evaluate(double[] a, int index, int num) {
                double sum = 0;
                final int max = index + num;
                for (int i = index; i < max; i++) {
                    sum += a[i];
                }
                return sum;
            }

            @Override
            public double evaluate(double[] a) {
                return evaluate(a, 0, a.length);
            }
        };

        final double sum = a.compute(add);
        assertEquals(0.5 * max * (max + 1), sum, 0);
    }

    private void verifyEquality(ResizableDoubleArray a,
                                ResizableDoubleArray b) {
        assertEquals(b, a);
        assertEquals(a, b);
        assertEquals(a.hashCode(), b.hashCode());
    }

    private void verifyInequality(ResizableDoubleArray a,
                                  ResizableDoubleArray b) {
        assertNotEquals(b, a);
        assertNotEquals(a, b);
        assertNotEquals(a.hashCode(), b.hashCode());
    }

}
