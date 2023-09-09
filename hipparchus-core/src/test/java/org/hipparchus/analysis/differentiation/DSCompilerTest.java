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

package org.hipparchus.analysis.differentiation;

import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.Stream;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.CombinatoricsUtils;
import org.junit.Assert;
import org.junit.Test;


/**
 * Test for class {@link DSCompiler}.
 */
public class DSCompilerTest {

    @Test
    public void testSize() {
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                long expected = CombinatoricsUtils.binomialCoefficient(i + j, i);
                Assert.assertEquals(expected, DSCompiler.getCompiler(i, j).getSize());
                Assert.assertEquals(expected, DSCompiler.getCompiler(j, i).getSize());
            }
        }
    }

    @Test
    public void testIndices() {

        DSCompiler c = DSCompiler.getCompiler(0, 0);
        checkIndices(c.getPartialDerivativeOrders(0), new int[0]);

        c = DSCompiler.getCompiler(0, 1);
        checkIndices(c.getPartialDerivativeOrders(0), new int[0]);

        c = DSCompiler.getCompiler(1, 0);
        checkIndices(c.getPartialDerivativeOrders(0), 0);

        c = DSCompiler.getCompiler(1, 1);
        checkIndices(c.getPartialDerivativeOrders(0), 0);
        checkIndices(c.getPartialDerivativeOrders(1), 1);

        c = DSCompiler.getCompiler(1, 2);
        checkIndices(c.getPartialDerivativeOrders(0), 0);
        checkIndices(c.getPartialDerivativeOrders(1), 1);
        checkIndices(c.getPartialDerivativeOrders(2), 2);

        c = DSCompiler.getCompiler(2, 1);
        checkIndices(c.getPartialDerivativeOrders(0), 0, 0);
        checkIndices(c.getPartialDerivativeOrders(1), 1, 0);
        checkIndices(c.getPartialDerivativeOrders(2), 0, 1);

        c = DSCompiler.getCompiler(1, 3);
        checkIndices(c.getPartialDerivativeOrders(0), 0);
        checkIndices(c.getPartialDerivativeOrders(1), 1);
        checkIndices(c.getPartialDerivativeOrders(2), 2);
        checkIndices(c.getPartialDerivativeOrders(3), 3);

        c = DSCompiler.getCompiler(2, 2);
        checkIndices(c.getPartialDerivativeOrders(0), 0, 0);
        checkIndices(c.getPartialDerivativeOrders(1), 1, 0);
        checkIndices(c.getPartialDerivativeOrders(2), 2, 0);
        checkIndices(c.getPartialDerivativeOrders(3), 0, 1);
        checkIndices(c.getPartialDerivativeOrders(4), 1, 1);
        checkIndices(c.getPartialDerivativeOrders(5), 0, 2);

        c = DSCompiler.getCompiler(3, 1);
        checkIndices(c.getPartialDerivativeOrders(0), 0, 0, 0);
        checkIndices(c.getPartialDerivativeOrders(1), 1, 0, 0);
        checkIndices(c.getPartialDerivativeOrders(2), 0, 1, 0);
        checkIndices(c.getPartialDerivativeOrders(3), 0, 0, 1);

        c = DSCompiler.getCompiler(1, 4);
        checkIndices(c.getPartialDerivativeOrders(0), 0);
        checkIndices(c.getPartialDerivativeOrders(1), 1);
        checkIndices(c.getPartialDerivativeOrders(2), 2);
        checkIndices(c.getPartialDerivativeOrders(3), 3);
        checkIndices(c.getPartialDerivativeOrders(4), 4);

        c = DSCompiler.getCompiler(2, 3);
        checkIndices(c.getPartialDerivativeOrders(0), 0, 0);
        checkIndices(c.getPartialDerivativeOrders(1), 1, 0);
        checkIndices(c.getPartialDerivativeOrders(2), 2, 0);
        checkIndices(c.getPartialDerivativeOrders(3), 3, 0);
        checkIndices(c.getPartialDerivativeOrders(4), 0, 1);
        checkIndices(c.getPartialDerivativeOrders(5), 1, 1);
        checkIndices(c.getPartialDerivativeOrders(6), 2, 1);
        checkIndices(c.getPartialDerivativeOrders(7), 0, 2);
        checkIndices(c.getPartialDerivativeOrders(8), 1, 2);
        checkIndices(c.getPartialDerivativeOrders(9), 0, 3);

        c = DSCompiler.getCompiler(3, 2);
        checkIndices(c.getPartialDerivativeOrders(0), 0, 0, 0);
        checkIndices(c.getPartialDerivativeOrders(1), 1, 0, 0);
        checkIndices(c.getPartialDerivativeOrders(2), 2, 0, 0);
        checkIndices(c.getPartialDerivativeOrders(3), 0, 1, 0);
        checkIndices(c.getPartialDerivativeOrders(4), 1, 1, 0);
        checkIndices(c.getPartialDerivativeOrders(5), 0, 2, 0);
        checkIndices(c.getPartialDerivativeOrders(6), 0, 0, 1);
        checkIndices(c.getPartialDerivativeOrders(7), 1, 0, 1);
        checkIndices(c.getPartialDerivativeOrders(8), 0, 1, 1);
        checkIndices(c.getPartialDerivativeOrders(9), 0, 0, 2);

        c = DSCompiler.getCompiler(4, 1);
        checkIndices(c.getPartialDerivativeOrders(0), 0, 0, 0, 0);
        checkIndices(c.getPartialDerivativeOrders(1), 1, 0, 0, 0);
        checkIndices(c.getPartialDerivativeOrders(2), 0, 1, 0, 0);
        checkIndices(c.getPartialDerivativeOrders(3), 0, 0, 1, 0);
        checkIndices(c.getPartialDerivativeOrders(4), 0, 0, 0, 1);

    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testIncompatibleParams() {
        DSCompiler.getCompiler(3, 2).checkCompatibility(DSCompiler.getCompiler(4, 2));
    }

    @Test(expected=MathIllegalArgumentException.class)
    public void testIncompatibleOrder() {
        DSCompiler.getCompiler(3, 3).checkCompatibility(DSCompiler.getCompiler(3, 2));
    }

    @Test
    public void testSymmetry() {
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                DSCompiler c = DSCompiler.getCompiler(i, j);
                for (int k = 0; k < c.getSize(); ++k) {
                    Assert.assertEquals(k, c.getPartialDerivativeIndex(c.getPartialDerivativeOrders(k)));
                }
            }
        }
    }

    @Test
    public void testMultiplicationRules()
        throws SecurityException, NoSuchFieldException, IllegalArgumentException,
               IllegalAccessException, NoSuchMethodException, InvocationTargetException {

        Map<String,String> referenceRules = new HashMap<String, String>();
        referenceRules.put("(f*g)",             "f * g");
        referenceRules.put("∂(f*g)/∂p₀",        "f * ∂g/∂p₀ + ∂f/∂p₀ * g");
        referenceRules.put("∂(f*g)/∂p₁",        referenceRules.get("∂(f*g)/∂p₀").replaceAll("p₀", "p₁"));
        referenceRules.put("∂(f*g)/∂p₂",        referenceRules.get("∂(f*g)/∂p₀").replaceAll("p₀", "p₂"));
        referenceRules.put("∂(f*g)/∂p₃",        referenceRules.get("∂(f*g)/∂p₀").replaceAll("p₀", "p₃"));
        referenceRules.put("∂²(f*g)/∂p₀²",      "f * ∂²g/∂p₀² + 2 * ∂f/∂p₀ * ∂g/∂p₀ + ∂²f/∂p₀² * g");
        referenceRules.put("∂²(f*g)/∂p₁²",      referenceRules.get("∂²(f*g)/∂p₀²").replaceAll("p₀", "p₁"));
        referenceRules.put("∂²(f*g)/∂p₂²",      referenceRules.get("∂²(f*g)/∂p₀²").replaceAll("p₀", "p₂"));
        referenceRules.put("∂²(f*g)/∂p₃²",      referenceRules.get("∂²(f*g)/∂p₀²").replaceAll("p₀", "p₃"));
        referenceRules.put("∂²(f*g)/∂p₀∂p₁",    "f * ∂²g/∂p₀∂p₁ + ∂f/∂p₁ * ∂g/∂p₀ + ∂f/∂p₀ * ∂g/∂p₁ + ∂²f/∂p₀∂p₁ * g");
        referenceRules.put("∂²(f*g)/∂p₀∂p₂",    referenceRules.get("∂²(f*g)/∂p₀∂p₁").replaceAll("p₁", "p₂"));
        referenceRules.put("∂²(f*g)/∂p₀∂p₃",    referenceRules.get("∂²(f*g)/∂p₀∂p₁").replaceAll("p₁", "p₃"));
        referenceRules.put("∂²(f*g)/∂p₁∂p₂",    referenceRules.get("∂²(f*g)/∂p₀∂p₂").replaceAll("p₀", "p₁"));
        referenceRules.put("∂²(f*g)/∂p₁∂p₃",    referenceRules.get("∂²(f*g)/∂p₀∂p₃").replaceAll("p₀", "p₁"));
        referenceRules.put("∂²(f*g)/∂p₂∂p₃",    referenceRules.get("∂²(f*g)/∂p₀∂p₃").replaceAll("p₀", "p₂"));
        referenceRules.put("∂³(f*g)/∂p₀³",      "f * ∂³g/∂p₀³ +" +
                                                " 3 * ∂f/∂p₀ * ∂²g/∂p₀² +" +
                                                " 3 * ∂²f/∂p₀² * ∂g/∂p₀ +" +
                                                " ∂³f/∂p₀³ * g");
        referenceRules.put("∂³(f*g)/∂p₁³",      referenceRules.get("∂³(f*g)/∂p₀³").replaceAll("p₀", "p₁"));
        referenceRules.put("∂³(f*g)/∂p₂³",      referenceRules.get("∂³(f*g)/∂p₀³").replaceAll("p₀", "p₂"));
        referenceRules.put("∂³(f*g)/∂p₃³",      referenceRules.get("∂³(f*g)/∂p₀³").replaceAll("p₀", "p₃"));
        referenceRules.put("∂³(f*g)/∂p₀²∂p₁",   "f * ∂³g/∂p₀²∂p₁ +" +
                                                " ∂f/∂p₁ * ∂²g/∂p₀² +" +
                                                " 2 * ∂f/∂p₀ * ∂²g/∂p₀∂p₁ +" +
                                                " 2 * ∂²f/∂p₀∂p₁ * ∂g/∂p₀ +" +
                                                " ∂²f/∂p₀² * ∂g/∂p₁ +" +
                                                " ∂³f/∂p₀²∂p₁ * g");
        referenceRules.put("∂³(f*g)/∂p₀∂p₁²",   "f * ∂³g/∂p₀∂p₁² +" +
                                                " 2 * ∂f/∂p₁ * ∂²g/∂p₀∂p₁ +" +
                                                " ∂²f/∂p₁² * ∂g/∂p₀ +" +
                                                " ∂f/∂p₀ * ∂²g/∂p₁² +" +
                                                " 2 * ∂²f/∂p₀∂p₁ * ∂g/∂p₁ +" +
                                                " ∂³f/∂p₀∂p₁² * g");
        referenceRules.put("∂³(f*g)/∂p₀²∂p₂",   referenceRules.get("∂³(f*g)/∂p₀²∂p₁").replaceAll("p₁", "p₂"));
        referenceRules.put("∂³(f*g)/∂p₁²∂p₂",   referenceRules.get("∂³(f*g)/∂p₀²∂p₂").replaceAll("p₀", "p₁"));
        referenceRules.put("∂³(f*g)/∂p₀∂p₂²",   referenceRules.get("∂³(f*g)/∂p₀∂p₁²").replaceAll("p₁", "p₂"));
        referenceRules.put("∂³(f*g)/∂p₁∂p₂²",   referenceRules.get("∂³(f*g)/∂p₀∂p₂²").replaceAll("p₀", "p₁"));
        referenceRules.put("∂³(f*g)/∂p₀²∂p₃",   referenceRules.get("∂³(f*g)/∂p₀²∂p₂").replaceAll("p₂", "p₃"));
        referenceRules.put("∂³(f*g)/∂p₁²∂p₃",   referenceRules.get("∂³(f*g)/∂p₀²∂p₃").replaceAll("p₀", "p₁"));
        referenceRules.put("∂³(f*g)/∂p₂²∂p₃",   referenceRules.get("∂³(f*g)/∂p₀²∂p₃").replaceAll("p₀", "p₂"));
        referenceRules.put("∂³(f*g)/∂p₀∂p₃²",   referenceRules.get("∂³(f*g)/∂p₀∂p₁²").replaceAll("p₁", "p₃"));
        referenceRules.put("∂³(f*g)/∂p₁∂p₃²",   referenceRules.get("∂³(f*g)/∂p₀∂p₃²").replaceAll("p₀", "p₁"));
        referenceRules.put("∂³(f*g)/∂p₂∂p₃²",   referenceRules.get("∂³(f*g)/∂p₀∂p₃²").replaceAll("p₀", "p₂"));
        referenceRules.put("∂³(f*g)/∂p₀∂p₁∂p₂", "f * ∂³g/∂p₀∂p₁∂p₂ +" +
                                                " ∂f/∂p₂ * ∂²g/∂p₀∂p₁ +" +
                                                " ∂f/∂p₁ * ∂²g/∂p₀∂p₂ +" +
                                                " ∂²f/∂p₁∂p₂ * ∂g/∂p₀ +" +
                                                " ∂f/∂p₀ * ∂²g/∂p₁∂p₂ +" +
                                                " ∂²f/∂p₀∂p₂ * ∂g/∂p₁ +" +
                                                " ∂²f/∂p₀∂p₁ * ∂g/∂p₂ +" +
                                                " ∂³f/∂p₀∂p₁∂p₂ * g");
        referenceRules.put("∂³(f*g)/∂p₀∂p₁∂p₃", referenceRules.get("∂³(f*g)/∂p₀∂p₁∂p₂").replaceAll("p₂", "p₃"));
        referenceRules.put("∂³(f*g)/∂p₀∂p₂∂p₃", referenceRules.get("∂³(f*g)/∂p₀∂p₁∂p₃").replaceAll("p₁", "p₂"));
        referenceRules.put("∂³(f*g)/∂p₁∂p₂∂p₃", referenceRules.get("∂³(f*g)/∂p₀∂p₂∂p₃").replaceAll("p₀", "p₁"));

        Field multFieldArrayField = DSCompiler.class.getDeclaredField("multIndirection");
        multFieldArrayField.setAccessible(true);
        Class<?> abstractMapperClass = Stream.
                        of(DSCompiler.class.getDeclaredClasses()).
                        filter(c -> c.getName().endsWith("AbstractMapper")).
                        findAny().
                        get();
        Class<?> multiplicationMapperClass = Stream.
                        of(DSCompiler.class.getDeclaredClasses()).
                        filter(c -> c.getName().endsWith("MultiplicationMapper")).
                        findAny().
                        get();
        Method coeffMethod = abstractMapperClass.getDeclaredMethod("getCoeff");
        Field lhsField = multiplicationMapperClass.getDeclaredField("lhsIndex");
        lhsField.setAccessible(true);
        Field rhsField = multiplicationMapperClass.getDeclaredField("rhsIndex");
        rhsField.setAccessible(true);
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 4; ++j) {
                DSCompiler compiler = DSCompiler.getCompiler(i, j);
                Object[][] multIndirection = (Object[][]) multFieldArrayField.get(compiler);
                for (int k = 0; k < multIndirection.length; ++k) {
                    String product = ordersToString(compiler.getPartialDerivativeOrders(k),
                                                    "(f*g)", variables("p"));
                    StringBuilder rule = new StringBuilder();
                    for (Object term : multIndirection[k]) {
                        if (rule.length() > 0) {
                            rule.append(" + ");
                        }
                        if (((Integer) coeffMethod.invoke(term)).intValue() > 1) {
                            rule.append(((Integer) coeffMethod.invoke(term)).intValue()).append(" * ");
                        }
                        rule.append(ordersToString(compiler.getPartialDerivativeOrders(((Integer) lhsField.get(term)).intValue()),
                                                   "f", variables("p")));
                        rule.append(" * ");
                        rule.append(ordersToString(compiler.getPartialDerivativeOrders(((Integer) rhsField.get(term)).intValue()),
                                                   "g", variables("p")));
                    }
                    Assert.assertEquals(product, referenceRules.get(product), rule.toString());
                }
            }
        }
    }

    @Test
    public void testCompositionRules()
        throws SecurityException, NoSuchFieldException, IllegalArgumentException,
               IllegalAccessException, InvocationTargetException, NoSuchMethodException {

        // the following reference rules have all been computed independently from the library,
        // using only pencil and paper and some search and replace to handle symmetries
        Map<String,String> referenceRules = new HashMap<String, String>();
        referenceRules.put("(f(g))",              "(f(g))");
        referenceRules.put("∂(f(g))/∂p₀",          "∂(f(g))/∂g * ∂g/∂p₀");
        referenceRules.put("∂(f(g))/∂p₁",          referenceRules.get("∂(f(g))/∂p₀").replaceAll("p₀", "p₁"));
        referenceRules.put("∂(f(g))/∂p₂",          referenceRules.get("∂(f(g))/∂p₀").replaceAll("p₀", "p₂"));
        referenceRules.put("∂(f(g))/∂p₃",          referenceRules.get("∂(f(g))/∂p₀").replaceAll("p₀", "p₃"));
        referenceRules.put("∂²(f(g))/∂p₀²",        "∂²(f(g))/∂g² * ∂g/∂p₀ * ∂g/∂p₀ + ∂(f(g))/∂g * ∂²g/∂p₀²");
        referenceRules.put("∂²(f(g))/∂p₁²",        referenceRules.get("∂²(f(g))/∂p₀²").replaceAll("p₀", "p₁"));
        referenceRules.put("∂²(f(g))/∂p₂²",        referenceRules.get("∂²(f(g))/∂p₀²").replaceAll("p₀", "p₂"));
        referenceRules.put("∂²(f(g))/∂p₃²",        referenceRules.get("∂²(f(g))/∂p₀²").replaceAll("p₀", "p₃"));
        referenceRules.put("∂²(f(g))/∂p₀∂p₁",       "∂²(f(g))/∂g² * ∂g/∂p₀ * ∂g/∂p₁ + ∂(f(g))/∂g * ∂²g/∂p₀∂p₁");
        referenceRules.put("∂²(f(g))/∂p₀∂p₂",       referenceRules.get("∂²(f(g))/∂p₀∂p₁").replaceAll("p₁", "p₂"));
        referenceRules.put("∂²(f(g))/∂p₀∂p₃",       referenceRules.get("∂²(f(g))/∂p₀∂p₁").replaceAll("p₁", "p₃"));
        referenceRules.put("∂²(f(g))/∂p₁∂p₂",       referenceRules.get("∂²(f(g))/∂p₀∂p₂").replaceAll("p₀", "p₁"));
        referenceRules.put("∂²(f(g))/∂p₁∂p₃",       referenceRules.get("∂²(f(g))/∂p₀∂p₃").replaceAll("p₀", "p₁"));
        referenceRules.put("∂²(f(g))/∂p₂∂p₃",       referenceRules.get("∂²(f(g))/∂p₀∂p₃").replaceAll("p₀", "p₂"));
        referenceRules.put("∂³(f(g))/∂p₀³",        "∂³(f(g))/∂g³ * ∂g/∂p₀ * ∂g/∂p₀ * ∂g/∂p₀ +" +
                                                  " 3 * ∂²(f(g))/∂g² * ∂g/∂p₀ * ∂²g/∂p₀² +" +
                                                  " ∂(f(g))/∂g * ∂³g/∂p₀³");
        referenceRules.put("∂³(f(g))/∂p₁³",        referenceRules.get("∂³(f(g))/∂p₀³").replaceAll("p₀", "p₁"));
        referenceRules.put("∂³(f(g))/∂p₂³",        referenceRules.get("∂³(f(g))/∂p₀³").replaceAll("p₀", "p₂"));
        referenceRules.put("∂³(f(g))/∂p₃³",        referenceRules.get("∂³(f(g))/∂p₀³").replaceAll("p₀", "p₃"));
        referenceRules.put("∂³(f(g))/∂p₀∂p₁²",      "∂³(f(g))/∂g³ * ∂g/∂p₀ * ∂g/∂p₁ * ∂g/∂p₁ +" +
                                                  " 2 * ∂²(f(g))/∂g² * ∂g/∂p₁ * ∂²g/∂p₀∂p₁ +" +
                                                  " ∂²(f(g))/∂g² * ∂g/∂p₀ * ∂²g/∂p₁² +" +
                                                  " ∂(f(g))/∂g * ∂³g/∂p₀∂p₁²");
        referenceRules.put("∂³(f(g))/∂p₀∂p₂²",      referenceRules.get("∂³(f(g))/∂p₀∂p₁²").replaceAll("p₁", "p₂"));
        referenceRules.put("∂³(f(g))/∂p₀∂p₃²",      referenceRules.get("∂³(f(g))/∂p₀∂p₁²").replaceAll("p₁", "p₃"));
        referenceRules.put("∂³(f(g))/∂p₁∂p₂²",      referenceRules.get("∂³(f(g))/∂p₀∂p₂²").replaceAll("p₀", "p₁"));
        referenceRules.put("∂³(f(g))/∂p₁∂p₃²",      referenceRules.get("∂³(f(g))/∂p₀∂p₃²").replaceAll("p₀", "p₁"));
        referenceRules.put("∂³(f(g))/∂p₂∂p₃²",      referenceRules.get("∂³(f(g))/∂p₀∂p₃²").replaceAll("p₀", "p₂"));
        referenceRules.put("∂³(f(g))/∂p₀²∂p₁",      "∂³(f(g))/∂g³ * ∂g/∂p₀ * ∂g/∂p₀ * ∂g/∂p₁ +" +
                                                  " 2 * ∂²(f(g))/∂g² * ∂g/∂p₀ * ∂²g/∂p₀∂p₁ +" +
                                                  " ∂²(f(g))/∂g² * ∂²g/∂p₀² * ∂g/∂p₁ +" +
                                                  " ∂(f(g))/∂g * ∂³g/∂p₀²∂p₁");
        referenceRules.put("∂³(f(g))/∂p₀²∂p₂",      referenceRules.get("∂³(f(g))/∂p₀²∂p₁").replaceAll("p₁", "p₂"));
        referenceRules.put("∂³(f(g))/∂p₀²∂p₃",      referenceRules.get("∂³(f(g))/∂p₀²∂p₁").replaceAll("p₁", "p₃"));
        referenceRules.put("∂³(f(g))/∂p₁²∂p₂",      referenceRules.get("∂³(f(g))/∂p₀²∂p₂").replaceAll("p₀", "p₁"));
        referenceRules.put("∂³(f(g))/∂p₁²∂p₃",      referenceRules.get("∂³(f(g))/∂p₀²∂p₃").replaceAll("p₀", "p₁"));
        referenceRules.put("∂³(f(g))/∂p₂²∂p₃",      referenceRules.get("∂³(f(g))/∂p₀²∂p₃").replaceAll("p₀", "p₂"));
        referenceRules.put("∂³(f(g))/∂p₀∂p₁∂p₂",     "∂³(f(g))/∂g³ * ∂g/∂p₀ * ∂g/∂p₁ * ∂g/∂p₂ +" +
                                                  " ∂²(f(g))/∂g² * ∂g/∂p₁ * ∂²g/∂p₀∂p₂ +" +
                                                  " ∂²(f(g))/∂g² * ∂g/∂p₀ * ∂²g/∂p₁∂p₂ +" +
                                                  " ∂²(f(g))/∂g² * ∂²g/∂p₀∂p₁ * ∂g/∂p₂ +" +
                                                  " ∂(f(g))/∂g * ∂³g/∂p₀∂p₁∂p₂");
        referenceRules.put("∂³(f(g))/∂p₀∂p₁∂p₃",     referenceRules.get("∂³(f(g))/∂p₀∂p₁∂p₂").replaceAll("p₂", "p₃"));
        referenceRules.put("∂³(f(g))/∂p₀∂p₂∂p₃",     referenceRules.get("∂³(f(g))/∂p₀∂p₁∂p₃").replaceAll("p₁", "p₂"));
        referenceRules.put("∂³(f(g))/∂p₁∂p₂∂p₃",     referenceRules.get("∂³(f(g))/∂p₀∂p₂∂p₃").replaceAll("p₀", "p₁"));
        referenceRules.put("∂⁴(f(g))/∂p₀⁴",        "∂⁴(f(g))/∂g⁴ * ∂g/∂p₀ * ∂g/∂p₀ * ∂g/∂p₀ * ∂g/∂p₀ +" +
                                                  " 6 * ∂³(f(g))/∂g³ * ∂g/∂p₀ * ∂g/∂p₀ * ∂²g/∂p₀² +" +
                                                  " 3 * ∂²(f(g))/∂g² * ∂²g/∂p₀² * ∂²g/∂p₀² +" +
                                                  " 4 * ∂²(f(g))/∂g² * ∂g/∂p₀ * ∂³g/∂p₀³ +" +
                                                  " ∂(f(g))/∂g * ∂⁴g/∂p₀⁴");
        referenceRules.put("∂⁴(f(g))/∂p₁⁴",        referenceRules.get("∂⁴(f(g))/∂p₀⁴").replaceAll("p₀", "p₁"));
        referenceRules.put("∂⁴(f(g))/∂p₂⁴",        referenceRules.get("∂⁴(f(g))/∂p₀⁴").replaceAll("p₀", "p₂"));
        referenceRules.put("∂⁴(f(g))/∂p₃⁴",        referenceRules.get("∂⁴(f(g))/∂p₀⁴").replaceAll("p₀", "p₃"));
        referenceRules.put("∂⁴(f(g))/∂p₀³∂p₁",      "∂⁴(f(g))/∂g⁴ * ∂g/∂p₀ * ∂g/∂p₀ * ∂g/∂p₀ * ∂g/∂p₁ +" +
                                                  " 3 * ∂³(f(g))/∂g³ * ∂g/∂p₀ * ∂g/∂p₀ * ∂²g/∂p₀∂p₁ +" +
                                                  " 3 * ∂³(f(g))/∂g³ * ∂g/∂p₀ * ∂²g/∂p₀² * ∂g/∂p₁ +" +
                                                  " 3 * ∂²(f(g))/∂g² * ∂²g/∂p₀² * ∂²g/∂p₀∂p₁ +" +
                                                  " 3 * ∂²(f(g))/∂g² * ∂g/∂p₀ * ∂³g/∂p₀²∂p₁ +" +
                                                  " ∂²(f(g))/∂g² * ∂³g/∂p₀³ * ∂g/∂p₁ +" +
                                                  " ∂(f(g))/∂g * ∂⁴g/∂p₀³∂p₁");
        referenceRules.put("∂⁴(f(g))/∂p₀³∂p₂",      referenceRules.get("∂⁴(f(g))/∂p₀³∂p₁").replaceAll("p₁", "p₂"));
        referenceRules.put("∂⁴(f(g))/∂p₀³∂p₃",      referenceRules.get("∂⁴(f(g))/∂p₀³∂p₁").replaceAll("p₁", "p₃"));
        referenceRules.put("∂⁴(f(g))/∂p₀∂p₁³",      "∂⁴(f(g))/∂g⁴ * ∂g/∂p₀ * ∂g/∂p₁ * ∂g/∂p₁ * ∂g/∂p₁ +" +
                                                  " 3 * ∂³(f(g))/∂g³ * ∂g/∂p₁ * ∂g/∂p₁ * ∂²g/∂p₀∂p₁ +" +
                                                  " 3 * ∂³(f(g))/∂g³ * ∂g/∂p₀ * ∂g/∂p₁ * ∂²g/∂p₁² +" +
                                                  " 3 * ∂²(f(g))/∂g² * ∂²g/∂p₀∂p₁ * ∂²g/∂p₁² +" +
                                                  " 3 * ∂²(f(g))/∂g² * ∂g/∂p₁ * ∂³g/∂p₀∂p₁² +" +
                                                  " ∂²(f(g))/∂g² * ∂g/∂p₀ * ∂³g/∂p₁³ +" +
                                                  " ∂(f(g))/∂g * ∂⁴g/∂p₀∂p₁³");
        referenceRules.put("∂⁴(f(g))/∂p₀∂p₂³",      referenceRules.get("∂⁴(f(g))/∂p₀∂p₁³").replaceAll("p₁", "p₂"));
        referenceRules.put("∂⁴(f(g))/∂p₀∂p₃³",      referenceRules.get("∂⁴(f(g))/∂p₀∂p₁³").replaceAll("p₁", "p₃"));
        referenceRules.put("∂⁴(f(g))/∂p₁³∂p₂",      referenceRules.get("∂⁴(f(g))/∂p₀³∂p₂").replaceAll("p₀", "p₁"));
        referenceRules.put("∂⁴(f(g))/∂p₁³∂p₃",      referenceRules.get("∂⁴(f(g))/∂p₀³∂p₃").replaceAll("p₀", "p₁"));
        referenceRules.put("∂⁴(f(g))/∂p₁∂p₂³",      referenceRules.get("∂⁴(f(g))/∂p₀∂p₂³").replaceAll("p₀", "p₁"));
        referenceRules.put("∂⁴(f(g))/∂p₁∂p₃³",      referenceRules.get("∂⁴(f(g))/∂p₀∂p₃³").replaceAll("p₀", "p₁"));
        referenceRules.put("∂⁴(f(g))/∂p₂³∂p₃",      referenceRules.get("∂⁴(f(g))/∂p₀³∂p₃").replaceAll("p₀", "p₂"));
        referenceRules.put("∂⁴(f(g))/∂p₂∂p₃³",      referenceRules.get("∂⁴(f(g))/∂p₀∂p₃³").replaceAll("p₀", "p₂"));
        referenceRules.put("∂⁴(f(g))/∂p₀²∂p₁²",     "∂⁴(f(g))/∂g⁴ * ∂g/∂p₀ * ∂g/∂p₀ * ∂g/∂p₁ * ∂g/∂p₁ +" +
                                                  " 4 * ∂³(f(g))/∂g³ * ∂g/∂p₀ * ∂g/∂p₁ * ∂²g/∂p₀∂p₁ +" +
                                                  " ∂³(f(g))/∂g³ * ∂g/∂p₀ * ∂g/∂p₀ * ∂²g/∂p₁² +" +
                                                  " 2 * ∂²(f(g))/∂g² * ∂²g/∂p₀∂p₁ * ∂²g/∂p₀∂p₁ +" +
                                                  " 2 * ∂²(f(g))/∂g² * ∂g/∂p₀ * ∂³g/∂p₀∂p₁² +" +
                                                  " ∂³(f(g))/∂g³ * ∂²g/∂p₀² * ∂g/∂p₁ * ∂g/∂p₁ +" +
                                                  " 2 * ∂²(f(g))/∂g² * ∂g/∂p₁ * ∂³g/∂p₀²∂p₁ +" +
                                                  " ∂²(f(g))/∂g² * ∂²g/∂p₀² * ∂²g/∂p₁² +" +
                                                  " ∂(f(g))/∂g * ∂⁴g/∂p₀²∂p₁²");
        referenceRules.put("∂⁴(f(g))/∂p₀²∂p₂²",     referenceRules.get("∂⁴(f(g))/∂p₀²∂p₁²").replaceAll("p₁", "p₂"));
        referenceRules.put("∂⁴(f(g))/∂p₀²∂p₃²",     referenceRules.get("∂⁴(f(g))/∂p₀²∂p₁²").replaceAll("p₁", "p₃"));
        referenceRules.put("∂⁴(f(g))/∂p₁²∂p₂²",     referenceRules.get("∂⁴(f(g))/∂p₀²∂p₂²").replaceAll("p₀", "p₁"));
        referenceRules.put("∂⁴(f(g))/∂p₁²∂p₃²",     referenceRules.get("∂⁴(f(g))/∂p₀²∂p₃²").replaceAll("p₀", "p₁"));
        referenceRules.put("∂⁴(f(g))/∂p₂²∂p₃²",     referenceRules.get("∂⁴(f(g))/∂p₀²∂p₃²").replaceAll("p₀", "p₂"));

        referenceRules.put("∂⁴(f(g))/∂p₀²∂p₁∂p₂",    "∂⁴(f(g))/∂g⁴ * ∂g/∂p₀ * ∂g/∂p₀ * ∂g/∂p₁ * ∂g/∂p₂ +" +
                                                  " 2 * ∂³(f(g))/∂g³ * ∂g/∂p₀ * ∂g/∂p₁ * ∂²g/∂p₀∂p₂ +" +
                                                  " ∂³(f(g))/∂g³ * ∂g/∂p₀ * ∂g/∂p₀ * ∂²g/∂p₁∂p₂ +" +
                                                  " 2 * ∂³(f(g))/∂g³ * ∂g/∂p₀ * ∂²g/∂p₀∂p₁ * ∂g/∂p₂ +" +
                                                  " 2 * ∂²(f(g))/∂g² * ∂²g/∂p₀∂p₁ * ∂²g/∂p₀∂p₂ +" +
                                                  " 2 * ∂²(f(g))/∂g² * ∂g/∂p₀ * ∂³g/∂p₀∂p₁∂p₂ +" +
                                                  " ∂³(f(g))/∂g³ * ∂²g/∂p₀² * ∂g/∂p₁ * ∂g/∂p₂ +" +
                                                  " ∂²(f(g))/∂g² * ∂g/∂p₁ * ∂³g/∂p₀²∂p₂ +" +
                                                  " ∂²(f(g))/∂g² * ∂²g/∂p₀² * ∂²g/∂p₁∂p₂ +" +
                                                  " ∂²(f(g))/∂g² * ∂³g/∂p₀²∂p₁ * ∂g/∂p₂ +" +
                                                  " ∂(f(g))/∂g * ∂⁴g/∂p₀²∂p₁∂p₂");
        referenceRules.put("∂⁴(f(g))/∂p₀²∂p₁∂p₃",    referenceRules.get("∂⁴(f(g))/∂p₀²∂p₁∂p₂").replaceAll("p₂", "p₃"));
        referenceRules.put("∂⁴(f(g))/∂p₀²∂p₂∂p₃",    referenceRules.get("∂⁴(f(g))/∂p₀²∂p₁∂p₃").replaceAll("p₁", "p₂"));
        referenceRules.put("∂⁴(f(g))/∂p₀∂p₁²∂p₂",    "∂⁴(f(g))/∂g⁴ * ∂g/∂p₀ * ∂g/∂p₁ * ∂g/∂p₁ * ∂g/∂p₂ +" +
                                                  " ∂³(f(g))/∂g³ * ∂g/∂p₁ * ∂g/∂p₁ * ∂²g/∂p₀∂p₂ +" +
                                                  " 2 * ∂³(f(g))/∂g³ * ∂g/∂p₀ * ∂g/∂p₁ * ∂²g/∂p₁∂p₂ +" +
                                                  " 2 * ∂³(f(g))/∂g³ * ∂g/∂p₁ * ∂²g/∂p₀∂p₁ * ∂g/∂p₂ +" +
                                                  " 2 * ∂²(f(g))/∂g² * ∂²g/∂p₀∂p₁ * ∂²g/∂p₁∂p₂ +" +
                                                  " 2 * ∂²(f(g))/∂g² * ∂g/∂p₁ * ∂³g/∂p₀∂p₁∂p₂ +" +
                                                  " ∂³(f(g))/∂g³ * ∂g/∂p₀ * ∂²g/∂p₁² * ∂g/∂p₂ +" +
                                                  " ∂²(f(g))/∂g² * ∂²g/∂p₁² * ∂²g/∂p₀∂p₂ +" +
                                                  " ∂²(f(g))/∂g² * ∂g/∂p₀ * ∂³g/∂p₁²∂p₂ +" +
                                                  " ∂²(f(g))/∂g² * ∂³g/∂p₀∂p₁² * ∂g/∂p₂ +" +
                                                  " ∂(f(g))/∂g * ∂⁴g/∂p₀∂p₁²∂p₂");
        referenceRules.put("∂⁴(f(g))/∂p₀∂p₁²∂p₃",    referenceRules.get("∂⁴(f(g))/∂p₀∂p₁²∂p₂").replaceAll("p₂", "p₃"));
        referenceRules.put("∂⁴(f(g))/∂p₁²∂p₂∂p₃",    referenceRules.get("∂⁴(f(g))/∂p₀²∂p₂∂p₃").replaceAll("p₀", "p₁"));
        referenceRules.put("∂⁴(f(g))/∂p₀∂p₁∂p₂²",    "∂⁴(f(g))/∂g⁴ * ∂g/∂p₀ * ∂g/∂p₁ * ∂g/∂p₂ * ∂g/∂p₂ +" +
                                                  " 2 * ∂³(f(g))/∂g³ * ∂g/∂p₁ * ∂g/∂p₂ * ∂²g/∂p₀∂p₂ +" +
                                                  " 2 * ∂³(f(g))/∂g³ * ∂g/∂p₀ * ∂g/∂p₂ * ∂²g/∂p₁∂p₂ +" +
                                                  " ∂³(f(g))/∂g³ * ∂g/∂p₀ * ∂g/∂p₁ * ∂²g/∂p₂² +" +
                                                  " 2 * ∂²(f(g))/∂g² * ∂²g/∂p₀∂p₂ * ∂²g/∂p₁∂p₂ +" +
                                                  " ∂²(f(g))/∂g² * ∂g/∂p₁ * ∂³g/∂p₀∂p₂² +" +
                                                  " ∂²(f(g))/∂g² * ∂g/∂p₀ * ∂³g/∂p₁∂p₂² +" +
                                                  " ∂³(f(g))/∂g³ * ∂²g/∂p₀∂p₁ * ∂g/∂p₂ * ∂g/∂p₂ +" +
                                                  " 2 * ∂²(f(g))/∂g² * ∂g/∂p₂ * ∂³g/∂p₀∂p₁∂p₂ +" +
                                                  " ∂²(f(g))/∂g² * ∂²g/∂p₀∂p₁ * ∂²g/∂p₂² +" +
                                                  " ∂(f(g))/∂g * ∂⁴g/∂p₀∂p₁∂p₂²");
        referenceRules.put("∂⁴(f(g))/∂p₀∂p₂²∂p₃",    referenceRules.get("∂⁴(f(g))/∂p₀∂p₁²∂p₃").replaceAll("p₁", "p₂"));
        referenceRules.put("∂⁴(f(g))/∂p₁∂p₂²∂p₃",    referenceRules.get("∂⁴(f(g))/∂p₀∂p₂²∂p₃").replaceAll("p₀", "p₁"));
        referenceRules.put("∂⁴(f(g))/∂p₀∂p₁∂p₃²",    referenceRules.get("∂⁴(f(g))/∂p₀∂p₁∂p₂²").replaceAll("p₂", "p₃"));
        referenceRules.put("∂⁴(f(g))/∂p₀∂p₂∂p₃²",    referenceRules.get("∂⁴(f(g))/∂p₀∂p₁∂p₃²").replaceAll("p₁", "p₂"));
        referenceRules.put("∂⁴(f(g))/∂p₁∂p₂∂p₃²",    referenceRules.get("∂⁴(f(g))/∂p₀∂p₂∂p₃²").replaceAll("p₀", "p₁"));
        referenceRules.put("∂⁴(f(g))/∂p₀∂p₁∂p₂∂p₃",   "∂⁴(f(g))/∂g⁴ * ∂g/∂p₀ * ∂g/∂p₁ * ∂g/∂p₂ * ∂g/∂p₃ +" +
                                                  " ∂³(f(g))/∂g³ * ∂g/∂p₁ * ∂g/∂p₂ * ∂²g/∂p₀∂p₃ +" +
                                                  " ∂³(f(g))/∂g³ * ∂g/∂p₀ * ∂g/∂p₂ * ∂²g/∂p₁∂p₃ +" +
                                                  " ∂³(f(g))/∂g³ * ∂g/∂p₀ * ∂g/∂p₁ * ∂²g/∂p₂∂p₃ +" +
                                                  " ∂³(f(g))/∂g³ * ∂g/∂p₁ * ∂²g/∂p₀∂p₂ * ∂g/∂p₃ +" +
                                                  " ∂²(f(g))/∂g² * ∂²g/∂p₀∂p₂ * ∂²g/∂p₁∂p₃ +" +
                                                  " ∂²(f(g))/∂g² * ∂g/∂p₁ * ∂³g/∂p₀∂p₂∂p₃ +" +
                                                  " ∂³(f(g))/∂g³ * ∂g/∂p₀ * ∂²g/∂p₁∂p₂ * ∂g/∂p₃ +" +
                                                  " ∂²(f(g))/∂g² * ∂²g/∂p₁∂p₂ * ∂²g/∂p₀∂p₃ +" +
                                                  " ∂²(f(g))/∂g² * ∂g/∂p₀ * ∂³g/∂p₁∂p₂∂p₃ +" +
                                                  " ∂³(f(g))/∂g³ * ∂²g/∂p₀∂p₁ * ∂g/∂p₂ * ∂g/∂p₃ +" +
                                                  " ∂²(f(g))/∂g² * ∂g/∂p₂ * ∂³g/∂p₀∂p₁∂p₃ +" +
                                                  " ∂²(f(g))/∂g² * ∂²g/∂p₀∂p₁ * ∂²g/∂p₂∂p₃ +" +
                                                  " ∂²(f(g))/∂g² * ∂³g/∂p₀∂p₁∂p₂ * ∂g/∂p₃ +" +
                                                  " ∂(f(g))/∂g * ∂⁴g/∂p₀∂p₁∂p₂∂p₃");

        Field compFieldArrayField = DSCompiler.class.getDeclaredField("compIndirection");
        compFieldArrayField.setAccessible(true);

        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                DSCompiler compiler = DSCompiler.getCompiler(i, j);
                Object[][] compIndirection = (Object[][]) compFieldArrayField.get(compiler);
                for (int k = 0; k < compIndirection.length; ++k) {
                    String product = ordersToString(compiler.getPartialDerivativeOrders(k),
                                                    "(f(g))", variables("p"));
                    String rule = univariateCompositionMappersToString(compiler, compIndirection[k]);
                    Assert.assertEquals(product, referenceRules.get(product), rule.toString());
                }
            }
        }

    }

    @Test
    public void testRebaserRules()
        throws IllegalAccessException, IllegalArgumentException, InvocationTargetException,
               NoSuchMethodException, SecurityException, NoSuchFieldException {

        // the following reference rules have all been computed independently from the library,
        // using only pencil and paper (which was really tedious) and using search and replace to handle symmetries
        Map<String,String> referenceRules = new HashMap<String, String>();
        referenceRules.put("f",              "f");
        referenceRules.put("∂f/∂q₀",         "∂f/∂p₀ ∂p₀/∂q₀ + ∂f/∂p₁ ∂p₁/∂q₀");
        referenceRules.put("∂f/∂q₁",         referenceRules.get("∂f/∂q₀").replaceAll("q₀", "q₁"));
        referenceRules.put("∂f/∂q₂",         referenceRules.get("∂f/∂q₀").replaceAll("q₀", "q₂"));
        referenceRules.put("∂²f/∂q₀²",       "∂²f/∂p₀² (∂p₀/∂q₀)² + 2 ∂²f/∂p₀∂p₁ ∂p₀/∂q₀ ∂p₁/∂q₀" +
                                             " + ∂f/∂p₀ ∂²p₀/∂q₀² + ∂²f/∂p₁² (∂p₁/∂q₀)² + ∂f/∂p₁ ∂²p₁/∂q₀²");
        referenceRules.put("∂²f/∂q₁²",       referenceRules.get("∂²f/∂q₀²").replaceAll("q₀", "q₁"));
        referenceRules.put("∂²f/∂q₂²",       referenceRules.get("∂²f/∂q₀²").replaceAll("q₀", "q₂"));
        referenceRules.put("∂²f/∂q₀∂q₁",     "∂²f/∂p₀² ∂p₀/∂q₀ ∂p₀/∂q₁ + ∂²f/∂p₀∂p₁ ∂p₀/∂q₁ ∂p₁/∂q₀ + ∂f/∂p₀ ∂²p₀/∂q₀∂q₁" +
                                             " + ∂²f/∂p₀∂p₁ ∂p₀/∂q₀ ∂p₁/∂q₁ + ∂²f/∂p₁² ∂p₁/∂q₀ ∂p₁/∂q₁ + ∂f/∂p₁ ∂²p₁/∂q₀∂q₁");
        referenceRules.put("∂²f/∂q₀∂q₂",     referenceRules.get("∂²f/∂q₀∂q₁").replaceAll("q₁", "q₂"));
        referenceRules.put("∂²f/∂q₁∂q₂",     referenceRules.get("∂²f/∂q₀∂q₂").replaceAll("q₀", "q₁"));
        referenceRules.put("∂³f/∂q₀³",       "∂³f/∂p₀³ (∂p₀/∂q₀)³ + 3 ∂³f/∂p₀²∂p₁ (∂p₀/∂q₀)² ∂p₁/∂q₀" +
                                             " + 3 ∂²f/∂p₀² ∂p₀/∂q₀ ∂²p₀/∂q₀² + 3 ∂³f/∂p₀∂p₁² ∂p₀/∂q₀ (∂p₁/∂q₀)²" +
                                             " + 3 ∂²f/∂p₀∂p₁ ∂²p₀/∂q₀² ∂p₁/∂q₀ + 3 ∂²f/∂p₀∂p₁ ∂p₀/∂q₀ ∂²p₁/∂q₀²" +
                                             " + ∂f/∂p₀ ∂³p₀/∂q₀³ + ∂³f/∂p₁³ (∂p₁/∂q₀)³" +
                                             " + 3 ∂²f/∂p₁² ∂p₁/∂q₀ ∂²p₁/∂q₀² + ∂f/∂p₁ ∂³p₁/∂q₀³");
        referenceRules.put("∂³f/∂q₁³",       referenceRules.get("∂³f/∂q₀³").replaceAll("q₀", "q₁"));
        referenceRules.put("∂³f/∂q₂³",       referenceRules.get("∂³f/∂q₀³").replaceAll("q₀", "q₂"));
        referenceRules.put("∂³f/∂q₀²∂q₁",    "∂³f/∂p₀³ (∂p₀/∂q₀)² ∂p₀/∂q₁ + 2 ∂³f/∂p₀²∂p₁ ∂p₀/∂q₀ ∂p₀/∂q₁ ∂p₁/∂q₀" +
                                             " + ∂²f/∂p₀² ∂²p₀/∂q₀² ∂p₀/∂q₁ + 2 ∂²f/∂p₀² ∂p₀/∂q₀ ∂²p₀/∂q₀∂q₁" +
                                             " + ∂³f/∂p₀∂p₁² ∂p₀/∂q₁ (∂p₁/∂q₀)² + 2 ∂²f/∂p₀∂p₁ ∂²p₀/∂q₀∂q₁ ∂p₁/∂q₀" +
                                             " + ∂²f/∂p₀∂p₁ ∂p₀/∂q₁ ∂²p₁/∂q₀² + ∂f/∂p₀ ∂³p₀/∂q₀²∂q₁" +
                                             " + ∂³f/∂p₀²∂p₁ (∂p₀/∂q₀)² ∂p₁/∂q₁ + 2 ∂³f/∂p₀∂p₁² ∂p₀/∂q₀ ∂p₁/∂q₀ ∂p₁/∂q₁" +
                                             " + ∂²f/∂p₀∂p₁ ∂²p₀/∂q₀² ∂p₁/∂q₁ + 2 ∂²f/∂p₀∂p₁ ∂p₀/∂q₀ ∂²p₁/∂q₀∂q₁" +
                                             " + ∂³f/∂p₁³ (∂p₁/∂q₀)² ∂p₁/∂q₁ + ∂²f/∂p₁² ∂²p₁/∂q₀² ∂p₁/∂q₁" +
                                             " + 2 ∂²f/∂p₁² ∂p₁/∂q₀ ∂²p₁/∂q₀∂q₁ + ∂f/∂p₁ ∂³p₁/∂q₀²∂q₁");
        referenceRules.put("∂³f/∂q₀²∂q₂",    referenceRules.get("∂³f/∂q₀²∂q₁").replaceAll("q₁", "q₂"));
        referenceRules.put("∂³f/∂q₁²∂q₂",    referenceRules.get("∂³f/∂q₀²∂q₂").replaceAll("q₀", "q₁"));
        referenceRules.put("∂³f/∂q₁²∂q₀",    referenceRules.get("∂³f/∂q₁²∂q₂").replaceAll("q₂", "q₀"));
        referenceRules.put("∂³f/∂q₀∂q₁²",    "∂³f/∂p₀³ ∂p₀/∂q₀ (∂p₀/∂q₁)² + ∂³f/∂p₀²∂p₁ (∂p₀/∂q₁)² ∂p₁/∂q₀" +
                                             " + 2 ∂²f/∂p₀² ∂p₀/∂q₁ ∂²p₀/∂q₀∂q₁ + 2 ∂³f/∂p₀²∂p₁ ∂p₀/∂q₀ ∂p₀/∂q₁ ∂p₁/∂q₁" +
                                             " + 2 ∂³f/∂p₀∂p₁² ∂p₀/∂q₁ ∂p₁/∂q₀ ∂p₁/∂q₁ + 2 ∂²f/∂p₀∂p₁ ∂²p₀/∂q₀∂q₁ ∂p₁/∂q₁" +
                                             " + 2 ∂²f/∂p₀∂p₁ ∂p₀/∂q₁ ∂²p₁/∂q₀∂q₁ + ∂²f/∂p₀² ∂p₀/∂q₀ ∂²p₀/∂q₁²" +
                                             " + ∂²f/∂p₀∂p₁ ∂²p₀/∂q₁² ∂p₁/∂q₀ + ∂f/∂p₀ ∂³p₀/∂q₀∂q₁²" +
                                             " + ∂³f/∂p₀∂p₁² ∂p₀/∂q₀ (∂p₁/∂q₁)² + ∂³f/∂p₁³ ∂p₁/∂q₀ (∂p₁/∂q₁)²" +
                                             " + 2 ∂²f/∂p₁² ∂p₁/∂q₁ ∂²p₁/∂q₀∂q₁ + ∂²f/∂p₀∂p₁ ∂p₀/∂q₀ ∂²p₁/∂q₁²" +
                                             " + ∂²f/∂p₁² ∂p₁/∂q₀ ∂²p₁/∂q₁²" +
                                             " + ∂f/∂p₁ ∂³p₁/∂q₀∂q₁²");
        referenceRules.put("∂³f/∂q₀∂q₂²",   referenceRules.get("∂³f/∂q₀∂q₁²").replaceAll("q₁", "q₂"));
        referenceRules.put("∂³f/∂q₁∂q₂²",   referenceRules.get("∂³f/∂q₀∂q₂²").replaceAll("q₀", "q₁"));
        referenceRules.put("∂³f/∂q₀∂q₁∂q₂", "∂³f/∂p₀³ ∂p₀/∂q₀ ∂p₀/∂q₁ ∂p₀/∂q₂ + ∂³f/∂p₀²∂p₁ ∂p₀/∂q₁ ∂p₀/∂q₂ ∂p₁/∂q₀" +
                                            " + ∂²f/∂p₀² ∂²p₀/∂q₀∂q₁ ∂p₀/∂q₂ + ∂²f/∂p₀² ∂p₀/∂q₁ ∂²p₀/∂q₀∂q₂" +
                                            " + ∂³f/∂p₀²∂p₁ ∂p₀/∂q₀ ∂p₀/∂q₂ ∂p₁/∂q₁ + ∂³f/∂p₀∂p₁² ∂p₀/∂q₂ ∂p₁/∂q₀ ∂p₁/∂q₁" +
                                            " + ∂²f/∂p₀∂p₁ ∂²p₀/∂q₀∂q₂ ∂p₁/∂q₁ + ∂²f/∂p₀∂p₁ ∂p₀/∂q₂ ∂²p₁/∂q₀∂q₁" +
                                            " + ∂²f/∂p₀² ∂p₀/∂q₀ ∂²p₀/∂q₁∂q₂ + ∂²f/∂p₀∂p₁ ∂²p₀/∂q₁∂q₂ ∂p₁/∂q₀" +
                                            " + ∂f/∂p₀ ∂³p₀/∂q₀∂q₁∂q₂ + ∂³f/∂p₀²∂p₁ ∂p₀/∂q₀ ∂p₀/∂q₁ ∂p₁/∂q₂" +
                                            " + ∂³f/∂p₀∂p₁² ∂p₀/∂q₁ ∂p₁/∂q₀ ∂p₁/∂q₂ + ∂²f/∂p₀∂p₁ ∂²p₀/∂q₀∂q₁ ∂p₁/∂q₂" +
                                            " + ∂²f/∂p₀∂p₁ ∂p₀/∂q₁ ∂²p₁/∂q₀∂q₂ + ∂³f/∂p₀∂p₁² ∂p₀/∂q₀ ∂p₁/∂q₁ ∂p₁/∂q₂" +
                                            " + ∂³f/∂p₁³ ∂p₁/∂q₀ ∂p₁/∂q₁ ∂p₁/∂q₂ + ∂²f/∂p₁² ∂²p₁/∂q₀∂q₁ ∂p₁/∂q₂" +
                                            " + ∂²f/∂p₁² ∂p₁/∂q₁ ∂²p₁/∂q₀∂q₂ + ∂²f/∂p₀∂p₁ ∂p₀/∂q₀ ∂²p₁/∂q₁∂q₂" +
                                            " + ∂²f/∂p₁² ∂p₁/∂q₀ ∂²p₁/∂q₁∂q₂ + ∂f/∂p₁ ∂³p₁/∂q₀∂q₁∂q₂");

        Method getterMethod = DSCompiler.class.getDeclaredMethod("getRebaser", DSCompiler.class);
        getterMethod.setAccessible(true);

        for (int order = 0; order < 4; ++order) {

            // assuming f = f(p₀, p₁)
            //          p₀ = p₀(q₀, q₁, q₂)
            //          p₁ = p₁(q₀, q₁, q₂)
            DSCompiler c2 = DSCompiler.getCompiler(2, order);
            DSCompiler c3 = DSCompiler.getCompiler(3, order);
            Object[][] rebaser = (Object[][]) getterMethod.invoke(c2, c3);

            Assert.assertEquals(c3.getSize(), rebaser.length);
            for (int k = 0; k < rebaser.length; ++k) {
                String key  = ordersToString(c3.getPartialDerivativeOrders(k), "f", variables("q"));
                String rule = multivariateCompositionMappersToString(c2, c3, rebaser[k]);
                Assert.assertEquals(referenceRules.get(key), rule);
            }

        }

    }

    private void checkIndices(int[] indices, int ... expected) {
        Assert.assertEquals(expected.length, indices.length);
        for (int i = 0; i < expected.length; ++i) {
            Assert.assertEquals(expected[i], indices[i]);
        }
    }

    private String orderToString(int order, String functionName, String parameterName) {
        if (order == 0) {
            return functionName;
        } else if (order == 1) {
            return "∂" + functionName + "/∂" + parameterName;
        } else {
            return "∂" + exponent(order) + functionName + "/∂" + parameterName + exponent(order);
        }
    }

    private String ordersToString(int[] orders, String functionName, String ... parametersNames) {

        int sumOrders = 0;
        for (int order : orders) {
            sumOrders += order;
        }

        if (sumOrders == 0) {
            return functionName;
        }

        StringBuilder builder = new StringBuilder();
        builder.append('∂');
        if (sumOrders > 1) {
            builder.append(exponent(sumOrders));
        }
        builder.append(functionName).append('/');
        for (int i = 0; i < orders.length; ++i) {
            if (orders[i] > 0) {
                builder.append('∂').append(parametersNames[i]);
                if (orders[i] > 1) {
                    builder.append(exponent(orders[i]));
                }
            }
        }
        return builder.toString();

    }

    private String univariateCompositionMappersToString(final DSCompiler compiler,  final Object[] mappers) {
         try {
             
             Class<?> abstractMapperClass = Stream.
                             of(DSCompiler.class.getDeclaredClasses()).
                             filter(c -> c.getName().endsWith("AbstractMapper")).
                             findAny().
                             get();
             Class<?> univariateCompositionMapperClass = Stream.
                             of(DSCompiler.class.getDeclaredClasses()).
                             filter(c -> c.getName().endsWith("UnivariateCompositionMapper")).
                             findAny().
                             get();
             Method coeffMethod = abstractMapperClass.getDeclaredMethod("getCoeff");
             Field fIndexField = univariateCompositionMapperClass.getDeclaredField("fIndex");
             fIndexField.setAccessible(true);
             Field dsIndicesField = univariateCompositionMapperClass.getDeclaredField("dsIndices");
             dsIndicesField.setAccessible(true);

             StringBuilder rule = new StringBuilder();
             for (Object term : mappers) {
                 if (rule.length() > 0) {
                     rule.append(" + ");
                 }
                 if (((Integer) coeffMethod.invoke(term)).intValue() > 1) {
                     rule.append(((Integer) coeffMethod.invoke(term)).intValue()).append(" * ");
                 }
                 rule.append(orderToString(((Integer) fIndexField.get(term)).intValue(), "(f(g))", "g"));
                 int[] dsIndex = (int[]) dsIndicesField.get(term);
                 for (int l = 0; l < dsIndex.length; ++l) {
                     rule.append(" * ");
                     rule.append(ordersToString(compiler.getPartialDerivativeOrders(dsIndex[l]),
                                                "g", "p₀", "p₁", "p₂", "p₃"));
                 }
             }
             return rule.toString();

         } catch (NoSuchMethodException | SecurityException | NoSuchFieldException | IllegalAccessException |
                  IllegalArgumentException | InvocationTargetException e) {
             Assert.fail(e.getLocalizedMessage());
             return null;
         }
     }

    private String multivariateCompositionMappersToString(final DSCompiler compiler, final DSCompiler baseCompiler,
                                                          final Object[] mappers) {
         try {
             Class<?> abstractMapperClass = Stream.
                             of(DSCompiler.class.getDeclaredClasses()).
                             filter(c -> c.getName().endsWith("AbstractMapper")).
                             findAny().
                             get();
             Class<?> multivariateCompositionMapperClass = Stream.
                             of(DSCompiler.class.getDeclaredClasses()).
                             filter(c -> c.getName().endsWith("MultivariateCompositionMapper")).
                             findAny().
                             get();
             Method coeffMethod = abstractMapperClass.getDeclaredMethod("getCoeff");
             Field dsIndexField = multivariateCompositionMapperClass.getDeclaredField("dsIndex");
             dsIndexField.setAccessible(true);
             Field productIndicesField = multivariateCompositionMapperClass.getDeclaredField("productIndices");
             productIndicesField.setAccessible(true);

             StringBuilder rule = new StringBuilder();
             for (int i = 0; i < mappers.length; ++i) {
                 if (i > 0) {
                     rule.append(" + ");
                 }
                 final int coeff = ((Integer) coeffMethod.invoke(mappers[i])).intValue();
                 if (coeff > 1) {
                     rule.append(coeff);
                     rule.append(' ');
                 }
                 final int dsIndex = dsIndexField.getInt(mappers[i]);
                 rule.append(ordersToString(compiler.getPartialDerivativeOrders(dsIndex),
                                            "f", variables("p")));
                 final int[] productIndices = (int[]) productIndicesField.get(mappers[i]);
                 int j = 0;
                 while (j < productIndices.length) {
                     int count = 1;
                     while (j + count < productIndices.length && productIndices[j + count] == productIndices[j]) {
                         ++count;
                     }
                     final int varIndex   = productIndices[j] / baseCompiler.getSize();
                     final int varDSIndex = productIndices[j] % baseCompiler.getSize();
                     rule.append(' ');
                     if (count > 1) {
                         rule.append('(');
                     }
                     rule.append(ordersToString(baseCompiler.getPartialDerivativeOrders(varDSIndex),
                                                variables("p")[varIndex], variables("q")));
                     if (count > 1) {
                         rule.append(')');
                         rule.append(exponent(count));
                     }
                     j += count;
                 }
             }
             return rule.toString();

         } catch (NoSuchMethodException | SecurityException | NoSuchFieldException | IllegalAccessException |
                  IllegalArgumentException | InvocationTargetException e) {
             Assert.fail(e.getLocalizedMessage());
             return null;
         }
     }

     private String[] variables(final String baseName) {
         return new String[] {
             baseName + "₀", baseName + "₁", baseName + "₂", baseName + "₃", baseName + "₄",
             baseName + "₅", baseName + "₆", baseName + "₇", baseName + "₈", baseName + "₉"
         };
     }

    private String exponent(int e) {
        switch (e) {
            case 0 : return "";
            case 1 : return "";
            case 2 : return "²";
            case 3 : return "³";
            case 4 : return "⁴";
            case 5 : return "⁵";
            case 6 : return "⁶";
            case 7 : return "⁷";
            case 8 : return "⁸";
            case 9 : return "⁹";
            default:
                Assert.fail("exponent out of range");
                return null;
        }
    }

}
