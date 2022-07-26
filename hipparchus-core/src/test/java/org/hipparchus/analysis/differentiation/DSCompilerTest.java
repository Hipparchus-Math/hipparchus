/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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
        referenceRules.put("(f*g)",          "f * g");
        referenceRules.put("∂(f*g)/∂x",      "f * ∂g/∂x + ∂f/∂x * g");
        referenceRules.put("∂(f*g)/∂y",      referenceRules.get("∂(f*g)/∂x").replaceAll("x", "y"));
        referenceRules.put("∂(f*g)/∂z",      referenceRules.get("∂(f*g)/∂x").replaceAll("x", "z"));
        referenceRules.put("∂(f*g)/∂t",      referenceRules.get("∂(f*g)/∂x").replaceAll("x", "t"));
        referenceRules.put("∂²(f*g)/∂x²",    "f * ∂²g/∂x² + 2 * ∂f/∂x * ∂g/∂x + ∂²f/∂x² * g");
        referenceRules.put("∂²(f*g)/∂y²",    referenceRules.get("∂²(f*g)/∂x²").replaceAll("x", "y"));
        referenceRules.put("∂²(f*g)/∂z²",    referenceRules.get("∂²(f*g)/∂x²").replaceAll("x", "z"));
        referenceRules.put("∂²(f*g)/∂t²",    referenceRules.get("∂²(f*g)/∂x²").replaceAll("x", "t"));
        referenceRules.put("∂²(f*g)/∂x∂y",   "f * ∂²g/∂x∂y + ∂f/∂y * ∂g/∂x + ∂f/∂x * ∂g/∂y + ∂²f/∂x∂y * g");
        referenceRules.put("∂²(f*g)/∂x∂z",   referenceRules.get("∂²(f*g)/∂x∂y").replaceAll("y", "z"));
        referenceRules.put("∂²(f*g)/∂x∂t",   referenceRules.get("∂²(f*g)/∂x∂y").replaceAll("y", "t"));
        referenceRules.put("∂²(f*g)/∂y∂z",   referenceRules.get("∂²(f*g)/∂x∂z").replaceAll("x", "y"));
        referenceRules.put("∂²(f*g)/∂y∂t",   referenceRules.get("∂²(f*g)/∂x∂t").replaceAll("x", "y"));
        referenceRules.put("∂²(f*g)/∂z∂t",   referenceRules.get("∂²(f*g)/∂x∂t").replaceAll("x", "z"));
        referenceRules.put("∂³(f*g)/∂x³",    "f * ∂³g/∂x³ +" +
                                             " 3 * ∂f/∂x * ∂²g/∂x² +" +
                                             " 3 * ∂²f/∂x² * ∂g/∂x +" +
                                             " ∂³f/∂x³ * g");
        referenceRules.put("∂³(f*g)/∂y³",   referenceRules.get("∂³(f*g)/∂x³").replaceAll("x", "y"));
        referenceRules.put("∂³(f*g)/∂z³",   referenceRules.get("∂³(f*g)/∂x³").replaceAll("x", "z"));
        referenceRules.put("∂³(f*g)/∂t³",   referenceRules.get("∂³(f*g)/∂x³").replaceAll("x", "t"));
        referenceRules.put("∂³(f*g)/∂x²∂y",  "f * ∂³g/∂x²∂y +" +
                                             " ∂f/∂y * ∂²g/∂x² +" +
                                             " 2 * ∂f/∂x * ∂²g/∂x∂y +" +
                                             " 2 * ∂²f/∂x∂y * ∂g/∂x +" +
                                             " ∂²f/∂x² * ∂g/∂y +" +
                                             " ∂³f/∂x²∂y * g");
        referenceRules.put("∂³(f*g)/∂x∂y²",  "f * ∂³g/∂x∂y² +" +
                                             " 2 * ∂f/∂y * ∂²g/∂x∂y +" +
                                             " ∂²f/∂y² * ∂g/∂x +" +
                                             " ∂f/∂x * ∂²g/∂y² +" +
                                             " 2 * ∂²f/∂x∂y * ∂g/∂y +" +
                                             " ∂³f/∂x∂y² * g");
        referenceRules.put("∂³(f*g)/∂x²∂z",   referenceRules.get("∂³(f*g)/∂x²∂y").replaceAll("y", "z"));
        referenceRules.put("∂³(f*g)/∂y²∂z",   referenceRules.get("∂³(f*g)/∂x²∂z").replaceAll("x", "y"));
        referenceRules.put("∂³(f*g)/∂x∂z²",   referenceRules.get("∂³(f*g)/∂x∂y²").replaceAll("y", "z"));
        referenceRules.put("∂³(f*g)/∂y∂z²",   referenceRules.get("∂³(f*g)/∂x∂z²").replaceAll("x", "y"));
        referenceRules.put("∂³(f*g)/∂x²∂t",   referenceRules.get("∂³(f*g)/∂x²∂z").replaceAll("z", "t"));
        referenceRules.put("∂³(f*g)/∂y²∂t",   referenceRules.get("∂³(f*g)/∂x²∂t").replaceAll("x", "y"));
        referenceRules.put("∂³(f*g)/∂z²∂t",   referenceRules.get("∂³(f*g)/∂x²∂t").replaceAll("x", "z"));
        referenceRules.put("∂³(f*g)/∂x∂t²",   referenceRules.get("∂³(f*g)/∂x∂y²").replaceAll("y", "t"));
        referenceRules.put("∂³(f*g)/∂y∂t²",   referenceRules.get("∂³(f*g)/∂x∂t²").replaceAll("x", "y"));
        referenceRules.put("∂³(f*g)/∂z∂t²",   referenceRules.get("∂³(f*g)/∂x∂t²").replaceAll("x", "z"));
        referenceRules.put("∂³(f*g)/∂x∂y∂z", "f * ∂³g/∂x∂y∂z +" +
                                             " ∂f/∂z * ∂²g/∂x∂y +" +
                                             " ∂f/∂y * ∂²g/∂x∂z +" +
                                             " ∂²f/∂y∂z * ∂g/∂x +" +
                                             " ∂f/∂x * ∂²g/∂y∂z +" +
                                             " ∂²f/∂x∂z * ∂g/∂y +" +
                                             " ∂²f/∂x∂y * ∂g/∂z +" +
                                             " ∂³f/∂x∂y∂z * g");
        referenceRules.put("∂³(f*g)/∂x∂y∂t", referenceRules.get("∂³(f*g)/∂x∂y∂z").replaceAll("z", "t"));
        referenceRules.put("∂³(f*g)/∂x∂z∂t", referenceRules.get("∂³(f*g)/∂x∂y∂t").replaceAll("y", "z"));
        referenceRules.put("∂³(f*g)/∂y∂z∂t", referenceRules.get("∂³(f*g)/∂x∂z∂t").replaceAll("x", "y"));

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
                                                    "(f*g)", "x", "y", "z", "t");
                    StringBuilder rule = new StringBuilder();
                    for (Object term : multIndirection[k]) {
                        if (rule.length() > 0) {
                            rule.append(" + ");
                        }
                        if (((Integer) coeffMethod.invoke(term)).intValue() > 1) {
                            rule.append(((Integer) coeffMethod.invoke(term)).intValue()).append(" * ");
                        }
                        rule.append(ordersToString(compiler.getPartialDerivativeOrders(((Integer) lhsField.get(term)).intValue()),
                                                   "f", "x", "y", "z", "t"));
                        rule.append(" * ");
                        rule.append(ordersToString(compiler.getPartialDerivativeOrders(((Integer) rhsField.get(term)).intValue()),
                                                   "g", "x", "y", "z", "t"));
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
        referenceRules.put("∂(f(g))/∂x",          "∂(f(g))/∂g * ∂g/∂x");
        referenceRules.put("∂(f(g))/∂y",          referenceRules.get("∂(f(g))/∂x").replaceAll("x", "y"));
        referenceRules.put("∂(f(g))/∂z",          referenceRules.get("∂(f(g))/∂x").replaceAll("x", "z"));
        referenceRules.put("∂(f(g))/∂t",          referenceRules.get("∂(f(g))/∂x").replaceAll("x", "t"));
        referenceRules.put("∂²(f(g))/∂x²",        "∂²(f(g))/∂g² * ∂g/∂x * ∂g/∂x + ∂(f(g))/∂g * ∂²g/∂x²");
        referenceRules.put("∂²(f(g))/∂y²",        referenceRules.get("∂²(f(g))/∂x²").replaceAll("x", "y"));
        referenceRules.put("∂²(f(g))/∂z²",        referenceRules.get("∂²(f(g))/∂x²").replaceAll("x", "z"));
        referenceRules.put("∂²(f(g))/∂t²",        referenceRules.get("∂²(f(g))/∂x²").replaceAll("x", "t"));
        referenceRules.put("∂²(f(g))/∂x∂y",       "∂²(f(g))/∂g² * ∂g/∂x * ∂g/∂y + ∂(f(g))/∂g * ∂²g/∂x∂y");
        referenceRules.put("∂²(f(g))/∂x∂z",       referenceRules.get("∂²(f(g))/∂x∂y").replaceAll("y", "z"));
        referenceRules.put("∂²(f(g))/∂x∂t",       referenceRules.get("∂²(f(g))/∂x∂y").replaceAll("y", "t"));
        referenceRules.put("∂²(f(g))/∂y∂z",       referenceRules.get("∂²(f(g))/∂x∂z").replaceAll("x", "y"));
        referenceRules.put("∂²(f(g))/∂y∂t",       referenceRules.get("∂²(f(g))/∂x∂t").replaceAll("x", "y"));
        referenceRules.put("∂²(f(g))/∂z∂t",       referenceRules.get("∂²(f(g))/∂x∂t").replaceAll("x", "z"));
        referenceRules.put("∂³(f(g))/∂x³",        "∂³(f(g))/∂g³ * ∂g/∂x * ∂g/∂x * ∂g/∂x +" +
                                                  " 3 * ∂²(f(g))/∂g² * ∂g/∂x * ∂²g/∂x² +" +
                                                  " ∂(f(g))/∂g * ∂³g/∂x³");
        referenceRules.put("∂³(f(g))/∂y³",        referenceRules.get("∂³(f(g))/∂x³").replaceAll("x", "y"));
        referenceRules.put("∂³(f(g))/∂z³",        referenceRules.get("∂³(f(g))/∂x³").replaceAll("x", "z"));
        referenceRules.put("∂³(f(g))/∂t³",        referenceRules.get("∂³(f(g))/∂x³").replaceAll("x", "t"));
        referenceRules.put("∂³(f(g))/∂x∂y²",      "∂³(f(g))/∂g³ * ∂g/∂x * ∂g/∂y * ∂g/∂y +" +
                                                  " 2 * ∂²(f(g))/∂g² * ∂g/∂y * ∂²g/∂x∂y +" +
                                                  " ∂²(f(g))/∂g² * ∂g/∂x * ∂²g/∂y² +" +
                                                  " ∂(f(g))/∂g * ∂³g/∂x∂y²");
        referenceRules.put("∂³(f(g))/∂x∂z²",      referenceRules.get("∂³(f(g))/∂x∂y²").replaceAll("y", "z"));
        referenceRules.put("∂³(f(g))/∂x∂t²",      referenceRules.get("∂³(f(g))/∂x∂y²").replaceAll("y", "t"));
        referenceRules.put("∂³(f(g))/∂y∂z²",      referenceRules.get("∂³(f(g))/∂x∂z²").replaceAll("x", "y"));
        referenceRules.put("∂³(f(g))/∂y∂t²",      referenceRules.get("∂³(f(g))/∂x∂t²").replaceAll("x", "y"));
        referenceRules.put("∂³(f(g))/∂z∂t²",      referenceRules.get("∂³(f(g))/∂x∂t²").replaceAll("x", "z"));
        referenceRules.put("∂³(f(g))/∂x²∂y",      "∂³(f(g))/∂g³ * ∂g/∂x * ∂g/∂x * ∂g/∂y +" +
                                                  " 2 * ∂²(f(g))/∂g² * ∂g/∂x * ∂²g/∂x∂y +" +
                                                  " ∂²(f(g))/∂g² * ∂²g/∂x² * ∂g/∂y +" +
                                                  " ∂(f(g))/∂g * ∂³g/∂x²∂y");
        referenceRules.put("∂³(f(g))/∂x²∂z",      referenceRules.get("∂³(f(g))/∂x²∂y").replaceAll("y", "z"));
        referenceRules.put("∂³(f(g))/∂x²∂t",      referenceRules.get("∂³(f(g))/∂x²∂y").replaceAll("y", "t"));
        referenceRules.put("∂³(f(g))/∂y²∂z",      referenceRules.get("∂³(f(g))/∂x²∂z").replaceAll("x", "y"));
        referenceRules.put("∂³(f(g))/∂y²∂t",      referenceRules.get("∂³(f(g))/∂x²∂t").replaceAll("x", "y"));
        referenceRules.put("∂³(f(g))/∂z²∂t",      referenceRules.get("∂³(f(g))/∂x²∂t").replaceAll("x", "z"));
        referenceRules.put("∂³(f(g))/∂x∂y∂z",     "∂³(f(g))/∂g³ * ∂g/∂x * ∂g/∂y * ∂g/∂z +" +
                                                  " ∂²(f(g))/∂g² * ∂g/∂y * ∂²g/∂x∂z +" +
                                                  " ∂²(f(g))/∂g² * ∂g/∂x * ∂²g/∂y∂z +" +
                                                  " ∂²(f(g))/∂g² * ∂²g/∂x∂y * ∂g/∂z +" +
                                                  " ∂(f(g))/∂g * ∂³g/∂x∂y∂z");
        referenceRules.put("∂³(f(g))/∂x∂y∂t",     referenceRules.get("∂³(f(g))/∂x∂y∂z").replaceAll("z", "t"));
        referenceRules.put("∂³(f(g))/∂x∂z∂t",     referenceRules.get("∂³(f(g))/∂x∂y∂t").replaceAll("y", "z"));
        referenceRules.put("∂³(f(g))/∂y∂z∂t",     referenceRules.get("∂³(f(g))/∂x∂z∂t").replaceAll("x", "y"));
        referenceRules.put("∂⁴(f(g))/∂x⁴",        "∂⁴(f(g))/∂g⁴ * ∂g/∂x * ∂g/∂x * ∂g/∂x * ∂g/∂x +" +
                                                  " 6 * ∂³(f(g))/∂g³ * ∂g/∂x * ∂g/∂x * ∂²g/∂x² +" +
                                                  " 3 * ∂²(f(g))/∂g² * ∂²g/∂x² * ∂²g/∂x² +" +
                                                  " 4 * ∂²(f(g))/∂g² * ∂g/∂x * ∂³g/∂x³ +" +
                                                  " ∂(f(g))/∂g * ∂⁴g/∂x⁴");
        referenceRules.put("∂⁴(f(g))/∂y⁴",        referenceRules.get("∂⁴(f(g))/∂x⁴").replaceAll("x", "y"));
        referenceRules.put("∂⁴(f(g))/∂z⁴",        referenceRules.get("∂⁴(f(g))/∂x⁴").replaceAll("x", "z"));
        referenceRules.put("∂⁴(f(g))/∂t⁴",        referenceRules.get("∂⁴(f(g))/∂x⁴").replaceAll("x", "t"));
        referenceRules.put("∂⁴(f(g))/∂x³∂y",      "∂⁴(f(g))/∂g⁴ * ∂g/∂x * ∂g/∂x * ∂g/∂x * ∂g/∂y +" +
                                                  " 3 * ∂³(f(g))/∂g³ * ∂g/∂x * ∂g/∂x * ∂²g/∂x∂y +" +
                                                  " 3 * ∂³(f(g))/∂g³ * ∂g/∂x * ∂²g/∂x² * ∂g/∂y +" +
                                                  " 3 * ∂²(f(g))/∂g² * ∂²g/∂x² * ∂²g/∂x∂y +" +
                                                  " 3 * ∂²(f(g))/∂g² * ∂g/∂x * ∂³g/∂x²∂y +" +
                                                  " ∂²(f(g))/∂g² * ∂³g/∂x³ * ∂g/∂y +" +
                                                  " ∂(f(g))/∂g * ∂⁴g/∂x³∂y");
        referenceRules.put("∂⁴(f(g))/∂x³∂z",      referenceRules.get("∂⁴(f(g))/∂x³∂y").replaceAll("y", "z"));
        referenceRules.put("∂⁴(f(g))/∂x³∂t",      referenceRules.get("∂⁴(f(g))/∂x³∂y").replaceAll("y", "t"));
        referenceRules.put("∂⁴(f(g))/∂x∂y³",      "∂⁴(f(g))/∂g⁴ * ∂g/∂x * ∂g/∂y * ∂g/∂y * ∂g/∂y +" +
                                                  " 3 * ∂³(f(g))/∂g³ * ∂g/∂y * ∂g/∂y * ∂²g/∂x∂y +" +
                                                  " 3 * ∂³(f(g))/∂g³ * ∂g/∂x * ∂g/∂y * ∂²g/∂y² +" +
                                                  " 3 * ∂²(f(g))/∂g² * ∂²g/∂x∂y * ∂²g/∂y² +" +
                                                  " 3 * ∂²(f(g))/∂g² * ∂g/∂y * ∂³g/∂x∂y² +" +
                                                  " ∂²(f(g))/∂g² * ∂g/∂x * ∂³g/∂y³ +" +
                                                  " ∂(f(g))/∂g * ∂⁴g/∂x∂y³");
        referenceRules.put("∂⁴(f(g))/∂x∂z³",      referenceRules.get("∂⁴(f(g))/∂x∂y³").replaceAll("y", "z"));
        referenceRules.put("∂⁴(f(g))/∂x∂t³",      referenceRules.get("∂⁴(f(g))/∂x∂y³").replaceAll("y", "t"));
        referenceRules.put("∂⁴(f(g))/∂y³∂z",      referenceRules.get("∂⁴(f(g))/∂x³∂z").replaceAll("x", "y"));
        referenceRules.put("∂⁴(f(g))/∂y³∂t",      referenceRules.get("∂⁴(f(g))/∂x³∂t").replaceAll("x", "y"));
        referenceRules.put("∂⁴(f(g))/∂y∂z³",      referenceRules.get("∂⁴(f(g))/∂x∂z³").replaceAll("x", "y"));
        referenceRules.put("∂⁴(f(g))/∂y∂t³",      referenceRules.get("∂⁴(f(g))/∂x∂t³").replaceAll("x", "y"));
        referenceRules.put("∂⁴(f(g))/∂z³∂t",      referenceRules.get("∂⁴(f(g))/∂x³∂t").replaceAll("x", "z"));
        referenceRules.put("∂⁴(f(g))/∂z∂t³",      referenceRules.get("∂⁴(f(g))/∂x∂t³").replaceAll("x", "z"));
        referenceRules.put("∂⁴(f(g))/∂x²∂y²",     "∂⁴(f(g))/∂g⁴ * ∂g/∂x * ∂g/∂x * ∂g/∂y * ∂g/∂y +" +
                                                  " 4 * ∂³(f(g))/∂g³ * ∂g/∂x * ∂g/∂y * ∂²g/∂x∂y +" +
                                                  " ∂³(f(g))/∂g³ * ∂g/∂x * ∂g/∂x * ∂²g/∂y² +" +
                                                  " 2 * ∂²(f(g))/∂g² * ∂²g/∂x∂y * ∂²g/∂x∂y +" +
                                                  " 2 * ∂²(f(g))/∂g² * ∂g/∂x * ∂³g/∂x∂y² +" +
                                                  " ∂³(f(g))/∂g³ * ∂²g/∂x² * ∂g/∂y * ∂g/∂y +" +
                                                  " 2 * ∂²(f(g))/∂g² * ∂g/∂y * ∂³g/∂x²∂y +" +
                                                  " ∂²(f(g))/∂g² * ∂²g/∂x² * ∂²g/∂y² +" +
                                                  " ∂(f(g))/∂g * ∂⁴g/∂x²∂y²");
        referenceRules.put("∂⁴(f(g))/∂x²∂z²",     referenceRules.get("∂⁴(f(g))/∂x²∂y²").replaceAll("y", "z"));
        referenceRules.put("∂⁴(f(g))/∂x²∂t²",     referenceRules.get("∂⁴(f(g))/∂x²∂y²").replaceAll("y", "t"));
        referenceRules.put("∂⁴(f(g))/∂y²∂z²",     referenceRules.get("∂⁴(f(g))/∂x²∂z²").replaceAll("x", "y"));
        referenceRules.put("∂⁴(f(g))/∂y²∂t²",     referenceRules.get("∂⁴(f(g))/∂x²∂t²").replaceAll("x", "y"));
        referenceRules.put("∂⁴(f(g))/∂z²∂t²",     referenceRules.get("∂⁴(f(g))/∂x²∂t²").replaceAll("x", "z"));

        referenceRules.put("∂⁴(f(g))/∂x²∂y∂z",    "∂⁴(f(g))/∂g⁴ * ∂g/∂x * ∂g/∂x * ∂g/∂y * ∂g/∂z +" +
                                                  " 2 * ∂³(f(g))/∂g³ * ∂g/∂x * ∂g/∂y * ∂²g/∂x∂z +" +
                                                  " ∂³(f(g))/∂g³ * ∂g/∂x * ∂g/∂x * ∂²g/∂y∂z +" +
                                                  " 2 * ∂³(f(g))/∂g³ * ∂g/∂x * ∂²g/∂x∂y * ∂g/∂z +" +
                                                  " 2 * ∂²(f(g))/∂g² * ∂²g/∂x∂y * ∂²g/∂x∂z +" +
                                                  " 2 * ∂²(f(g))/∂g² * ∂g/∂x * ∂³g/∂x∂y∂z +" +
                                                  " ∂³(f(g))/∂g³ * ∂²g/∂x² * ∂g/∂y * ∂g/∂z +" +
                                                  " ∂²(f(g))/∂g² * ∂g/∂y * ∂³g/∂x²∂z +" +
                                                  " ∂²(f(g))/∂g² * ∂²g/∂x² * ∂²g/∂y∂z +" +
                                                  " ∂²(f(g))/∂g² * ∂³g/∂x²∂y * ∂g/∂z +" +
                                                  " ∂(f(g))/∂g * ∂⁴g/∂x²∂y∂z");
        referenceRules.put("∂⁴(f(g))/∂x²∂y∂t",    referenceRules.get("∂⁴(f(g))/∂x²∂y∂z").replaceAll("z", "t"));
        referenceRules.put("∂⁴(f(g))/∂x²∂z∂t",    referenceRules.get("∂⁴(f(g))/∂x²∂y∂t").replaceAll("y", "z"));
        referenceRules.put("∂⁴(f(g))/∂x∂y²∂z",    "∂⁴(f(g))/∂g⁴ * ∂g/∂x * ∂g/∂y * ∂g/∂y * ∂g/∂z +" +
                                                  " ∂³(f(g))/∂g³ * ∂g/∂y * ∂g/∂y * ∂²g/∂x∂z +" +
                                                  " 2 * ∂³(f(g))/∂g³ * ∂g/∂x * ∂g/∂y * ∂²g/∂y∂z +" +
                                                  " 2 * ∂³(f(g))/∂g³ * ∂g/∂y * ∂²g/∂x∂y * ∂g/∂z +" +
                                                  " 2 * ∂²(f(g))/∂g² * ∂²g/∂x∂y * ∂²g/∂y∂z +" +
                                                  " 2 * ∂²(f(g))/∂g² * ∂g/∂y * ∂³g/∂x∂y∂z +" +
                                                  " ∂³(f(g))/∂g³ * ∂g/∂x * ∂²g/∂y² * ∂g/∂z +" +
                                                  " ∂²(f(g))/∂g² * ∂²g/∂y² * ∂²g/∂x∂z +" +
                                                  " ∂²(f(g))/∂g² * ∂g/∂x * ∂³g/∂y²∂z +" +
                                                  " ∂²(f(g))/∂g² * ∂³g/∂x∂y² * ∂g/∂z +" +
                                                  " ∂(f(g))/∂g * ∂⁴g/∂x∂y²∂z");
        referenceRules.put("∂⁴(f(g))/∂x∂y²∂t",    referenceRules.get("∂⁴(f(g))/∂x∂y²∂z").replaceAll("z", "t"));
        referenceRules.put("∂⁴(f(g))/∂y²∂z∂t",    referenceRules.get("∂⁴(f(g))/∂x²∂z∂t").replaceAll("x", "y"));
        referenceRules.put("∂⁴(f(g))/∂x∂y∂z²",    "∂⁴(f(g))/∂g⁴ * ∂g/∂x * ∂g/∂y * ∂g/∂z * ∂g/∂z +" +
                                                  " 2 * ∂³(f(g))/∂g³ * ∂g/∂y * ∂g/∂z * ∂²g/∂x∂z +" +
                                                  " 2 * ∂³(f(g))/∂g³ * ∂g/∂x * ∂g/∂z * ∂²g/∂y∂z +" +
                                                  " ∂³(f(g))/∂g³ * ∂g/∂x * ∂g/∂y * ∂²g/∂z² +" +
                                                  " 2 * ∂²(f(g))/∂g² * ∂²g/∂x∂z * ∂²g/∂y∂z +" +
                                                  " ∂²(f(g))/∂g² * ∂g/∂y * ∂³g/∂x∂z² +" +
                                                  " ∂²(f(g))/∂g² * ∂g/∂x * ∂³g/∂y∂z² +" +
                                                  " ∂³(f(g))/∂g³ * ∂²g/∂x∂y * ∂g/∂z * ∂g/∂z +" +
                                                  " 2 * ∂²(f(g))/∂g² * ∂g/∂z * ∂³g/∂x∂y∂z +" +
                                                  " ∂²(f(g))/∂g² * ∂²g/∂x∂y * ∂²g/∂z² +" +
                                                  " ∂(f(g))/∂g * ∂⁴g/∂x∂y∂z²");
        referenceRules.put("∂⁴(f(g))/∂x∂z²∂t",    referenceRules.get("∂⁴(f(g))/∂x∂y²∂t").replaceAll("y", "z"));
        referenceRules.put("∂⁴(f(g))/∂y∂z²∂t",    referenceRules.get("∂⁴(f(g))/∂x∂z²∂t").replaceAll("x", "y"));
        referenceRules.put("∂⁴(f(g))/∂x∂y∂t²",    referenceRules.get("∂⁴(f(g))/∂x∂y∂z²").replaceAll("z", "t"));
        referenceRules.put("∂⁴(f(g))/∂x∂z∂t²",    referenceRules.get("∂⁴(f(g))/∂x∂y∂t²").replaceAll("y", "z"));
        referenceRules.put("∂⁴(f(g))/∂y∂z∂t²",    referenceRules.get("∂⁴(f(g))/∂x∂z∂t²").replaceAll("x", "y"));
        referenceRules.put("∂⁴(f(g))/∂x∂y∂z∂t",   "∂⁴(f(g))/∂g⁴ * ∂g/∂x * ∂g/∂y * ∂g/∂z * ∂g/∂t +" +
                                                  " ∂³(f(g))/∂g³ * ∂g/∂y * ∂g/∂z * ∂²g/∂x∂t +" +
                                                  " ∂³(f(g))/∂g³ * ∂g/∂x * ∂g/∂z * ∂²g/∂y∂t +" +
                                                  " ∂³(f(g))/∂g³ * ∂g/∂x * ∂g/∂y * ∂²g/∂z∂t +" +
                                                  " ∂³(f(g))/∂g³ * ∂g/∂y * ∂²g/∂x∂z * ∂g/∂t +" +
                                                  " ∂²(f(g))/∂g² * ∂²g/∂x∂z * ∂²g/∂y∂t +" +
                                                  " ∂²(f(g))/∂g² * ∂g/∂y * ∂³g/∂x∂z∂t +" +
                                                  " ∂³(f(g))/∂g³ * ∂g/∂x * ∂²g/∂y∂z * ∂g/∂t +" +
                                                  " ∂²(f(g))/∂g² * ∂²g/∂y∂z * ∂²g/∂x∂t +" +
                                                  " ∂²(f(g))/∂g² * ∂g/∂x * ∂³g/∂y∂z∂t +" +
                                                  " ∂³(f(g))/∂g³ * ∂²g/∂x∂y * ∂g/∂z * ∂g/∂t +" +
                                                  " ∂²(f(g))/∂g² * ∂g/∂z * ∂³g/∂x∂y∂t +" +
                                                  " ∂²(f(g))/∂g² * ∂²g/∂x∂y * ∂²g/∂z∂t +" +
                                                  " ∂²(f(g))/∂g² * ∂³g/∂x∂y∂z * ∂g/∂t +" +
                                                  " ∂(f(g))/∂g * ∂⁴g/∂x∂y∂z∂t");

        Field compFieldArrayField = DSCompiler.class.getDeclaredField("compIndirection");
        compFieldArrayField.setAccessible(true);
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
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                DSCompiler compiler = DSCompiler.getCompiler(i, j);
                Object[][] compIndirection = (Object[][]) compFieldArrayField.get(compiler);
                for (int k = 0; k < compIndirection.length; ++k) {
                    String product = ordersToString(compiler.getPartialDerivativeOrders(k),
                                                    "(f(g))", "x", "y", "z", "t");
                    StringBuilder rule = new StringBuilder();
                    for (Object term : compIndirection[k]) {
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
                                                       "g", "x", "y", "z", "t"));
                        }
                    }
                    Assert.assertEquals(product, referenceRules.get(product), rule.toString());
                }
            }
        }
    }

    @Test
    public void testRebaserRules()
        throws IllegalAccessException, IllegalArgumentException, InvocationTargetException,
               NoSuchMethodException, SecurityException, NoSuchFieldException {

        Method getterMethod = DSCompiler.class.getDeclaredMethod("getRebaser", DSCompiler.class);
        getterMethod.setAccessible(true);

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

        // assuming f = f(p₀, p₁)
        //          p₀ = p₀(q₀, q₁, q₂)
        //          p₁ = p₁(q₀, q₁, q₂)

        for (int order = 0; order < 4; ++order) {

            DSCompiler c2 = DSCompiler.getCompiler(2, order);
            DSCompiler c3 = DSCompiler.getCompiler(3, order);
            int baseSize = c3.getSize();
            Object[][] rebaser = (Object[][]) getterMethod.invoke(c2, c3);

            Assert.assertEquals(c3.getSize(), rebaser.length);

            // composition rule for function value
            Object[] fRule = rebaser[c3.getPartialDerivativeIndex(0, 0, 0)];
            Assert.assertEquals(1, fRule.length);
            Assert.assertEquals(1, ((Integer) coeffMethod.invoke(fRule[0])).intValue());
            Assert.assertEquals(0, ((Integer) dsIndexField.get(fRule[0])).intValue());
            Assert.assertEquals(0, ((int[]) productIndicesField.get(fRule[0])).length);

            if (order > 0) {
                // composition rules for first derivatives

                // ∂f/∂q₀        = ∂f/∂p₀ ∂p₀/∂q₀ + ∂f/∂p₁ ∂p₁/∂q₀
                Object[] dFdQ0Rule = rebaser[c3.getPartialDerivativeIndex(1, 0, 0)];
                Assert.assertEquals(2, dFdQ0Rule.length);
                Assert.assertEquals(1, ((Integer) coeffMethod.invoke(dFdQ0Rule[0])).intValue());
                Assert.assertEquals(c2.getPartialDerivativeIndex(1, 0), ((Integer) dsIndexField.get(dFdQ0Rule[0])).intValue());
                Assert.assertEquals(1, ((int[]) productIndicesField.get(dFdQ0Rule[0])).length);
                Assert.assertEquals(c3.getPartialDerivativeIndex(1, 0, 0), ((int[]) productIndicesField.get(dFdQ0Rule[0]))[0]);
                Assert.assertEquals(1, ((Integer) coeffMethod.invoke(dFdQ0Rule[1])).intValue());
                Assert.assertEquals(c2.getPartialDerivativeIndex(0, 1), ((Integer) dsIndexField.get(dFdQ0Rule[1])).intValue());
                Assert.assertEquals(1, ((int[]) productIndicesField.get(dFdQ0Rule[1])).length);
                Assert.assertEquals(baseSize + c3.getPartialDerivativeIndex(1, 0, 0), ((int[]) productIndicesField.get(dFdQ0Rule[1]))[0]);

                // ∂f/∂q₁        = ∂f/∂p₀ ∂p₀/∂q₁ + ∂f/∂p₁ ∂p₁/∂q₁
                Object[] dFdQ1Rule = rebaser[c3.getPartialDerivativeIndex(0, 1, 0)];
                Assert.assertEquals(2, dFdQ1Rule.length);
                Assert.assertEquals(1, ((Integer) coeffMethod.invoke(dFdQ1Rule[0])).intValue());
                Assert.assertEquals(c2.getPartialDerivativeIndex(1, 0), ((Integer) dsIndexField.get(dFdQ1Rule[0])).intValue());
                Assert.assertEquals(1, ((int[]) productIndicesField.get(dFdQ1Rule[0])).length);
                Assert.assertEquals(c3.getPartialDerivativeIndex(0, 1, 0), ((int[]) productIndicesField.get(dFdQ1Rule[0]))[0]);
                Assert.assertEquals(1, ((Integer) coeffMethod.invoke(dFdQ1Rule[1])).intValue());
                Assert.assertEquals(c2.getPartialDerivativeIndex(0, 1), ((Integer) dsIndexField.get(dFdQ1Rule[1])).intValue());
                Assert.assertEquals(1, ((int[]) productIndicesField.get(dFdQ1Rule[1])).length);
                Assert.assertEquals(baseSize + c3.getPartialDerivativeIndex(0, 1, 0), ((int[]) productIndicesField.get(dFdQ1Rule[1]))[0]);

                // ∂f/∂q₂        = ∂f/∂p₀ ∂p₀/∂q₂ + ∂f/∂p₁ ∂p₁/∂q₂
                Object[] dFdQ2Rule = rebaser[c3.getPartialDerivativeIndex(0, 0, 1)];
                Assert.assertEquals(2, dFdQ2Rule.length);
                Assert.assertEquals(1, ((Integer) coeffMethod.invoke(dFdQ2Rule[0])).intValue());
                Assert.assertEquals(c2.getPartialDerivativeIndex(1, 0), ((Integer) dsIndexField.get(dFdQ2Rule[0])).intValue());
                Assert.assertEquals(1, ((int[]) productIndicesField.get(dFdQ2Rule[0])).length);
                Assert.assertEquals(c3.getPartialDerivativeIndex(0, 0, 1), ((int[]) productIndicesField.get(dFdQ2Rule[0]))[0]);
                Assert.assertEquals(1, ((Integer) coeffMethod.invoke(dFdQ2Rule[1])).intValue());
                Assert.assertEquals(c2.getPartialDerivativeIndex(0, 1), ((Integer) dsIndexField.get(dFdQ2Rule[1])).intValue());
                Assert.assertEquals(1, ((int[]) productIndicesField.get(dFdQ2Rule[1])).length);
                Assert.assertEquals(baseSize + c3.getPartialDerivativeIndex(0, 0, 1), ((int[]) productIndicesField.get(dFdQ2Rule[1]))[0]);

            }

            if (order > 1) {
                // composition rules for second derivatives

                // ∂²f/∂q₀²      = ∂²f/∂p₀² (∂p₀/∂q₀)² + 2 ∂²f/∂p₀∂p₁ ∂p₀/∂q₀ ∂p₁/∂q₀ + ∂²f/∂p₁² (∂p₁/∂q₀)²
                //               + ∂f/∂p₀ ∂²p₀/∂q₀² + ∂f/∂p₁ ∂²p₁/∂q₀²
                Object[] d2FdQ02Rule = rebaser[c3.getPartialDerivativeIndex(2, 0, 0)];
                Assert.assertEquals(5, d2FdQ02Rule.length);
                Assert.assertEquals(1, ((Integer) coeffMethod.invoke(d2FdQ02Rule[0])).intValue());
                Assert.assertEquals(c2.getPartialDerivativeIndex(2, 0), ((Integer) dsIndexField.get(d2FdQ02Rule[0])).intValue());
                Assert.assertEquals(2, ((int[]) productIndicesField.get(d2FdQ02Rule[0])).length);
                Assert.assertEquals(c3.getPartialDerivativeIndex(1, 0, 0), ((int[]) productIndicesField.get(d2FdQ02Rule[0]))[0]);
                Assert.assertEquals(c3.getPartialDerivativeIndex(1, 0, 0), ((int[]) productIndicesField.get(d2FdQ02Rule[0]))[1]);
                Assert.assertEquals(2, ((Integer) coeffMethod.invoke(d2FdQ02Rule[1])).intValue());
                Assert.assertEquals(c2.getPartialDerivativeIndex(1, 1), ((Integer) dsIndexField.get(d2FdQ02Rule[1])).intValue());
                Assert.assertEquals(2, ((int[]) productIndicesField.get(d2FdQ02Rule[1])).length);
                Assert.assertEquals(c3.getPartialDerivativeIndex(1, 0, 0), ((int[]) productIndicesField.get(d2FdQ02Rule[1]))[0]);
                Assert.assertEquals(baseSize + c3.getPartialDerivativeIndex(1, 0, 0), ((int[]) productIndicesField.get(d2FdQ02Rule[1]))[1]);
                Assert.assertEquals(1, ((Integer) coeffMethod.invoke(d2FdQ02Rule[2])).intValue());
                Assert.assertEquals(c2.getPartialDerivativeIndex(0, 2), ((Integer) dsIndexField.get(d2FdQ02Rule[2])).intValue());
                Assert.assertEquals(2, ((int[]) productIndicesField.get(d2FdQ02Rule[2])).length);
                Assert.assertEquals(baseSize + c3.getPartialDerivativeIndex(1, 0, 0), ((int[]) productIndicesField.get(d2FdQ02Rule[2]))[0]);
                Assert.assertEquals(baseSize + c3.getPartialDerivativeIndex(1, 0, 0), ((int[]) productIndicesField.get(d2FdQ02Rule[2]))[1]);
                Assert.assertEquals(1, ((Integer) coeffMethod.invoke(d2FdQ02Rule[3])).intValue());
                Assert.assertEquals(c2.getPartialDerivativeIndex(1, 0), ((Integer) dsIndexField.get(d2FdQ02Rule[3])).intValue());
                Assert.assertEquals(1, ((int[]) productIndicesField.get(d2FdQ02Rule[3])).length);
                Assert.assertEquals(c3.getPartialDerivativeIndex(2, 0, 0), ((int[]) productIndicesField.get(d2FdQ02Rule[3]))[0]);
                Assert.assertEquals(1, ((Integer) coeffMethod.invoke(d2FdQ02Rule[4])).intValue());
                Assert.assertEquals(c2.getPartialDerivativeIndex(0, 1), ((Integer) dsIndexField.get(d2FdQ02Rule[4])).intValue());
                Assert.assertEquals(1, ((int[]) productIndicesField.get(d2FdQ02Rule[4])).length);
                Assert.assertEquals(baseSize + c3.getPartialDerivativeIndex(2, 0, 0), ((int[]) productIndicesField.get(d2FdQ02Rule[4]))[0]);

                // ∂²f/∂q₁²      = ∂²f/∂p₀² (∂p₀/∂q₁)² + 2 ∂²f/∂p₀∂p₁ ∂p₀/∂q₁ ∂p₁/∂q₁ + ∂²f/∂p₁² (∂p₁/∂q₁)²
                //               + ∂f/∂p₀ ∂²p₀/∂q₁² + ∂f/∂p₁ ∂²p₁/∂q₁²
                // ∂²f/∂q₂²      = ∂²f/∂p₀² (∂p₀/∂q₂)² + 2 ∂²f/∂p₀∂p₁ ∂p₀/∂q₂ ∂p₁/∂q₂ + ∂²f/∂p₁² (∂p₁/∂q₂)²
                //               + ∂f/∂p₀ ∂²p₀/∂q₂² + ∂f/∂p₁ ∂²p₁/∂q₂²
                // ∂²f/∂q₀∂q₁    = ∂²f/∂p₀² ∂p₀/∂q₀ ∂p₀/∂q₁ + ∂²f/∂p₀∂p₁ (∂p₁/∂q₀ ∂p₀/∂q₁ + ∂p₀/∂q₀ ∂p₁/∂q₁) + ∂²f/∂p₁² ∂p₁/∂q₀ ∂p₁/∂q₁
                //               + ∂f/∂p₀ ∂²p₀/∂q₀∂q₁ + ∂f/∂p₁ ∂²p₁/∂q₀∂q₁
                // ∂²f/∂q₀∂q₂    = ∂²f/∂p₀² ∂p₀/∂q₀ ∂p₀/∂q₂ + ∂²f/∂p₀∂p₁ (∂p₁/∂q₀ ∂p₀/∂q₂ + ∂p₀/∂q₀ ∂p₁/∂q₂) + ∂²f/∂p₁² ∂p₁/∂q₀ ∂p₁/∂q₂
                //               + ∂f/∂p₀ ∂²p₀/∂q₀∂q₂ + ∂f/∂p₁ ∂²p₁/∂q₀∂q₂
                // ∂²f/∂q₁∂q₂    = ∂²f/∂p₀² ∂p₀/∂q₁ ∂p₀/∂q₂ + ∂²f/∂p₀∂p₁ (∂p₁/∂q₁ ∂p₀/∂q₂ + ∂p₀/∂q₁ ∂p₁/∂q₂) + ∂²f/∂p₁² ∂p₁/∂q₁ ∂p₁/∂q₂
                //               + ∂f/∂p₀ ∂²p₀/∂q₁∂q₂ + ∂f/∂p₁ ∂²p₁/∂q₁∂q₂

            }

            if (order > 2) {
                // composition rules for third derivatives

                // ∂³f/∂q₀³      = ∂³f/∂p₀³    (∂p₀/∂q₀)³ + 3 ∂³f/∂p₀²∂p₁ (∂p₀/∂q₀)² ∂p₁/∂q₀ + 3 ∂³f/∂p₀∂p₁² ∂p₀/∂q₀ (∂p₁/∂q₀)² + ∂³f/∂p₁³ (∂p₁/∂q₀)³
                //               + …
                // ∂³f/∂q₁³      = ∂³f/∂p₀³    (∂p₀/∂q₁)³ + 3 ∂³f/∂p₀²∂p₁ (∂p₀/∂q₁)² ∂p₁/∂q₁ + 3 ∂³f/∂p₀∂p₁² ∂p₀/∂q₁ (∂p₁/∂q₁)² + ∂³f/∂p₁³ (∂p₁/∂q₁)³
                //               + …
                // ∂³f/∂q₂³      = ∂³f/∂p₀³    (∂p₀/∂q₂)³ + 3 ∂³f/∂p₀²∂p₁ (∂p₀/∂q₂)² ∂p₁/∂q₂ + 3 ∂³f/∂p₀∂p₁² ∂p₀/∂q₂ (∂p₁/∂q₂)² + ∂³f/∂p₁³ (∂p₁/∂q₂)³
                //               + …
                // ∂³f/∂q₀²∂q₁   = ∂³f/∂p₀³    (∂p₀/∂q₀)² ∂p₀/∂q₁ +
                //                 ∂³f/∂p₀²∂p₁ ((∂p₀/∂q₀)² ∂p₁/∂q₁ + ∂p₀/∂q₀ ∂p₀/∂q₁ ∂p₁/∂q₀ + ∂p₀/∂q₀ ∂p₀/∂q₁ ∂p₁/∂q₀) +
                //                 ∂³f/∂p₀∂p₁² (∂p₀/∂q₀ ∂p₁/∂q₀ ∂p₁/∂q₁ + ∂p₀/∂q₁ (∂p₁/∂q₀)² + ∂p₀/∂q₀ ∂p₁/∂q₀ ∂p₁/∂q₁) +
                //                 ∂³f/∂p₁³    (∂p₁/∂q₀)² ∂p₁/∂q₁
                //               + …
                // ∂³f/∂q₀²∂q₂   = ∂³f/∂p₀³    (∂p₀/∂q₀)² ∂p₀/∂q₂ +
                //                 ∂³f/∂p₀²∂p₁ ((∂p₀/∂q₀)² ∂p₁/∂q₂ + ∂p₀/∂q₀ ∂p₀/∂q₂ ∂p₁/∂q₀ + ∂p₀/∂q₀ ∂p₀/∂q₂ ∂p₁/∂q₀) +
                //                 ∂³f/∂p₀∂p₁² (∂p₀/∂q₀ ∂p₁/∂q₀ ∂p₁/∂q₂ + ∂p₀/∂q₂ (∂p₁/∂q₀)² + ∂p₀/∂q₀ ∂p₁/∂q₀ ∂p₁/∂q₂) +
                //                 ∂³f/∂p₁³    (∂p₁/∂q₀)² ∂p₁/∂q₂
                //               + …
                // ∂³f/∂q₁²∂q₂   = ∂³f/∂p₀³    (∂p₀/∂q₁)² ∂p₀/∂q₂ +
                //                 ∂³f/∂p₀²∂p₁ ((∂p₀/∂q₁)² ∂p₁/∂q₂ + ∂p₀/∂q₁ ∂p₀/∂q₂ ∂p₁/∂q₁ + ∂p₀/∂q₁ ∂p₀/∂q₂ ∂p₁/∂q₁) +
                //                 ∂³f/∂p₀∂p₁² (∂p₀/∂q₁ ∂p₁/∂q₁ ∂p₁/∂q₂ + ∂p₀/∂q₂ (∂p₁/∂q₁)² + ∂p₀/∂q₁ ∂p₁/∂q₁ ∂p₁/∂q₂) +
                //                 ∂³f/∂p₁³    (∂p₁/∂q₁)² ∂p₁/∂q₂
                //               + …
                // ∂³f/∂q₀∂q₁∂q₂ = ∂³f/∂p₀³    ∂p₀/∂q₀ ∂p₀/∂q₁ ∂p₀/∂q₂ +
                //                 ∂³f/∂p₀²∂p₁ (∂p₀/∂q₂ ∂p₀/∂q₀ ∂p₁/∂q₁ + ∂p₀/∂q₀ ∂p₀/∂q₁ ∂p₁/∂q₂ + ∂p₀/∂q₁ ∂p₀/∂q₂ ∂p₁/∂q₀) +
                //                 ∂³f/∂p₀∂p₁² (∂p₀/∂q₂ ∂p₁/∂q₀ ∂p₁/∂q₁ + ∂p₀/∂q₁ ∂p₁/∂q₀ ∂p₁/∂q₂ + ∂p₀/∂q₀ ∂p₁/∂q₁ ∂p₁/∂q₂) +
                //                 ∂³f/∂p₁³    ∂p₁/∂q₀ ∂p₁/∂q₁ ∂p₁/∂q₂
                //               + …

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
