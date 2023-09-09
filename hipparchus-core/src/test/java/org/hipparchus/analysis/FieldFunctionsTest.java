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
package org.hipparchus.analysis;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.util.Binary64;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.MathArrays;
import org.junit.Assert;
import org.junit.Test;

public class FieldFunctionsTest {

    @Test
    public void testScalarUnivariateFunctionConversion() {
        FieldUnivariateFunction f1 = new FieldUnivariateFunction() {
            public <T extends CalculusFieldElement<T>> T value(T x) {
                return x.multiply(2);
            }
        };
        CalculusFieldUnivariateFunction<Binary64> f1Converted = f1.toCalculusFieldUnivariateFunction(Binary64Field.getInstance());
        CalculusFieldUnivariateFunction<Binary64> f2 = x -> x.multiply(2);

        for (double x = 0; x < 1; x += 0.01) {
            Assert.assertEquals(f2.value(new Binary64(x)).getReal(),
                                f1Converted.value(new Binary64(x)).getReal(),
                                1.0e-15);
        }
    }

    @Test
    public void testScalarMultivariateFunctionConversion() {
        FieldMultivariateFunction f1 = new FieldMultivariateFunction() {
            public <T extends CalculusFieldElement<T>> T value(@SuppressWarnings("unchecked") T... x) {
                return x[0].multiply(2).add(x[1]);
            }
        };
        CalculusFieldMultivariateFunction<Binary64> f1Converted = f1.toCalculusFieldMultivariateFunction(Binary64Field.getInstance());
        CalculusFieldMultivariateFunction<Binary64> f2 = x -> x[0].multiply(2).add(x[1]);

        for (double x0 = 0; x0 < 1; x0 += 0.01) {
            for (double x1 = 0; x1 < 1; x1 += 0.01) {
                Assert.assertEquals(f2.value(new Binary64(x0), new Binary64(x1)).getReal(),
                                    f1Converted.value(new Binary64(x0), new Binary64(x1)).getReal(),
                                    1.0e-15);
            }
        }
    }

    @Test
    public void testVectorUnivariateFunctionConversion() {
        FieldUnivariateVectorFunction f1 = new FieldUnivariateVectorFunction() {
            public <T extends CalculusFieldElement<T>> T[] value(T x) {
                T[] y = MathArrays.buildArray(x.getField(), 3);
                y[0] = x.add(1);
                y[1] = x.multiply(2);
                y[2] = x.multiply(x);
                return y;
            }
        };
        CalculusFieldUnivariateVectorFunction<Binary64> f1Converted = f1.toCalculusFieldUnivariateVectorFunction(Binary64Field.getInstance());
        CalculusFieldUnivariateVectorFunction<Binary64> f2 = x -> new Binary64[] {
            x.add(1), x.multiply(2), x.multiply(x)
        };

        for (double x = 0; x < 1; x += 0.01) {
            for (int i = 0; i < 3; ++i) {
                Assert.assertEquals(f2.value(new Binary64(x))[i].getReal(),
                                    f1Converted.value(new Binary64(x))[i].getReal(),
                                    1.0e-15);
            }
        }
    }

    @Test
    public void testVectorMultivariateFunctionConversion() {
        FieldMultivariateVectorFunction f1 = new FieldMultivariateVectorFunction() {
            public <T extends CalculusFieldElement<T>> T[] value(@SuppressWarnings("unchecked") T... x) {
                T[] y = MathArrays.buildArray(x[0].getField(), 3);
                y[0] = x[0].add(1);
                y[1] = x[1].multiply(2);
                y[2] = x[0].multiply(x[1]);
                return y;
            }
        };
        CalculusFieldMultivariateVectorFunction<Binary64> f1Converted = f1.toCalculusFieldMultivariateVectorFunction(Binary64Field.getInstance());
        CalculusFieldMultivariateVectorFunction<Binary64> f2 = x -> new Binary64[] {
            x[0].add(1), x[1].multiply(2), x[0].multiply(x[1])
        };

        for (double x0 = 0; x0 < 1; x0 += 0.01) {
            for (double x1 = 0; x1 < 1; x1 += 0.01) {
                for (int i = 0; i < 3; ++i) {
                    Assert.assertEquals(f2.value(new Binary64(x0), new Binary64(x1))[i].getReal(),
                                        f1Converted.value(new Binary64(x0), new Binary64(x1))[i].getReal(),
                                        1.0e-15);
                }
            }
        }
    }

    @Test
    public void testMatrixUnivariateFunctionConversion() {
        FieldUnivariateMatrixFunction f1 = new FieldUnivariateMatrixFunction() {
            public <T extends CalculusFieldElement<T>> T[][] value(T x) {
                T[][] y = MathArrays.buildArray(x.getField(), 2, 2);
                y[0][0] = x.add(1);
                y[0][1] = x.multiply(2);
                y[1][0] = x.multiply(x);
                y[1][1] = x.sin();
                return y;
            }
        };
        CalculusFieldUnivariateMatrixFunction<Binary64> f1Converted = f1.toCalculusFieldUnivariateMatrixFunction(Binary64Field.getInstance());
        CalculusFieldUnivariateMatrixFunction<Binary64> f2 = x -> new Binary64[][] {
            { x.add(1), x.multiply(2) },
            { x.multiply(x), x.sin() }
        };

        for (double x = 0; x < 1; x += 0.01) {
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    Assert.assertEquals(f2.value(new Binary64(x))[i][j].getReal(),
                                        f1Converted.value(new Binary64(x))[i][j].getReal(),
                                        1.0e-15);
                }
            }
        }
    }

    @Test
    public void testMatrixMultivariateFunctionConversion() {
        FieldMultivariateMatrixFunction f1 = new FieldMultivariateMatrixFunction() {
            public <T extends CalculusFieldElement<T>> T[][] value(@SuppressWarnings("unchecked") T... x) {
                T[][] y = MathArrays.buildArray(x[0].getField(), 2, 2);
                y[0][0] = x[0].add(1);
                y[0][1] = x[1].multiply(2);
                y[1][0] = x[0].multiply(x[1]);
                y[1][1] = x[1].sin();
                return y;
            }
        };
        CalculusFieldMultivariateMatrixFunction<Binary64> f1Converted = f1.toCalculusFieldMultivariateMatrixFunction(Binary64Field.getInstance());
        CalculusFieldMultivariateMatrixFunction<Binary64> f2 = x -> new Binary64[][] {
            { x[0].add(1), x[1].multiply(2) },
            { x[0].multiply(x[1]), x[1].sin() }
        };

        for (double x0 = 0; x0 < 1; x0 += 0.01) {
            for (double x1 = 0; x1 < 1; x1 += 0.01) {
                for (int i = 0; i < 2; ++i) {
                    for (int j = 0; j < 2; ++j) {
                        Assert.assertEquals(f2.value(new Binary64(x0), new Binary64(x1))[i][j].getReal(),
                                            f1Converted.value(new Binary64(x0), new Binary64(x1))[i][j].getReal(),
                                            1.0e-15);
                    }
                }
            }
        }
    }

}
