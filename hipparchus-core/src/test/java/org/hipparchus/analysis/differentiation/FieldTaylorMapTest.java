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

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.linear.FieldQRDecomposer;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;
import org.junit.jupiter.api.Test;

import java.lang.reflect.Array;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test for class {@linkFieldTaylorMap<T>}.
 */
class FieldTaylorMapTest {

    @Test
    void testNullPoint() {
        doTestNullPoint(Binary64Field.getInstance());
    }

    @Test
    void testDim0Point() {
        doTestDim0Point(Binary64Field.getInstance());
    }

    @Test
    void testNullFunctions() {
        doTestNullFunctions(Binary64Field.getInstance());
    }

    @Test
    void testNoFunctions() {
        doTestNoFunctions(Binary64Field.getInstance());
    }

    @Test
    void testIncompatiblePointAndFunctions() {
        doTestIncompatiblePointAndFunctions(Binary64Field.getInstance());
    }

    @Test
    void testIncompatible() {
        doTestIncompatible(Binary64Field.getInstance());
    }

    @Test
    void testNbParameters() {
        doTestNbParameters(Binary64Field.getInstance());
    }

    @Test
    void testNbFunctions() {
        doTestNbFunctions(Binary64Field.getInstance());
    }

    @Test
    void testIdentity() {
        doTestIdentity(Binary64Field.getInstance());
    }

    @Test
    void testValue() {
        doTestValue(Binary64Field.getInstance());
    }

    @Test
    void testCompose() {
        doTestCompose(Binary64Field.getInstance());
    }

    @Test
    void testInvertNonSquare() {
        doTestInvertNonSquare(Binary64Field.getInstance());
    }

    @Test
    void testInvertMonoDimensional() {
        doTestInvertMonoDimensional(Binary64Field.getInstance());
    }

    @Test
    void testInvertBiDimensional() {
        doTestInvertBiDimensional(Binary64Field.getInstance());
    }

    @SuppressWarnings("unchecked")
    private <T extends CalculusFieldElement<T>> void doTestNullPoint(final Field<T> field) {
        try {
            new FieldTaylorMap<>(null, (FieldDerivativeStructure<T>[]) Array.newInstance(FieldDerivativeStructure.class, 2));
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.NUMBER_OF_ELEMENTS_SHOULD_BE_POSITIVE, miae.getSpecifier());
            assertEquals(0, ((Integer) miae.getParts()[0]).intValue());
        }
    }

    @SuppressWarnings("unchecked")
    private <T extends CalculusFieldElement<T>> void doTestDim0Point(final Field<T> field) {
        try {
            new FieldTaylorMap<>(MathArrays.buildArray(field, 0), (FieldDerivativeStructure<T>[]) Array.newInstance(FieldDerivativeStructure.class, 2));
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.NUMBER_OF_ELEMENTS_SHOULD_BE_POSITIVE, miae.getSpecifier());
            assertEquals(0, ((Integer) miae.getParts()[0]).intValue());
        }
    }

    private <T extends CalculusFieldElement<T>> void doTestNullFunctions(final Field<T> field) {
        try {
            new FieldTaylorMap<>(MathArrays.buildArray(field, 2), null);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.NUMBER_OF_ELEMENTS_SHOULD_BE_POSITIVE, miae.getSpecifier());
            assertEquals(0, ((Integer) miae.getParts()[0]).intValue());
        }
    }

    @SuppressWarnings("unchecked")
    private <T extends CalculusFieldElement<T>> void doTestNoFunctions(final Field<T> field) {
        try {
            new FieldTaylorMap<>(MathArrays.buildArray(field, 2), (FieldDerivativeStructure<T>[]) Array.newInstance(FieldDerivativeStructure.class, 0));
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.NUMBER_OF_ELEMENTS_SHOULD_BE_POSITIVE, miae.getSpecifier());
            assertEquals(0, ((Integer) miae.getParts()[0]).intValue());
        }
    }

    private <T extends CalculusFieldElement<T>> void doTestIncompatiblePointAndFunctions(final Field<T> field) {
        FDSFactory<T> factory = new FDSFactory<>(field, 6, 6);
        @SuppressWarnings("unchecked")
        FieldDerivativeStructure<T>[] functions = (FieldDerivativeStructure<T>[]) Array.newInstance(FieldDerivativeStructure.class,
                                                                                                    factory.getCompiler().getFreeParameters());
        for (int i = 0; i < functions.length; ++i) {
            functions[i] = factory.constant(0);
        }
        try {
            new FieldTaylorMap<>(MathArrays.buildArray(field, functions.length - 1), functions);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
            assertEquals(5, ((Integer) miae.getParts()[0]).intValue());
            assertEquals(6, ((Integer) miae.getParts()[1]).intValue());
        }
    }

    private <T extends CalculusFieldElement<T>> void doTestIncompatible(final Field<T> field) {
        FDSFactory<T> factory = new FDSFactory<>(field, 6, 6);
        @SuppressWarnings("unchecked")
        FieldDerivativeStructure<T>[] functions = (FieldDerivativeStructure<T>[]) Array.newInstance(FieldDerivativeStructure.class,
                                                                                                    factory.getCompiler().getFreeParameters());
        for (int i = 0; i < functions.length - 1; ++i) {
            functions[i] = factory.constant(0);
        }
        functions[functions.length - 1] = new FDSFactory<>(field,
                                                           factory.getCompiler().getFreeParameters(),
                                                           factory.getCompiler().getOrder() - 1).
                                          constant(1.0);
        try {
            new FieldTaylorMap<>(MathArrays.buildArray(field, functions.length), functions);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
            assertEquals(6, ((Integer) miae.getParts()[0]).intValue());
            assertEquals(5, ((Integer) miae.getParts()[1]).intValue());
        }
    }

    private <T extends CalculusFieldElement<T>> void doTestNbParameters(final Field<T> field) {
        final int nbParameters = 6;
        final int nbFunctions  = 3;
        FDSFactory<T> factory = new FDSFactory<>(field, nbParameters, 6);
        @SuppressWarnings("unchecked")
        FieldDerivativeStructure<T>[] functions = (FieldDerivativeStructure<T>[]) Array.newInstance(FieldDerivativeStructure.class,
                                                                                                    nbFunctions);
        for (int i = 0; i < functions.length; ++i) {
            functions[i] = factory.constant(0);
        }
        assertEquals(nbParameters,
                            new FieldTaylorMap<>(MathArrays.buildArray(field, nbParameters), functions).getFreeParameters());
    }

    private <T extends CalculusFieldElement<T>> void doTestNbFunctions(final Field<T> field) {
        final int nbParameters = 6;
        final int nbFunctions  = 3;
        FDSFactory<T> factory = new FDSFactory<>(field, nbParameters, 6);
        @SuppressWarnings("unchecked")
        FieldDerivativeStructure<T>[] functions = (FieldDerivativeStructure<T>[]) Array.newInstance(FieldDerivativeStructure.class,
                                                                                                    nbFunctions);
        for (int i = 0; i < functions.length; ++i) {
            functions[i] = factory.constant(0);
        }
        assertEquals(nbFunctions,
                            new FieldTaylorMap<>(MathArrays.buildArray(field, nbParameters), functions).getNbFunctions());
    }

    private <T extends CalculusFieldElement<T>> void doTestIdentity(final Field<T> field) {
        final FieldTaylorMap<T> map = new FieldTaylorMap<>(field, 7, 3, 4);
        for (int i = 0; i < map.getNbFunctions(); ++i) {

            final FieldDerivativeStructure<T> mi = map.getFunction(i);

            assertEquals(0.0, mi.getValue().getReal(), 1.0e-15);

            int[] orders = new int[7];
            orders[i] = 1;
            int expectedOne = mi.getFactory().getCompiler().getPartialDerivativeIndex(orders);

            for (int j = 0; j < mi.getFactory().getCompiler().getSize(); ++j) {
                assertEquals(j == expectedOne ? 1.0 : 0.0, mi.getAllDerivatives()[j].getReal(), 1.0e-15);
            }

        }
    }

    @SuppressWarnings("unchecked")
    private <T extends CalculusFieldElement<T>> void doTestValue(final Field<T> field) {

        final FDSFactory<T>               factory = new FDSFactory<>(field, 2, 3);
        final FieldDerivativeStructure<T> p0      = factory.variable(0, field.getZero().newInstance( 1.0));
        final FieldDerivativeStructure<T> p1      = factory.variable(1, field.getZero().newInstance(-3.0));
        final FieldDerivativeStructure<T> f0      = p0.sin();
        final FieldDerivativeStructure<T> f1      = p0.add(p1);
        final T[] p = MathArrays.buildArray(field, 2);
        p[0] = p0.getValue();
        p[1] = p1.getValue();
        final FieldDerivativeStructure<T>[] f = (FieldDerivativeStructure<T>[]) Array.newInstance(FieldDerivativeStructure.class, 2);
        f[0] = f0;
        f[1] = f1;
        final FieldTaylorMap<T> map = new FieldTaylorMap<>(p, f);

        for (double dp0 = -0.1; dp0 < 0.1; dp0 += 0.01) {
            final T dp0T = field.getZero().newInstance(dp0);
            for (double dp1 = -0.1; dp1 < 0.1; dp1 += 0.01) {
                final T dp1T = field.getZero().newInstance(dp1);
                assertEquals(f0.taylor(dp0, dp1).getReal(), map.value(dp0, dp1)[0].getReal(),   1.0e-15);
                assertEquals(f1.taylor(dp0, dp1).getReal(), map.value(dp0, dp1)[1].getReal(),   1.0e-15);
                assertEquals(f0.taylor(dp0, dp1).getReal(), map.value(dp0T, dp1T)[0].getReal(), 1.0e-15);
                assertEquals(f1.taylor(dp0, dp1).getReal(), map.value(dp0T, dp1T)[1].getReal(), 1.0e-15);
            }
        }

    }

    private <T extends CalculusFieldElement<T>> void doTestCompose(final Field<T> field) {

        final FDSFactory<T>               factory2 = new FDSFactory<>(field, 2, 2);
        final FDSFactory<T>               factory3 = new FDSFactory<>(field, 3, factory2.getCompiler().getOrder());
        final FieldDerivativeStructure<T> p0       = factory2.variable(0,  1.0);
        final FieldDerivativeStructure<T> p1       = factory2.variable(1, -3.0);
        final FieldDerivativeStructure<T> g0       = p0.sin();
        final FieldDerivativeStructure<T> g1       = p0.add(p1);
        final FieldDerivativeStructure<T> g2       = p1.multiply(p0);
        final FieldDerivativeStructure<T> f0       = factory3.variable(0,  g0.getValue()).
                                                     add(factory3.variable(1, g1.getValue()));
        final FieldDerivativeStructure<T> f1       = factory3.variable(0,  g0.getValue()).
                                                     subtract(factory3.variable(1, g1.getValue())).
                                                     add(factory3.variable(2, g2.getValue()));
        final T[] p = MathArrays.buildArray(field, 2);
        p[0] = p0.getValue();
        p[1] = p1.getValue();
        @SuppressWarnings("unchecked")
        final FieldDerivativeStructure<T>[] g = (FieldDerivativeStructure<T>[]) Array.newInstance(FieldDerivativeStructure.class, 3);
        g[0] = g0;
        g[1] = g1;
        g[2] = g2;
        final T[] gT = MathArrays.buildArray(field, 3);
        gT[0] = g0.getValue();
        gT[1] = g1.getValue();
        gT[2] = g2.getValue();
        @SuppressWarnings("unchecked")
        final FieldDerivativeStructure<T>[] f = (FieldDerivativeStructure<T>[]) Array.newInstance(FieldDerivativeStructure.class, 2);
        f[0] = f0;
        f[1] = f1;
        final FieldTaylorMap<T>           mapG     = new FieldTaylorMap<>(p,  g);
        final FieldTaylorMap<T>           mapF     = new FieldTaylorMap<>(gT, f);
        final FieldTaylorMap<T>           composed = mapF.compose(mapG);

        for (double dp0 = -0.1; dp0 < 0.1; dp0 += 0.01) {
            for (double dp1 = -0.1; dp1 < 0.1; dp1 += 0.01) {
                assertEquals(g0.taylor(dp0, dp1).add(g1.taylor(dp0, dp1)).getReal(),
                                    composed.value(dp0, dp1)[0].getReal(),
                                    1.0e-15);
                assertEquals(g0.taylor(dp0, dp1).subtract(g1.taylor(dp0, dp1)).add(g2.taylor(dp0, dp1)).getReal(),
                                    composed.value(dp0, dp1)[1].getReal(),
                                    1.0e-15);
            }
        }

        assertEquals(p0.getValue().getReal(), mapG.getPoint()[0].getReal(),     1.0e-15);
        assertEquals(p1.getValue().getReal(), mapG.getPoint()[1].getReal(),     1.0e-15);
        assertEquals(g0.getValue().getReal(), mapF.getPoint()[0].getReal(),     1.0e-15);
        assertEquals(g1.getValue().getReal(), mapF.getPoint()[1].getReal(),     1.0e-15);
        assertEquals(g2.getValue().getReal(), mapF.getPoint()[2].getReal(),     1.0e-15);
        assertEquals(p0.getValue().getReal(), composed.getPoint()[0].getReal(), 1.0e-15);
        assertEquals(p1.getValue().getReal(), composed.getPoint()[1].getReal(), 1.0e-15);

        // the partial derivatives of f are only (∂f/∂g₀, ∂f/∂g₁, ∂f/∂g₂)
        assertEquals(+1.0, mapF.getFunction(0).getPartialDerivative(1, 0, 0).getReal(), 1.0e-15);
        assertEquals(+1.0, mapF.getFunction(0).getPartialDerivative(0, 1, 0).getReal(), 1.0e-15);
        assertEquals( 0.0, mapF.getFunction(0).getPartialDerivative(0, 0, 1).getReal(), 1.0e-15);
        assertEquals(+1.0, mapF.getFunction(1).getPartialDerivative(1, 0, 0).getReal(), 1.0e-15);
        assertEquals(-1.0, mapF.getFunction(1).getPartialDerivative(0, 1, 0).getReal(), 1.0e-15);
        assertEquals(+1.0, mapF.getFunction(1).getPartialDerivative(0, 0, 1).getReal(), 1.0e-15);

        assertEquals( 0.0, mapF.getFunction(0).getPartialDerivative(2, 0, 0).getReal(), 1.0e-15);
        assertEquals( 0.0, mapF.getFunction(0).getPartialDerivative(1, 1, 0).getReal(), 1.0e-15);
        assertEquals( 0.0, mapF.getFunction(0).getPartialDerivative(1, 0, 1).getReal(), 1.0e-15);
        assertEquals( 0.0, mapF.getFunction(0).getPartialDerivative(0, 2, 0).getReal(), 1.0e-15);
        assertEquals( 0.0, mapF.getFunction(0).getPartialDerivative(0, 1, 1).getReal(), 1.0e-15);
        assertEquals( 0.0, mapF.getFunction(0).getPartialDerivative(0, 0, 2).getReal(), 1.0e-15);
        assertEquals( 0.0, mapF.getFunction(1).getPartialDerivative(2, 0, 0).getReal(), 1.0e-15);
        assertEquals( 0.0, mapF.getFunction(1).getPartialDerivative(1, 1, 0).getReal(), 1.0e-15);
        assertEquals( 0.0, mapF.getFunction(1).getPartialDerivative(1, 0, 1).getReal(), 1.0e-15);
        assertEquals( 0.0, mapF.getFunction(1).getPartialDerivative(0, 2, 0).getReal(), 1.0e-15);
        assertEquals( 0.0, mapF.getFunction(1).getPartialDerivative(0, 1, 1).getReal(), 1.0e-15);
        assertEquals( 0.0, mapF.getFunction(1).getPartialDerivative(0, 0, 2).getReal(), 1.0e-15);

        // the partial derivatives of the composed map are (∂f/∂p₀, ∂f/∂p₁)
        assertEquals(FastMath.cos(p0.getValue()).getReal() + 1.0,                           composed.getFunction(0).getPartialDerivative(1, 0).getReal(), 1.0e-15);
        assertEquals(+1.0,                                                                  composed.getFunction(0).getPartialDerivative(0, 1).getReal(), 1.0e-15);
        assertEquals(FastMath.cos(p0.getValue()).getReal() - 1.0 + p1.getValue().getReal(), composed.getFunction(1).getPartialDerivative(1, 0).getReal(), 1.0e-15);
        assertEquals(-1.0 + p0.getValue().getReal(),                                        composed.getFunction(1).getPartialDerivative(0, 1).getReal(), 1.0e-15);
        assertEquals(-FastMath.sin(p0.getValue()).getReal(),                                composed.getFunction(0).getPartialDerivative(2, 0).getReal(), 1.0e-15);
        assertEquals( 0.0,                                                                  composed.getFunction(0).getPartialDerivative(1, 1).getReal(), 1.0e-15);
        assertEquals( 0.0,                                                                  composed.getFunction(0).getPartialDerivative(0, 2).getReal(), 1.0e-15);
        assertEquals(-FastMath.sin(p0.getValue()).getReal(),                                composed.getFunction(1).getPartialDerivative(2, 0).getReal(), 1.0e-15);
        assertEquals(+1.0,                                                                  composed.getFunction(1).getPartialDerivative(1, 1).getReal(), 1.0e-15);
        assertEquals( 0.0,                                                                  composed.getFunction(1).getPartialDerivative(0, 2).getReal(), 1.0e-15);

    }

    private <T extends CalculusFieldElement<T>> void doTestInvertNonSquare(final Field<T> field) {
        final FDSFactory<T>           factory   = new FDSFactory<>(field, 2, 2);
        final FieldDerivativeStructure<T> p0        = factory.variable(0,  1.0);
        final FieldDerivativeStructure<T> p1        = factory.variable(1, -3.0);
        final T[] p = MathArrays.buildArray(field, 2);
        p[0] = p0.getValue();
        p[1] = p1.getValue();
        @SuppressWarnings("unchecked")
        final FieldDerivativeStructure<T>[] f = (FieldDerivativeStructure<T>[]) Array.newInstance(FieldDerivativeStructure.class, 3);
        f[0] = p0;
        f[1] = p1;
        f[2] = p0.add(p1);
        final FieldTaylorMap<T> nonSquare = new FieldTaylorMap<>(p, f);
        assertEquals(2, nonSquare.getFreeParameters());
        assertEquals(3, nonSquare.getNbFunctions());
        try {
            nonSquare.invert(new FieldQRDecomposer<>(field.getZero().newInstance(1.0e-10)));
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
            assertEquals(3, ((Integer) miae.getParts()[0]).intValue());
            assertEquals(2, ((Integer) miae.getParts()[1]).intValue());
        }
    }

    private <T extends CalculusFieldElement<T>> void doTestInvertMonoDimensional(final Field<T> field) {
        final FDSFactory<T> factory = new FDSFactory<>(field, 1, 6);
        for (double x = 0.0; x < 3.0; x += 0.01) {
            final FieldDerivativeStructure<T> xDS = factory.variable(0, x);
            final T[] p = MathArrays.buildArray(field, 1);
            p[0] = xDS.getValue();
            @SuppressWarnings("unchecked")
            final FieldDerivativeStructure<T>[] f = (FieldDerivativeStructure<T>[]) Array.newInstance(FieldDerivativeStructure.class, 1);
            f[0] = xDS.exp();
            final FieldTaylorMap<T> expMap = new FieldTaylorMap<>(p, f);
            final FieldTaylorMap<T> inverse = expMap.invert(new FieldQRDecomposer<>(field.getZero().newInstance(1.0e-10)));
            final FieldDerivativeStructure<T> log = factory.variable(0, expMap.getFunction(0).getValue()).log();
            FieldDerivativeStructureAbstractTest.checkEquals(log, inverse.getFunction(0), 4.7e-13);
        }
    }

    private <T extends CalculusFieldElement<T>> void doTestInvertBiDimensional(final Field<T> field) {
        final FDSFactory<T> factory = new FDSFactory<>(field, 2, 4);
        for (double x = -2.0 + FastMath.scalb(1.0, -6); x < 2.0; x += FastMath.scalb(1.0, -5)) {
            final FieldDerivativeStructure<T> xDS = factory.variable(0, x);
            for (double y = -2.0 + FastMath.scalb(1.0, -6); y < 2.0; y += FastMath.scalb(1.0, -5)) {
                final FieldDerivativeStructure<T> yDS = factory.variable(1, y);
                final T[] p = MathArrays.buildArray(field, 2);
                p[0] = xDS.getValue();
                p[1] = yDS.getValue();
                @SuppressWarnings("unchecked")
                final FieldDerivativeStructure<T>[] f = (FieldDerivativeStructure<T>[]) Array.newInstance(FieldDerivativeStructure.class, 2);
                f[0] = FastMath.hypot(xDS, yDS);
                f[1] = FastMath.atan2(yDS, xDS);
                final FieldTaylorMap<T> polarMap = new FieldTaylorMap<>(p, f);
                final FieldTaylorMap<T> cartMap  = polarMap.invert(new FieldQRDecomposer<>(field.getZero().newInstance(1.0e-10)));
                final FieldTaylorMap<T> idMap    = cartMap.compose(polarMap);
                FieldDerivativeStructureAbstractTest.checkEquals(xDS, idMap.getFunction(0), 2.8e-9);
                FieldDerivativeStructureAbstractTest.checkEquals(yDS, idMap.getFunction(1), 2.8e-9);
            }
        }
    }

}
