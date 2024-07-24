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

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.linear.QRDecomposer;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Test for class {@link TaylorMap}.
 */
class TaylorMapTest {

    @Test
    void testNullPoint() {
        try {
            new TaylorMap(null, new DerivativeStructure[2]);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.NUMBER_OF_ELEMENTS_SHOULD_BE_POSITIVE, miae.getSpecifier());
            assertEquals(0, ((Integer) miae.getParts()[0]).intValue());
        }
    }

    @Test
    void testDim0Point() {
        try {
            new TaylorMap(new double[0], new DerivativeStructure[2]);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.NUMBER_OF_ELEMENTS_SHOULD_BE_POSITIVE, miae.getSpecifier());
            assertEquals(0, ((Integer) miae.getParts()[0]).intValue());
        }
    }

    @Test
    void testNullFunctions() {
        try {
            new TaylorMap(new double[2], null);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.NUMBER_OF_ELEMENTS_SHOULD_BE_POSITIVE, miae.getSpecifier());
            assertEquals(0, ((Integer) miae.getParts()[0]).intValue());
        }
    }

    @Test
    void testNoFunctions() {
        try {
            new TaylorMap(new double[2], new DerivativeStructure[0]);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.NUMBER_OF_ELEMENTS_SHOULD_BE_POSITIVE, miae.getSpecifier());
            assertEquals(0, ((Integer) miae.getParts()[0]).intValue());
        }
    }

    @Test
    void testIncompatiblePointAndFunctions() {
        DSFactory factory = new DSFactory(6, 6);
        DerivativeStructure[] functions = new DerivativeStructure[factory.getCompiler().getFreeParameters()];
        for (int i = 0; i < functions.length; ++i) {
            functions[i] = factory.constant(0);
        }
        try {
            new TaylorMap(new double[functions.length - 1], functions);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
            assertEquals(5, ((Integer) miae.getParts()[0]).intValue());
            assertEquals(6, ((Integer) miae.getParts()[1]).intValue());
        }
    }

    @Test
    void testIncompatible() {
        DSFactory factory = new DSFactory(6, 6);
        DerivativeStructure[] functions = new DerivativeStructure[factory.getCompiler().getFreeParameters()];
        for (int i = 0; i < functions.length - 1; ++i) {
            functions[i] = factory.constant(0);
        }
        functions[functions.length - 1] = new DSFactory(factory.getCompiler().getFreeParameters(),
                                                        factory.getCompiler().getOrder() - 1).
                                          constant(1.0);
        try {
            new TaylorMap(new double[functions.length], functions);
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
            assertEquals(6, ((Integer) miae.getParts()[0]).intValue());
            assertEquals(5, ((Integer) miae.getParts()[1]).intValue());
        }
    }

    @Test
    void testNbParameters() {
        final int nbParameters = 6;
        final int nbFunctions  = 3;
        DSFactory factory = new DSFactory(nbParameters, 6);
        DerivativeStructure[] functions = new DerivativeStructure[nbFunctions];
        for (int i = 0; i < functions.length; ++i) {
            functions[i] = factory.constant(0);
        }
        assertEquals(nbParameters,
                            new TaylorMap(new double[nbParameters], functions).getFreeParameters());
    }

    @Test
    void testNbFunctions() {
        final int nbParameters = 6;
        final int nbFunctions  = 3;
        DSFactory factory = new DSFactory(nbParameters, 6);
        DerivativeStructure[] functions = new DerivativeStructure[nbFunctions];
        for (int i = 0; i < functions.length; ++i) {
            functions[i] = factory.constant(0);
        }
        assertEquals(nbFunctions,
                            new TaylorMap(new double[nbParameters], functions).getNbFunctions());
    }

    @Test
    void testIdentity() {
        final TaylorMap map = new TaylorMap(7, 3, 4);
        for (int i = 0; i < map.getNbFunctions(); ++i) {

            final DerivativeStructure mi = map.getFunction(i);

            assertEquals(0.0, mi.getValue(), 1.0e-15);

            int[] orders = new int[7];
            orders[i] = 1;
            int expectedOne = mi.getFactory().getCompiler().getPartialDerivativeIndex(orders);

            for (int j = 0; j < mi.getFactory().getCompiler().getSize(); ++j) {
                assertEquals(j == expectedOne ? 1.0 : 0.0, mi.getAllDerivatives()[j], 1.0e-15);
            }

        }
    }

    @Test
    void testValue() {

        final DSFactory           factory = new DSFactory(2, 3);
        final DerivativeStructure p0      = factory.variable(0,  1.0);
        final DerivativeStructure p1      = factory.variable(1, -3.0);
        final DerivativeStructure f0      = p0.sin();
        final DerivativeStructure f1      = p0.add(p1);
        final TaylorMap           map     = new TaylorMap(new double[] { p0.getValue(), p1.getValue() },
                                                          new DerivativeStructure[] { f0, f1 });

        for (double dp0 = -0.1; dp0 < 0.1; dp0 += 0.01) {
            for (double dp1 = -0.1; dp1 < 0.1; dp1 += 0.01) {
                assertEquals(f0.taylor(dp0, dp1), map.value(dp0, dp1)[0], 1.0e-15);
                assertEquals(f1.taylor(dp0, dp1), map.value(dp0, dp1)[1], 1.0e-15);
            }
        }

    }

    @Test
    void testCompose() {

        final DSFactory           factory2 = new DSFactory(2, 2);
        final DSFactory           factory3 = new DSFactory(3, factory2.getCompiler().getOrder());
        final DerivativeStructure p0       = factory2.variable(0,  1.0);
        final DerivativeStructure p1       = factory2.variable(1, -3.0);
        final DerivativeStructure g0       = p0.sin();
        final DerivativeStructure g1       = p0.add(p1);
        final DerivativeStructure g2       = p1.multiply(p0);
        final DerivativeStructure f0       = factory3.variable(0,  g0.getValue()).
                                             add(factory3.variable(1, g1.getValue()));
        final DerivativeStructure f1       = factory3.variable(0,  g0.getValue()).
                                             subtract(factory3.variable(1, g1.getValue())).
                                             add(factory3.variable(2, g2.getValue()));
        final TaylorMap           mapG     = new TaylorMap(new double[] { p0.getValue(), p1.getValue() },
                                                           new DerivativeStructure[] { g0, g1, g2 });
        final TaylorMap           mapF     = new TaylorMap(new double[] { g0.getValue(), g1.getValue(), g2.getValue() },
                                                           new DerivativeStructure[] { f0, f1 });
        final TaylorMap           composed = mapF.compose(mapG);

        for (double dp0 = -0.1; dp0 < 0.1; dp0 += 0.01) {
            for (double dp1 = -0.1; dp1 < 0.1; dp1 += 0.01) {
                assertEquals(g0.taylor(dp0, dp1) + g1.taylor(dp0, dp1),
                                    composed.value(dp0, dp1)[0],
                                    1.0e-15);
                assertEquals(g0.taylor(dp0, dp1) - g1.taylor(dp0, dp1) + g2.taylor(dp0, dp1),
                                    composed.value(dp0, dp1)[1],
                                    1.0e-15);
            }
        }

        assertEquals(p0.getValue(), mapG.getPoint()[0],     1.0e-15);
        assertEquals(p1.getValue(), mapG.getPoint()[1],     1.0e-15);
        assertEquals(g0.getValue(), mapF.getPoint()[0],     1.0e-15);
        assertEquals(g1.getValue(), mapF.getPoint()[1],     1.0e-15);
        assertEquals(g2.getValue(), mapF.getPoint()[2],     1.0e-15);
        assertEquals(p0.getValue(), composed.getPoint()[0], 1.0e-15);
        assertEquals(p1.getValue(), composed.getPoint()[1], 1.0e-15);

        // the partial derivatives of f are only (∂f/∂g₀, ∂f/∂g₁, ∂f/∂g₂)
        assertEquals(+1.0, mapF.getFunction(0).getPartialDerivative(1, 0, 0), 1.0e-15);
        assertEquals(+1.0, mapF.getFunction(0).getPartialDerivative(0, 1, 0), 1.0e-15);
        assertEquals( 0.0, mapF.getFunction(0).getPartialDerivative(0, 0, 1), 1.0e-15);
        assertEquals(+1.0, mapF.getFunction(1).getPartialDerivative(1, 0, 0), 1.0e-15);
        assertEquals(-1.0, mapF.getFunction(1).getPartialDerivative(0, 1, 0), 1.0e-15);
        assertEquals(+1.0, mapF.getFunction(1).getPartialDerivative(0, 0, 1), 1.0e-15);

        assertEquals( 0.0, mapF.getFunction(0).getPartialDerivative(2, 0, 0), 1.0e-15);
        assertEquals( 0.0, mapF.getFunction(0).getPartialDerivative(1, 1, 0), 1.0e-15);
        assertEquals( 0.0, mapF.getFunction(0).getPartialDerivative(1, 0, 1), 1.0e-15);
        assertEquals( 0.0, mapF.getFunction(0).getPartialDerivative(0, 2, 0), 1.0e-15);
        assertEquals( 0.0, mapF.getFunction(0).getPartialDerivative(0, 1, 1), 1.0e-15);
        assertEquals( 0.0, mapF.getFunction(0).getPartialDerivative(0, 0, 2), 1.0e-15);
        assertEquals( 0.0, mapF.getFunction(1).getPartialDerivative(2, 0, 0), 1.0e-15);
        assertEquals( 0.0, mapF.getFunction(1).getPartialDerivative(1, 1, 0), 1.0e-15);
        assertEquals( 0.0, mapF.getFunction(1).getPartialDerivative(1, 0, 1), 1.0e-15);
        assertEquals( 0.0, mapF.getFunction(1).getPartialDerivative(0, 2, 0), 1.0e-15);
        assertEquals( 0.0, mapF.getFunction(1).getPartialDerivative(0, 1, 1), 1.0e-15);
        assertEquals( 0.0, mapF.getFunction(1).getPartialDerivative(0, 0, 2), 1.0e-15);

        // the partial derivatives of the composed map are (∂f/∂p₀, ∂f/∂p₁)
        assertEquals(FastMath.cos(p0.getValue()) + 1.0,                 composed.getFunction(0).getPartialDerivative(1, 0), 1.0e-15);
        assertEquals(+1.0,                                              composed.getFunction(0).getPartialDerivative(0, 1), 1.0e-15);
        assertEquals(FastMath.cos(p0.getValue()) - 1.0 + p1.getValue(), composed.getFunction(1).getPartialDerivative(1, 0), 1.0e-15);
        assertEquals(-1.0 + p0.getValue(),                              composed.getFunction(1).getPartialDerivative(0, 1), 1.0e-15);
        assertEquals(-FastMath.sin(p0.getValue()),                      composed.getFunction(0).getPartialDerivative(2, 0), 1.0e-15);
        assertEquals( 0.0,                                              composed.getFunction(0).getPartialDerivative(1, 1), 1.0e-15);
        assertEquals( 0.0,                                              composed.getFunction(0).getPartialDerivative(0, 2), 1.0e-15);
        assertEquals(-FastMath.sin(p0.getValue()),                      composed.getFunction(1).getPartialDerivative(2, 0), 1.0e-15);
        assertEquals(+1.0,                                              composed.getFunction(1).getPartialDerivative(1, 1), 1.0e-15);
        assertEquals( 0.0,                                              composed.getFunction(1).getPartialDerivative(0, 2), 1.0e-15);

    }

    @Test
    void testInvertNonSquare() {
        final DSFactory           factory   = new DSFactory(2, 2);
        final DerivativeStructure p0        = factory.variable(0,  1.0);
        final DerivativeStructure p1        = factory.variable(1, -3.0);
        final TaylorMap           nonSquare = new TaylorMap(new double[] { p0.getValue(), p1.getValue() },
                                                            new DerivativeStructure[] { p0, p1, p0.add(p1) });
        assertEquals(2, nonSquare.getFreeParameters());
        assertEquals(3, nonSquare.getNbFunctions());
        try {
            nonSquare.invert(new QRDecomposer(1.0e-10));
            fail("an exception should have been thrown");
        } catch (MathIllegalArgumentException miae) {
            assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, miae.getSpecifier());
            assertEquals(3, ((Integer) miae.getParts()[0]).intValue());
            assertEquals(2, ((Integer) miae.getParts()[1]).intValue());
        }
    }

    @Test
    void testInvertMonoDimensional() {
        final DSFactory factory = new DSFactory(1, 6);
        for (double x = 0.0; x < 3.0; x += 0.01) {
            final DerivativeStructure xDS = factory.variable(0, x);
            final TaylorMap expMap = new TaylorMap(new double[] { xDS.getValue() },
                                                   new DerivativeStructure[]  { xDS.exp() });
            final TaylorMap inverse = expMap.invert(new QRDecomposer(1.0e-10));
            final DerivativeStructure log = factory.variable(0, expMap.getFunction(0).getValue()).log();
            DerivativeStructureTest.checkEquals(log, inverse.getFunction(0), 4.7e-13);
        }
    }

    @Test
    void testInvertBiDimensional() {
        final DSFactory factory = new DSFactory(2, 4);
        for (double x = -2.0 + FastMath.scalb(1.0, -6); x < 2.0; x += FastMath.scalb(1.0, -5)) {
            final DerivativeStructure xDS = factory.variable(0, x);
            for (double y = -2.0 + FastMath.scalb(1.0, -6); y < 2.0; y += FastMath.scalb(1.0, -5)) {
                final DerivativeStructure yDS = factory.variable(1, y);
                final TaylorMap polarMap = new TaylorMap(new double[] { xDS.getValue(), yDS.getValue() },
                                                         new DerivativeStructure[]  { FastMath.hypot(xDS, yDS), FastMath.atan2(yDS, xDS)});
                final TaylorMap cartMap  = polarMap.invert(new QRDecomposer(1.0e-10));
                final TaylorMap idMap    = cartMap.compose(polarMap);
                DerivativeStructureTest.checkEquals(xDS, idMap.getFunction(0), 2.8e-9);
                DerivativeStructureTest.checkEquals(yDS, idMap.getFunction(1), 2.8e-9);
            }
        }
    }

}
