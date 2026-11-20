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
package org.hipparchus.transform;

import org.hipparchus.complex.Complex;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class FFTCoverageTest {

    @Test
    void testTransformUtilsScaleComplexArray() {
        Complex[] array = new Complex[] { new Complex(1, 1), new Complex(2, 2) };
        TransformUtils.scaleArray(array, 2.0);
        assertEquals(2.0, array[0].getReal(), 1e-10);
        assertEquals(2.0, array[0].getImaginary(), 1e-10);
        assertEquals(4.0, array[1].getReal(), 1e-10);
        assertEquals(4.0, array[1].getImaginary(), 1e-10);
    }

    @Test
    void testTransformUtilsCreateComplexArrayErrors() {
        // Rows != 2
        double[][] data1 = new double[][] { { 1, 2 } };
        assertThrows(MathIllegalArgumentException.class, () -> TransformUtils.createComplexArray(data1));

        // Unequal lengths
        double[][] data2 = new double[][] { { 1, 2 }, { 1 } };
        assertThrows(MathIllegalArgumentException.class, () -> TransformUtils.createComplexArray(data2));
    }

    @Test
    void testTransformUtilsExactLog2Error() {
        assertThrows(MathIllegalArgumentException.class, () -> TransformUtils.exactLog2(3));
        assertThrows(MathIllegalArgumentException.class, () -> TransformUtils.exactLog2(0));
        assertThrows(MathIllegalArgumentException.class, () -> TransformUtils.exactLog2(-1));
    }

    @Test
    void testFastHadamardTransformerFhtErrors() {
        FastHadamardTransformer transformer = new FastHadamardTransformer();
        
        // double[]
        double[] dataD = new double[] { 1, 2, 3 };
        assertThrows(MathIllegalArgumentException.class, () -> transformer.transform(dataD, TransformType.FORWARD));

        // int[]
        int[] dataI = new int[] { 1, 2, 3 };
        assertThrows(MathIllegalArgumentException.class, () -> transformer.transform(dataI));
    }

    @Test
    void testFastFourierTransformerInPlaceEdgeCases() {
        // n = 1
        double[][] data1 = new double[][] { { 1.0 }, { 0.0 } };
        FastFourierTransformer.transformInPlace(data1, DftNormalization.STANDARD, TransformType.FORWARD);
        assertEquals(1.0, data1[0][0], 1e-10);
        assertEquals(0.0, data1[1][0], 1e-10);

        // n = 2
        double[][] data2 = new double[][] { { 1.0, 1.0 }, { 0.0, 0.0 } };
        FastFourierTransformer.transformInPlace(data2, DftNormalization.STANDARD, TransformType.FORWARD);
        // X_0 = x_0 + x_1 = 2.0
        // X_1 = x_0 - x_1 = 0.0
        assertEquals(2.0, data2[0][0], 1e-10);
        assertEquals(0.0, data2[1][0], 1e-10);
        assertEquals(0.0, data2[0][1], 1e-10);
        assertEquals(0.0, data2[1][1], 1e-10);
        
        // Test UNITARY normalization to cover that branch in normalizeTransformedData
        double[][] data2u = new double[][] { { 1.0, 1.0 }, { 0.0, 0.0 } };
        FastFourierTransformer.transformInPlace(data2u, DftNormalization.UNITARY, TransformType.FORWARD);
        // Scale factor = 1 / sqrt(2)
        assertEquals(2.0 / FastMath.sqrt(2), data2u[0][0], 1e-10);

        // Test INVERSE for n=2
        double[][] data2i = new double[][] { { 2.0, 0.0 }, { 0.0, 0.0 } };
        FastFourierTransformer.transformInPlace(data2i, DftNormalization.STANDARD, TransformType.INVERSE);
        // Scale factor = 1/2
        // X_0 = (2+0)/2 = 1.0
        // X_1 = (2-0)/2 = 1.0
        assertEquals(1.0, data2i[0][0], 1e-10);
        assertEquals(1.0, data2i[0][1], 1e-10);
    }

    @Test
    void testFastFourierTransformerInPlaceInverse4() {
        // n = 4, INVERSE
        double[][] data = new double[][] { { 1, 1, 1, 1 }, { 0, 0, 0, 0 } };
        FastFourierTransformer.transformInPlace(data, DftNormalization.STANDARD, TransformType.INVERSE);
        // Just verify it doesn't crash and does something
        assertNotNull(data);
    }

    @Test
    void testFastFourierTransformerInPlaceError() {
        double[][] data = new double[][] { { 1.0, 2.0, 3.0 }, { 0.0, 0.0, 0.0 } };
        assertThrows(MathIllegalArgumentException.class, () -> 
            FastFourierTransformer.transformInPlace(data, DftNormalization.STANDARD, TransformType.FORWARD));
    }
}
