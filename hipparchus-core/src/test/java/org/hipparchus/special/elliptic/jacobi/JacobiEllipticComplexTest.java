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
package org.hipparchus.special.elliptic.jacobi;

import java.io.IOException;

import org.hipparchus.complex.Complex;
import org.junit.Assert;
import org.junit.Test;

public class JacobiEllipticComplexTest {

    @Test
    public void testComplex() throws IOException {
        final Complex z = Complex.valueOf(5.2, -2.5);
        final FieldJacobiElliptic<Complex> je = JacobiEllipticBuilder.build(Complex.valueOf(0.3, 1.0));
        final FieldCopolarC<Complex> valuesC = je.valuesC(z);
        Assert.assertEquals(-0.24609405083573348938, valuesC.sc().getRealPart(),      1.0e-15);
        Assert.assertEquals( 0.74202229271111558523, valuesC.sc().getImaginaryPart(), 1.0e-15);
    }

}
