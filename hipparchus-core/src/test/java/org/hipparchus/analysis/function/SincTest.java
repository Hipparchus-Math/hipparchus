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
package org.hipparchus.analysis.function;

import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.differentiation.DSFactory;
import org.hipparchus.analysis.differentiation.DerivativeStructure;
import org.hipparchus.dfp.Dfp;
import org.hipparchus.dfp.DfpField;
import org.hipparchus.dfp.DfpMath;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

class SincTest {

    @Test
    void testShortcut() {
       final Sinc s = new Sinc();
       final UnivariateFunction f = new UnivariateFunction() {
           @Override
           public double value(double x) {
               Dfp dfpX = new DfpField(25).newDfp(x);
               return DfpMath.sin(dfpX).divide(dfpX).toDouble();
           }
       };

       for (double x = 1e-30; x < 1e10; x *= 2) {
           final double fX = f.value(x);
           final double sX = s.value(x);
           assertEquals(fX, sX, 2.0e-16, "x=" + x);
       }
   }

    @Test
    void testCrossings() {
       final Sinc s = new Sinc(true);
       final int numCrossings = 1000;
       final double tol = 2e-16;
       for (int i = 1; i <= numCrossings; i++) {
           assertEquals(0, s.value(i), tol, "i=" + i);
       }
   }

    @Test
    void testZero() {
       final Sinc s = new Sinc();
       assertEquals(1d, s.value(0), 0);
   }

    @Test
    void testEuler() {
       final Sinc s = new Sinc();
       final double x = 123456.789;
       double prod = 1;
       double xOverPow2 = x / 2;
       while (xOverPow2 > 0) {
           prod *= FastMath.cos(xOverPow2);
           xOverPow2 /= 2;
       }
       assertEquals(prod, s.value(x), 1e-13);
   }

    @Test
    void testDerivativeZero() {
       final DerivativeStructure s0 = new Sinc(true).value(new DSFactory(1, 1).variable(0, 0.0));
       assertEquals(0, s0.getPartialDerivative(1), 0);
   }

    @Test
    void testDerivatives1Dot2Unnormalized() {
       DerivativeStructure s = new Sinc(false).value(new DSFactory(1, 5).variable(0, 1.2));
       assertEquals( 0.77669923830602195806, s.getPartialDerivative(0), 1.0e-16);
       assertEquals(-0.34528456985779031701, s.getPartialDerivative(1), 1.0e-16);
       assertEquals(-0.2012249552097047631,  s.getPartialDerivative(2), 1.0e-16);
       assertEquals( 0.2010975926270339262,  s.getPartialDerivative(3), 4.0e-16);
       assertEquals( 0.106373929549242204,   s.getPartialDerivative(4), 1.0e-15);
       assertEquals(-0.1412599110579478695,  s.getPartialDerivative(5), 3.0e-15);
   }

    @Test
    void testDerivatives1Dot2Normalized() {
       DerivativeStructure s = new Sinc(true).value(new DSFactory(1, 5).variable(0, 1.2));
       assertEquals(-0.15591488063143983888, s.getPartialDerivative(0), 6.0e-17);
       assertEquals(-0.54425176145292298767, s.getPartialDerivative(1), 2.0e-16);
       assertEquals(2.4459044611635856107,   s.getPartialDerivative(2), 9.0e-16);
       assertEquals(0.5391369206235909586,   s.getPartialDerivative(3), 7.0e-16);
       assertEquals(-16.984649869728849865,  s.getPartialDerivative(4), 8.0e-15);
       assertEquals(5.0980327462666316586,   s.getPartialDerivative(5), 9.0e-15);
   }

    @Test
    void testDerivativeShortcut() {
       final Sinc sinc = new Sinc();
       final UnivariateFunction f = new UnivariateFunction() {
               @Override
            public double value(double x) {
                   Dfp dfpX = new DfpField(25).newDfp(x);
                   return DfpMath.cos(dfpX).subtract(DfpMath.sin(dfpX).divide(dfpX)).divide(dfpX).toDouble();
               }
           };

       DSFactory factory = new DSFactory(1, 1);
       for (double x = 1e-30; x < 1e10; x *= 2) {
           final double fX = f.value(x);
           final DerivativeStructure sX = sinc.value(factory.variable(0, x));
           assertEquals(fX, sX.getPartialDerivative(1), 3.0e-13, "x=" + x);
       }
   }
}
