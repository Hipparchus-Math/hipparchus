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
import org.hipparchus.analysis.differentiation.UnivariateDifferentiableFunction;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

public class SqrtTest {
   @Test
   public void testComparison() {
       final Sqrt s = new Sqrt();
       final UnivariateFunction f = new UnivariateFunction() {
           @Override
           public double value(double x) {
               return FastMath.sqrt(x);
           }
       };

       for (double x = 1e-30; x < 1e10; x *= 2) {
           final double fX = f.value(x);
           final double sX = s.value(x);
           Assertions.assertEquals(fX, sX, 0, "x=" + x);
       }
   }

   @Test
   public void testDerivativeComparison() {
       final UnivariateDifferentiableFunction sPrime = new Sqrt();
       final UnivariateFunction f = new UnivariateFunction() {
               @Override
            public double value(double x) {
                   return 1 / (2 * FastMath.sqrt(x));
               }
           };

       DSFactory factory = new DSFactory(1, 1);
       for (double x = 1e-30; x < 1e10; x *= 2) {
           final double fX = f.value(x);
           final double sX = sPrime.value(factory.variable(0, x)).getPartialDerivative(1);
           Assertions.assertEquals(fX, sX, FastMath.ulp(fX), "x=" + x);
       }
   }

   @Test
   public void testDerivativesHighOrder() {
       DerivativeStructure s = new Sqrt().value(new DSFactory(1, 5).variable(0, 1.2));
       Assertions.assertEquals(1.0954451150103322269, s.getPartialDerivative(0), 1.0e-16);
       Assertions.assertEquals(0.45643546458763842789, s.getPartialDerivative(1), 1.0e-16);
       Assertions.assertEquals(-0.1901814435781826783,  s.getPartialDerivative(2), 1.0e-16);
       Assertions.assertEquals(0.23772680447272834785,  s.getPartialDerivative(3), 1.0e-16);
       Assertions.assertEquals(-0.49526417598485072465,   s.getPartialDerivative(4), 5.0e-16);
       Assertions.assertEquals(1.4445205132891479465,  s.getPartialDerivative(5), 7.0e-16);
   }

}
