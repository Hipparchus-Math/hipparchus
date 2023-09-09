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

package org.hipparchus.ode.nonstiff;


import org.hipparchus.CalculusFieldElement;
import org.hipparchus.ode.FieldOrdinaryDifferentialEquation;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.FieldSinCos;
import org.hipparchus.util.MathArrays;

public class Ellipse<T extends CalculusFieldElement<T>> implements FieldOrdinaryDifferentialEquation<T> {

    private final T a;
    private final T b;
    private final T omega;

    Ellipse(final T a, final T b, final T omega) {
        this.a     = a;
        this.b     = b;
        this.omega = omega;
    }

    @Override
    public int getDimension() {
        return 4;
    }

    @Override
    public T[] computeDerivatives(T t, T[] y) {
        final T[] yDot = MathArrays.buildArray(t.getField(), getDimension());
        yDot[0] = y[2];
        yDot[1] = y[3];
        yDot[2] = y[0].multiply(omega.multiply(omega).negate());
        yDot[3] = y[1].multiply(omega.multiply(omega).negate());
        return yDot;
    }

    public T[] computeTheoreticalState(T t) {
        final T[] y = MathArrays.buildArray(t.getField(), getDimension());
        final FieldSinCos<T> sc = FastMath.sinCos(t.multiply(omega));
        y[0] = a.multiply(sc.cos());
        y[1] = b.multiply(sc.sin());
        y[2] = a.multiply(omega).multiply(sc.sin()).negate();
        y[3] = b.multiply(omega).multiply(sc.cos());
        return y;
    }

}
