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

package org.hipparchus.ode;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.special.elliptic.jacobi.FieldCopolarN;
import org.hipparchus.special.elliptic.jacobi.FieldJacobiElliptic;
import org.hipparchus.special.elliptic.jacobi.JacobiEllipticBuilder;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;

/**
 * This class is used in the junit tests for the ODE integrators.

 * <p>This specific problem correspond to torque-free motion of a solid body
 * with moments of inertia I₁, I₂, and I₃ with respect to body axes x, y, and z.
 * We use here the notations from Landau and Lifchitz Course of Theoretical Physics,
 * Mechanics vol 1.
 * </p>
 * <p>
 * The equations of torque-free motion are given in the solid body frame by
 * equation 36.5:
 * <pre>
 *    I₁ dΩ₁/dt + (I₃ - I₂) Ω₂ Ω₃ = 0
 *    I₂ dΩ₂/dt + (I₁ - I₃) Ω₃ Ω₁ = 0
 *    I₃ dΩ₃/dt + (I₂ - I₁) Ω₁ Ω₂ = 0
 * </pre>
 * <p>
 * The moments of inertia and initial conditions are: I₁ = 3/8, I₂ = 1/2, I₃ = 5/8,
 * Ω₁ = 5, Ω₂ = 0, Ω₃ = 4. This corresponds to a motion with angular velocity
 * describing a large polhode around Z axis in solid body frame. The motion is almost
 * unstable as M² is only slightly greater than 2EI₂. Increasing Ω₁ to √(80/3) ≈ 5.16398
 * would imply M²=2EI₂, and the polhode would degenerate to two intersecting ellipses.
 * </p>
 * <p>
 * The torque-free motion can be solved analytically using Jacobi elliptic functions
 * (in the rotating body frame)
 * <pre>
 *   τ      = t √([I₃-I₂][M²-2EI₁]/[I₁I₂I₃])
 *   Ω₁ (τ) = √([2EI₃-M²]/[I₁(I₃-I₁)]) cn(τ)
 *   Ω₂ (τ) = √([2EI₃-M²]/[I₂(I₃-I₂)]) sn(τ)
 *   Ω₃ (τ) = √([M²-2EI₁]/[I₃(I₃-I₁)]) dn(τ)
 * </pre>
 * </p>
 * <p>
 * This problem solves only solves the rotation rate part, whereas {@code FieldTestProblem8}
 * solves the full motion (rotation rate and rotation).
 * </p>
 * @param <T> the type of the field elements
 */
public class TestFieldProblem7<T extends CalculusFieldElement<T>>
    extends TestFieldProblemAbstract<T> {

    /** Moments of inertia. */
    final T i1;
    final T i2;
    final T i3;

    /** Twice the angular kinetic energy. */
    final T twoE;

    final T m2;

    /** Time scaling factor. */
    final T tScale;

    /** State scaling factors. */
    final T o1Scale;
    final T o2Scale;
    final T o3Scale;

    /** Jacobi elliptic function. */
    final FieldJacobiElliptic<T> jacobi;

    /**
     * Simple constructor.
     * @param field field to which elements belong
     */
    public TestFieldProblem7(Field<T> field) {
        super(convert(field, 0.0),
              convert(field, new double[] { 5.0, 0.0, 4.0 }),
              convert(field, 4.0),
              convert(field, new double[] { 1.0, 1.0, 1.0 }));
        i1 = convert(field, 3.0 / 8.0);
        i2 = convert(field, 1.0 / 2.0);
        i3 = convert(field, 5.0 / 8.0);

        final T[] s0 = getInitialState().getPrimaryState();
        final T o12 = s0[0].multiply(s0[0]);
        final T o22 = s0[1].multiply(s0[1]);
        final T o32 = s0[2].multiply(s0[2]);
        twoE    =  i1.multiply(o12).add(i2.multiply(o22)).add(i3.multiply(o32));
        m2      =  i1.multiply(i1.multiply(o12)).add(i2.multiply(i2.multiply(o22))).add(i3.multiply(i3.multiply(o32)));
        tScale  = FastMath.sqrt(i3.subtract(i2).multiply(m2.subtract(twoE.multiply(i1))).divide(i1.multiply(i2).multiply(i3)));
        o1Scale = FastMath.sqrt(twoE.multiply(i3).subtract(m2).divide(i1.multiply(i3.subtract(i1))));
        o2Scale = FastMath.sqrt(twoE.multiply(i3).subtract(m2).divide(i2.multiply(i3.subtract(i2))));
        o3Scale = FastMath.sqrt(m2.subtract(twoE.multiply(i1)).divide(i3.multiply(i3.subtract(i1))));

        final T k2 = i2.subtract(i1).multiply(twoE.multiply(i3).subtract(m2)).divide(i3.subtract(i2).multiply(m2.subtract(twoE.multiply(i1))));
        jacobi = JacobiEllipticBuilder.build(k2);
    }

    @Override
    public T[] doComputeDerivatives(T t, T[] y) {

        final T[] yDot = MathArrays.buildArray(getField(), getDimension());

        // compute the derivatives using Euler equations
        yDot[0] = y[1].multiply(y[2]).multiply(i2.subtract(i3)).divide(i1);
        yDot[1] = y[2].multiply(y[0]).multiply(i3.subtract(i1)).divide(i2);
        yDot[2] = y[0].multiply(y[1]).multiply(i1.subtract(i2)).divide(i3);

        return yDot;

    }

    @Override
    public T[] computeTheoreticalState(T t) {

        final T[] y = MathArrays.buildArray(getField(), getDimension());

        final FieldCopolarN<T> valuesN = jacobi.valuesN(t.multiply(tScale));
        y[0] = o1Scale.multiply(valuesN.cn());
        y[1] = o2Scale.multiply(valuesN.sn());
        y[2] = o3Scale.multiply(valuesN.dn());

        return y;

    }

}
