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

import org.hipparchus.special.elliptic.jacobi.CopolarN;
import org.hipparchus.special.elliptic.jacobi.JacobiElliptic;
import org.hipparchus.special.elliptic.jacobi.JacobiEllipticBuilder;
import org.hipparchus.util.FastMath;

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
 * This problem solves only solves the rotation rate part, whereas {@code TestProblem8}
 * solves the full motion (rotation rate and rotation).
 * </p>

 */
public class TestProblem7 extends TestProblemAbstract {

    /** Moments of inertia. */
    final double i1;
    final double i2;
    final double i3;

    /** Twice the angular kinetic energy. */
    final double twoE;

    final double m2;

    /** Time scaling factor. */
    final double tScale;

    /** State scaling factors. */
    final double o1Scale;
    final double o2Scale;
    final double o3Scale;

    /** Jacobi elliptic function. */
    final JacobiElliptic jacobi;

    /**
     * Simple constructor.
     */
    public TestProblem7() {

        super(0.0, new double[] { 5.0, 0.0, 4.0 }, 4.0, new double[] { 1.0, 1.0, 1.0 });
        i1 = 3.0 / 8.0;
        i2 = 1.0 / 2.0;
        i3 = 5.0 / 8.0;

        final double[] s0 = getInitialState().getPrimaryState();
        final double o12 = s0[0] * s0[0];
        final double o22 = s0[1] * s0[1];
        final double o32 = s0[2] * s0[2];
        twoE    =  i1 * o12 + i2 * o22 + i3 * o32;
        m2      =  i1 * i1 * o12 + i2 * i2 * o22 + i3 * i3 * o32;
        tScale  = FastMath.sqrt((i3 - i2) * (m2 - twoE * i1) / (i1 * i2 * i3));
        o1Scale = FastMath.sqrt((twoE * i3 - m2) / (i1 * (i3 - i1)));
        o2Scale = FastMath.sqrt((twoE * i3 - m2) / (i2 * (i3 - i2)));
        o3Scale = FastMath.sqrt((m2 - twoE * i1) / (i3 * (i3 - i1)));

        final double k2 = (i2 - i1) * (twoE * i3 - m2) / ((i3 - i2) * (m2 - twoE * i1));
        jacobi = JacobiEllipticBuilder.build(k2);

    }

    @Override
    public double[] doComputeDerivatives(double t, double[] y) {

        final  double[] yDot = new double[getDimension()];

        // compute the derivatives using Euler equations
        yDot[0] = y[1] * y[2] * (i2 - i3) / i1;
        yDot[1] = y[2] * y[0] * (i3 - i1) / i2;
        yDot[2] = y[0] * y[1] * (i1 - i2) / i3;

        return yDot;

    }

    @Override
    public double[] computeTheoreticalState(double t) {
        final CopolarN valuesN = jacobi.valuesN(t * tScale);
        return new double[] {
            o1Scale * valuesN.cn(),
            o2Scale * valuesN.sn(),
            o3Scale * valuesN.dn()
        };
    }

}
