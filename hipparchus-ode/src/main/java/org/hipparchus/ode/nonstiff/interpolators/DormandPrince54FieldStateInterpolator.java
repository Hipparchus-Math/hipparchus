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

package org.hipparchus.ode.nonstiff.interpolators;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.ode.FieldEquationsMapper;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.ode.nonstiff.DormandPrince54Integrator;

/**
 * This class represents an interpolator over the last step during an
 * ODE integration for the 5(4) Dormand-Prince integrator.
 *
 * @see DormandPrince54Integrator
 *
 * @param <T> the type of the field elements
 */

public class DormandPrince54FieldStateInterpolator<T extends CalculusFieldElement<T>>
    extends RungeKuttaFieldStateInterpolator<T> {

    /** Last row of the Butcher-array internal weights, element 0. */
    private final T a70;

    // element 1 is zero, so it is neither stored nor used

    /** Last row of the Butcher-array internal weights, element 2. */
    private final T a72;

    /** Last row of the Butcher-array internal weights, element 3. */
    private final T a73;

    /** Last row of the Butcher-array internal weights, element 4. */
    private final T a74;

    /** Last row of the Butcher-array internal weights, element 5. */
    private final T a75;

    /** Shampine (1986) Dense output, element 0. */
    private final T d0;

    // element 1 is zero, so it is neither stored nor used

    /** Shampine (1986) Dense output, element 2. */
    private final T d2;

    /** Shampine (1986) Dense output, element 3. */
    private final T d3;

    /** Shampine (1986) Dense output, element 4. */
    private final T d4;

    /** Shampine (1986) Dense output, element 5. */
    private final T d5;

    /** Shampine (1986) Dense output, element 6. */
    private final T d6;

    /** Simple constructor.
     * @param field field to which the time and state vector elements belong
     * @param forward integration direction indicator
     * @param yDotK slopes at the intermediate points
     * @param globalPreviousState start of the global step
     * @param globalCurrentState end of the global step
     * @param softPreviousState start of the restricted step
     * @param softCurrentState end of the restricted step
     * @param mapper equations mapper for the all equations
     */
    public DormandPrince54FieldStateInterpolator(final Field<T> field, final boolean forward,
                                                 final T[][] yDotK,
                                                 final FieldODEStateAndDerivative<T> globalPreviousState,
                                                 final FieldODEStateAndDerivative<T> globalCurrentState,
                                                 final FieldODEStateAndDerivative<T> softPreviousState,
                                                 final FieldODEStateAndDerivative<T> softCurrentState,
                                                 final FieldEquationsMapper<T> mapper) {
        super(field, forward, yDotK, globalPreviousState, globalCurrentState, softPreviousState, softCurrentState,
                mapper);
        final T one = field.getOne();
        a70 = one.newInstance(   35.0 /  384.0);
        a72 = one.newInstance(  500.0 / 1113.0);
        a73 = one.newInstance(  125.0 /  192.0);
        a74 = one.newInstance(-2187.0 / 6784.0);
        a75 = one.newInstance(   11.0 /   84.0);
        d0  = one.newInstance(-12715105075.0 /  11282082432.0);
        d2  = one.newInstance( 87487479700.0 /  32700410799.0);
        d3  = one.newInstance(-10690763975.0 /   1880347072.0);
        d4  = one.newInstance(701980252875.0 / 199316789632.0);
        d5  = one.newInstance( -1453857185.0 /    822651844.0);
        d6  = one.newInstance(    69997945.0 /     29380423.0);
    }

    /** {@inheritDoc} */
    @Override
    protected DormandPrince54FieldStateInterpolator<T> create(final Field<T> newField, final boolean newForward, final T[][] newYDotK,
                                                              final FieldODEStateAndDerivative<T> newGlobalPreviousState,
                                                              final FieldODEStateAndDerivative<T> newGlobalCurrentState,
                                                              final FieldODEStateAndDerivative<T> newSoftPreviousState,
                                                              final FieldODEStateAndDerivative<T> newSoftCurrentState,
                                                              final FieldEquationsMapper<T> newMapper) {
        return new DormandPrince54FieldStateInterpolator<>(newField, newForward, newYDotK,
                                                            newGlobalPreviousState, newGlobalCurrentState,
                                                            newSoftPreviousState, newSoftCurrentState,
                                                            newMapper);
    }
    /** {@inheritDoc} */
    @SuppressWarnings("unchecked")
    @Override
    protected FieldODEStateAndDerivative<T> computeInterpolatedStateAndDerivatives(final FieldEquationsMapper<T> mapper,
                                                                                   final T time, final T theta,
                                                                                   final T thetaH, final T oneMinusThetaH) {

        // interpolate
        final T one      = time.getField().getOne();
        final T eta      = one.subtract(theta);
        final T twoTheta = theta.multiply(2);
        final T dot2     = one.subtract(twoTheta);
        final T dot3     = theta.multiply(theta.multiply(-3).add(2));
        final T dot4     = twoTheta.multiply(theta.multiply(twoTheta.subtract(3)).add(1));
        final T[] interpolatedState;
        final T[] interpolatedDerivatives;
        if (getGlobalPreviousState() != null && theta.getReal() <= 0.5) {
            final T f1        = thetaH;
            final T f2        = f1.multiply(eta);
            final T f3        = f2.multiply(theta);
            final T f4        = f3.multiply(eta);
            final T coeff0    = f1.multiply(a70).
                                subtract(f2.multiply(a70.subtract(1))).
                                add(f3.multiply(a70.multiply(2).subtract(1))).
                                add(f4.multiply(d0));
            final T coeff1    = time.getField().getZero();
            final T coeff2    = f1.multiply(a72).
                                subtract(f2.multiply(a72)).
                                add(f3.multiply(a72.multiply(2))).
                                add(f4.multiply(d2));
            final T coeff3    = f1.multiply(a73).
                                subtract(f2.multiply(a73)).
                                add(f3.multiply(a73.multiply(2))).
                                add(f4.multiply(d3));
            final T coeff4    = f1.multiply(a74).
                                subtract(f2.multiply(a74)).
                                add(f3.multiply(a74.multiply(2))).
                                add(f4.multiply(d4));
            final T coeff5    = f1.multiply(a75).
                                subtract(f2.multiply(a75)).
                                add(f3.multiply(a75.multiply(2))).
                                add(f4.multiply(d5));
            final T coeff6    = f4.multiply(d6).subtract(f3);
            final T coeffDot0 = a70.
                                subtract(dot2.multiply(a70.subtract(1))).
                                add(dot3.multiply(a70.multiply(2).subtract(1))).
                                add(dot4.multiply(d0));
            final T coeffDot1 = time.getField().getZero();
            final T coeffDot2 = a72.
                                subtract(dot2.multiply(a72)).
                                add(dot3.multiply(a72.multiply(2))).
                                add(dot4.multiply(d2));
            final T coeffDot3 = a73.
                                subtract(dot2.multiply(a73)).
                                add(dot3.multiply(a73.multiply(2))).
                                add(dot4.multiply(d3));
            final T coeffDot4 = a74.
                                subtract(dot2.multiply(a74)).
                                add(dot3.multiply(a74.multiply(2))).
                                add(dot4.multiply(d4));
            final T coeffDot5 = a75.
                                subtract(dot2.multiply(a75)).
                                add(dot3.multiply(a75.multiply(2))).
                                add(dot4.multiply(d5));
            final T coeffDot6 = dot4.multiply(d6).subtract(dot3);
            interpolatedState       = previousStateLinearCombination(coeff0, coeff1, coeff2, coeff3,
                                                                     coeff4, coeff5, coeff6);
            interpolatedDerivatives = derivativeLinearCombination(coeffDot0, coeffDot1, coeffDot2, coeffDot3,
                                                                  coeffDot4, coeffDot5, coeffDot6);
        } else {
            final T f1        = oneMinusThetaH.negate();
            final T f2        = oneMinusThetaH.multiply(theta);
            final T f3        = f2.multiply(theta);
            final T f4        = f3.multiply(eta);
            final T coeff0    = f1.multiply(a70).
                                subtract(f2.multiply(a70.subtract(1))).
                                add(f3.multiply(a70.multiply(2).subtract(1))).
                                add(f4.multiply(d0));
            final T coeff1    = time.getField().getZero();
            final T coeff2    = f1.multiply(a72).
                                subtract(f2.multiply(a72)).
                                add(f3.multiply(a72.multiply(2))).
                                add(f4.multiply(d2));
            final T coeff3    = f1.multiply(a73).
                                subtract(f2.multiply(a73)).
                                add(f3.multiply(a73.multiply(2))).
                                add(f4.multiply(d3));
            final T coeff4    = f1.multiply(a74).
                                subtract(f2.multiply(a74)).
                                add(f3.multiply(a74.multiply(2))).
                                add(f4.multiply(d4));
            final T coeff5    = f1.multiply(a75).
                                subtract(f2.multiply(a75)).
                                add(f3.multiply(a75.multiply(2))).
                                add(f4.multiply(d5));
            final T coeff6    = f4.multiply(d6).subtract(f3);
            final T coeffDot0 = a70.
                                subtract(dot2.multiply(a70.subtract(1))).
                                add(dot3.multiply(a70.multiply(2).subtract(1))).
                                add(dot4.multiply(d0));
            final T coeffDot1 = time.getField().getZero();
            final T coeffDot2 = a72.
                                subtract(dot2.multiply(a72)).
                                add(dot3.multiply(a72.multiply(2))).
                                add(dot4.multiply(d2));
            final T coeffDot3 = a73.
                                subtract(dot2.multiply(a73)).
                                add(dot3.multiply(a73.multiply(2))).
                                add(dot4.multiply(d3));
            final T coeffDot4 = a74.
                                subtract(dot2.multiply(a74)).
                                add(dot3.multiply(a74.multiply(2))).
                                add(dot4.multiply(d4));
            final T coeffDot5 = a75.
                                subtract(dot2.multiply(a75)).
                                add(dot3.multiply(a75.multiply(2))).
                                add(dot4.multiply(d5));
            final T coeffDot6 = dot4.multiply(d6).subtract(dot3);
            interpolatedState       = currentStateLinearCombination(coeff0, coeff1, coeff2, coeff3,
                                                                    coeff4, coeff5, coeff6);
            interpolatedDerivatives = derivativeLinearCombination(coeffDot0, coeffDot1, coeffDot2, coeffDot3,
                                                                  coeffDot4, coeffDot5, coeffDot6);
        }
        return mapper.mapStateAndDerivative(time, interpolatedState, interpolatedDerivatives);

    }

}
