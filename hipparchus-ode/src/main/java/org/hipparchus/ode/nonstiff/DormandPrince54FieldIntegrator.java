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

package org.hipparchus.ode.nonstiff;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.ode.FieldEquationsMapper;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.ode.nonstiff.interpolators.DormandPrince54FieldStateInterpolator;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;


/**
 * This class implements the 5(4) Dormand-Prince integrator for Ordinary
 * Differential Equations.

 * <p>This integrator is an embedded Runge-Kutta integrator
 * of order 5(4) used in local extrapolation mode (i.e. the solution
 * is computed using the high order formula) with stepsize control
 * (and automatic step initialization) and continuous output. This
 * method uses 7 functions evaluations per step. However, since this
 * is an <i>fsal</i>, the last evaluation of one step is the same as
 * the first evaluation of the next step and hence can be avoided. So
 * the cost is really 6 functions evaluations per step.</p>
 *
 * <p>This method has been published (whithout the continuous output
 * that was added by Shampine in 1986) in the following article :</p>
 * <pre>
 *  A family of embedded Runge-Kutta formulae
 *  J. R. Dormand and P. J. Prince
 *  Journal of Computational and Applied Mathematics
 *  volume 6, no 1, 1980, pp. 19-26
 * </pre>
 *
 * @param <T> the type of the field elements
 */

public class DormandPrince54FieldIntegrator<T extends CalculusFieldElement<T>>
    extends EmbeddedRungeKuttaFieldIntegrator<T> {

    /** Name of integration scheme. */
    public static final String METHOD_NAME = DormandPrince54Integrator.METHOD_NAME;

    /** Simple constructor.
     * Build a fifth order Dormand-Prince integrator with the given step bounds
     * @param field field to which the time and state vector elements belong
     * @param minStep minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param maxStep maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param scalAbsoluteTolerance allowed absolute error
     * @param scalRelativeTolerance allowed relative error
     */
    public DormandPrince54FieldIntegrator(final Field<T> field,
                                          final double minStep, final double maxStep,
                                          final double scalAbsoluteTolerance,
                                          final double scalRelativeTolerance) {
        super(field, METHOD_NAME, 6,
              minStep, maxStep, scalAbsoluteTolerance, scalRelativeTolerance);
    }

    /** Simple constructor.
     * Build a fifth order Dormand-Prince integrator with the given step bounds
     * @param field field to which the time and state vector elements belong
     * @param minStep minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param maxStep maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param vecAbsoluteTolerance allowed absolute error
     * @param vecRelativeTolerance allowed relative error
     */
    public DormandPrince54FieldIntegrator(final Field<T> field,
                                          final double minStep, final double maxStep,
                                          final double[] vecAbsoluteTolerance,
                                          final double[] vecRelativeTolerance) {
        super(field, DormandPrince54Integrator.METHOD_NAME, 6,
              minStep, maxStep, vecAbsoluteTolerance, vecRelativeTolerance);
    }

    /** {@inheritDoc} */
    @Override
    public T[] getC() {
        final T[] c = MathArrays.buildArray(getField(), 6);
        c[0] = FieldExplicitRungeKuttaIntegrator.fraction(getField(), 1,  5);
        c[1] = FieldExplicitRungeKuttaIntegrator.fraction(getField(), 3, 10);
        c[2] = FieldExplicitRungeKuttaIntegrator.fraction(getField(),4,  5);
        c[3] = FieldExplicitRungeKuttaIntegrator.fraction(getField(),8,  9);
        c[4] = getField().getOne();
        c[5] = getField().getOne();
        return c;
    }

    /** {@inheritDoc} */
    @Override
    public T[][] getA() {
        final T[][] a = MathArrays.buildArray(getField(), 6, -1);
        for (int i = 0; i < a.length; ++i) {
            a[i] = MathArrays.buildArray(getField(), i + 1);
        }
        a[0][0] = FieldExplicitRungeKuttaIntegrator.fraction(getField(),      1,     5);
        a[1][0] = FieldExplicitRungeKuttaIntegrator.fraction(getField(),      3,    40);
        a[1][1] = FieldExplicitRungeKuttaIntegrator.fraction(getField(),      9,    40);
        a[2][0] = FieldExplicitRungeKuttaIntegrator.fraction(getField(),     44,    45);
        a[2][1] = FieldExplicitRungeKuttaIntegrator.fraction(getField(),    -56,    15);
        a[2][2] = FieldExplicitRungeKuttaIntegrator.fraction(getField(),     32,     9);
        a[3][0] = FieldExplicitRungeKuttaIntegrator.fraction(getField(),  19372,  6561);
        a[3][1] = FieldExplicitRungeKuttaIntegrator.fraction(getField(), -25360,  2187);
        a[3][2] = FieldExplicitRungeKuttaIntegrator.fraction(getField(),  64448,  6561);
        a[3][3] = FieldExplicitRungeKuttaIntegrator.fraction(getField(),   -212,   729);
        a[4][0] = FieldExplicitRungeKuttaIntegrator.fraction(getField(),   9017,  3168);
        a[4][1] = FieldExplicitRungeKuttaIntegrator.fraction(getField(),   -355,    33);
        a[4][2] = FieldExplicitRungeKuttaIntegrator.fraction(getField(),  46732,  5247);
        a[4][3] = FieldExplicitRungeKuttaIntegrator.fraction(getField(),     49,   176);
        a[4][4] = FieldExplicitRungeKuttaIntegrator.fraction(getField(),  -5103, 18656);
        a[5][0] = FieldExplicitRungeKuttaIntegrator.fraction(getField(),     35,   384);
        a[5][1] = getField().getZero();
        a[5][2] = FieldExplicitRungeKuttaIntegrator.fraction(getField(),    500,  1113);
        a[5][3] = FieldExplicitRungeKuttaIntegrator.fraction(getField(),    125,   192);
        a[5][4] = FieldExplicitRungeKuttaIntegrator.fraction(getField(),  -2187,  6784);
        a[5][5] = FieldExplicitRungeKuttaIntegrator.fraction(getField(),     11,    84);
        return a;
    }

    /** {@inheritDoc} */
    @Override
    public T[] getB() {
        final T[] b = MathArrays.buildArray(getField(), 7);
        b[0] = FieldExplicitRungeKuttaIntegrator.fraction(getField(),    35,   384);
        b[1] = getField().getZero();
        b[2] = FieldExplicitRungeKuttaIntegrator.fraction(getField(),   500, 1113);
        b[3] = FieldExplicitRungeKuttaIntegrator.fraction(getField(),   125,  192);
        b[4] = FieldExplicitRungeKuttaIntegrator.fraction(getField(), -2187, 6784);
        b[5] = FieldExplicitRungeKuttaIntegrator.fraction(getField(),    11,   84);
        b[6] = getField().getZero();
        return b;
    }

    /** {@inheritDoc} */
    @Override
    protected DormandPrince54FieldStateInterpolator<T>
        createInterpolator(final boolean forward, T[][] yDotK,
                           final FieldODEStateAndDerivative<T> globalPreviousState,
                           final FieldODEStateAndDerivative<T> globalCurrentState, final FieldEquationsMapper<T> mapper) {
        return new DormandPrince54FieldStateInterpolator<>(getField(), forward, yDotK,
                                                           globalPreviousState, globalCurrentState,
                                                           globalPreviousState, globalCurrentState,
                                                           mapper);
    }

    /** {@inheritDoc} */
    @Override
    public int getOrder() {
        return 5;
    }

    /** {@inheritDoc} */
    @Override
    protected double estimateError(final T[][] yDotK, final T[] y0, final T[] y1, final T h) {

        final StepsizeHelper helper = getStepSizeHelper();
        double error = 0;

        for (int j = 0; j < helper.getMainSetDimension(); ++j) {
            final double errSum = DormandPrince54Integrator.E1 * yDotK[0][j].getReal() +  DormandPrince54Integrator.E3 * yDotK[2][j].getReal() +
                                  DormandPrince54Integrator.E4 * yDotK[3][j].getReal() +  DormandPrince54Integrator.E5 * yDotK[4][j].getReal() +
                                  DormandPrince54Integrator.E6 * yDotK[5][j].getReal() +  DormandPrince54Integrator.E7 * yDotK[6][j].getReal();
            final double tol = helper.getTolerance(j, FastMath.max(FastMath.abs(y0[j].getReal()), FastMath.abs(y1[j].getReal())));
            final double ratio  = h.getReal() * errSum / tol;
            error += ratio * ratio;
        }

        return FastMath.sqrt(error / helper.getMainSetDimension());

    }

}
