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
import org.hipparchus.Field;
import org.hipparchus.ode.FieldExpandableODE;
import org.hipparchus.ode.FieldODEIntegrator;
import org.hipparchus.ode.FieldOrdinaryDifferentialEquation;
import org.hipparchus.util.MathArrays;

/**
 * This interface implements the part of Runge-Kutta
 * Field integrators for Ordinary Differential Equations
 * common to fixed- and adaptive steps.
 *
 * <p>These methods are explicit Runge-Kutta methods, their Butcher
 * arrays are as follows :</p>
 * <pre>
 *    0  |
 *   c2  | a21
 *   c3  | a31  a32
 *   ... |        ...
 *   cs  | as1  as2  ...  ass-1
 *       |--------------------------
 *       |  b1   b2  ...   bs-1  bs
 * </pre>
 *
 * @see FieldButcherArrayProvider
 * @see FixedStepRungeKuttaFieldIntegrator
 * @see EmbeddedRungeKuttaFieldIntegrator
 * @param <T> the type of the field elements
 * @since 3.1
 */

public interface FieldExplicitRungeKuttaIntegrator<T extends CalculusFieldElement<T>>
    extends FieldButcherArrayProvider<T>, FieldODEIntegrator<T> {

    /** Get the time steps from Butcher array (without the first zero). Real version (non-Field).
     * @return time steps from Butcher array (without the first zero).
     */
    default double[] getRealC() {
        final T[] c = getC();
        final double[] cReal = new double[c.length];
        for (int i = 0; i < c.length; i++) {
            cReal[i] = c[i].getReal();
        }
        return cReal;
    }

    /** Get the internal weights from Butcher array (without the first empty row). Real version (non-Field).
     * @return internal weights from Butcher array (without the first empty row)
     */
    default double[][] getRealA() {
        final T[][] a = getA();
        final double[][] aReal = new double[a.length][];
        for (int i = 0; i < a.length; i++) {
            aReal[i] = new double[a[i].length];
            for (int j = 0; j < aReal[i].length; j++) {
                aReal[i][j] = a[i][j].getReal();
            }
        }
        return aReal;
    }

    /** Get the external weights for the high order method from Butcher array. Real version (non-Field).
     * @return external weights for the high order method from Butcher array
     */
    default double[] getRealB() {
        final T[] b = getB();
        final double[] bReal = new double[b.length];
        for (int i = 0; i < b.length; i++) {
            bReal[i] = b[i].getReal();
        }
        return bReal;
    }

    /**
     * Getter for the flag between real or Field coefficients in the Butcher array.
     *
     * @return flag
     */
    boolean isUsingFieldCoefficients();

    /**
     * Getter for the number of stages corresponding to the Butcher array.
     *
     * @return number of stages
     */
    default int getNumberOfStages() {
        return getB().length;
    }

    /** Fast computation of a single step of ODE integration.
     * <p>This method is intended for the limited use case of
     * very fast computation of only one step without using any of the
     * rich features of general integrators that may take some time
     * to set up (i.e. no step handlers, no events handlers, no additional
     * states, no interpolators, no error control, no evaluations count,
     * no sanity checks ...). It handles the strict minimum of computation,
     * so it can be embedded in outer loops.</p>
     * <p>
     * This method is <em>not</em> used at all by the {@link #integrate(FieldExpandableODE,
     * org.hipparchus.ode.FieldODEState, CalculusFieldElement)} method. It also completely ignores the step set at
     * construction time, and uses only a single step to go from {@code t0} to {@code t}.
     * </p>
     * <p>
     * As this method does not use any of the state-dependent features of the integrator,
     * it should be reasonably thread-safe <em>if and only if</em> the provided differential
     * equations are themselves thread-safe.
     * </p>
     * @param equations differential equations to integrate
     * @param t0 initial time
     * @param y0 initial value of the state vector at t0
     * @param t target time for the integration
     * (can be set to a value smaller than {@code t0} for backward integration)
     * @return state vector at {@code t}
     */
    default T[] singleStep(final FieldOrdinaryDifferentialEquation<T> equations, final T t0, final T[] y0, final T t) {

        // create some internal working arrays
        final int stages  = getNumberOfStages();
        final T[][] yDotK = MathArrays.buildArray(t0.getField(), stages, -1);

        // first stage
        final T h = t.subtract(t0);
        final FieldExpandableODE<T> fieldExpandableODE = new FieldExpandableODE<>(equations);
        yDotK[0] = fieldExpandableODE.computeDerivatives(t0, y0);

        if (isUsingFieldCoefficients()) {
            applyInternalButcherWeights(fieldExpandableODE, t0, y0, h, getA(), getC(), yDotK);
            return applyExternalButcherWeights(y0, yDotK, h, getB());
        } else {
            applyInternalButcherWeights(fieldExpandableODE, t0, y0, h, getRealA(), getRealC(), yDotK);
            return applyExternalButcherWeights(y0, yDotK, h, getRealB());
        }
    }

    /**
     * Create a fraction from integers.
     *
     * @param <T> the type of the field elements
     * @param field field to which elements belong
     * @param p numerator
     * @param q denominator
     * @return p/q computed in the instance field
     */
    static <T extends CalculusFieldElement<T>> T fraction(final Field<T> field, final int p, final int q) {
        final T zero = field.getZero();
        return zero.newInstance(p).divide(zero.newInstance(q));
    }

    /**
     * Create a fraction from doubles.
     * @param <T> the type of the field elements
     * @param field field to which elements belong
     * @param p numerator
     * @param q denominator
     * @return p/q computed in the instance field
     */
    static <T extends CalculusFieldElement<T>> T fraction(final Field<T> field, final double p, final double q) {
        final T zero = field.getZero();
        return zero.newInstance(p).divide(zero.newInstance(q));
    }

    /**
     * Apply internal weights of Butcher array, with corresponding times.
     * @param <T> the type of the field elements
     * @param equations differential equations to integrate
     * @param t0        initial time
     * @param y0        initial value of the state vector at t0
     * @param h         step size
     * @param a         internal weights of Butcher array
     * @param c         times of Butcher array
     * @param yDotK     array where to store result
     */
    static <T extends CalculusFieldElement<T>> void applyInternalButcherWeights(final FieldExpandableODE<T> equations,
                                                                                final T t0, final T[] y0, final T h,
                                                                                final T[][] a, final T[] c,
                                                                                final T[][] yDotK) {
        // create some internal working arrays
        final int stages = c.length + 1;
        final T[] yTmp = y0.clone();

        for (int k = 1; k < stages; ++k) {

            for (int j = 0; j < y0.length; ++j) {
                T sum = yDotK[0][j].multiply(a[k - 1][0]);
                for (int l = 1; l < k; ++l) {
                    sum = sum.add(yDotK[l][j].multiply(a[k - 1][l]));
                }
                yTmp[j] = y0[j].add(h.multiply(sum));
            }

            yDotK[k] = equations.computeDerivatives(t0.add(h.multiply(c[k - 1])), yTmp);
        }
    }

    /** Apply internal weights of Butcher array, with corresponding times. Version with real Butcher array (non-Field).
     * @param <T> the type of the field elements
     * @param equations differential equations to integrate
     * @param t0 initial time
     * @param y0 initial value of the state vector at t0
     * @param h step size
     * @param a internal weights of Butcher array
     * @param c times of Butcher array
     * @param yDotK array where to store result
     */
    static <T extends CalculusFieldElement<T>> void applyInternalButcherWeights(final FieldExpandableODE<T> equations,
                                                                                final T t0, final T[] y0, final T h,
                                                                                final double[][] a, final double[] c,
                                                                                final T[][] yDotK) {
        // create some internal working arrays
        final int stages = c.length + 1;
        final T[] yTmp = y0.clone();

        for (int k = 1; k < stages; ++k) {

            for (int j = 0; j < y0.length; ++j) {
                T sum = yDotK[0][j].multiply(a[k - 1][0]);
                for (int l = 1; l < k; ++l) {
                    sum = sum.add(yDotK[l][j].multiply(a[k - 1][l]));
                }
                yTmp[j] = y0[j].add(h.multiply(sum));
            }

            yDotK[k] = equations.computeDerivatives(t0.add(h.multiply(c[k - 1])), yTmp);
        }
    }

    /** Apply external weights of Butcher array, assuming internal ones have been applied.
     * @param <T> the type of the field elements
     * @param yDotK output of stages
     * @param y0 initial value of the state vector at t0
     * @param h step size
     * @param b external weights of Butcher array
     * @return state vector
     */
    static <T extends CalculusFieldElement<T>> T[] applyExternalButcherWeights(final T[] y0, final T[][] yDotK,
                                                                               final T h, final T[] b) {
        final T[] y = y0.clone();
        final int stages = b.length;
        for (int j = 0; j < y0.length; ++j) {
            T sum = yDotK[0][j].multiply(b[0]);
            for (int l = 1; l < stages; ++l) {
                sum = sum.add(yDotK[l][j].multiply(b[l]));
            }
            y[j] = y[j].add(h.multiply(sum));
        }
        return y;
    }

    /** Apply external weights of Butcher array, assuming internal ones have been applied. Version with real Butcher
     * array (non-Field version).
     * @param <T> the type of the field elements
     * @param yDotK output of stages
     * @param y0 initial value of the state vector at t0
     * @param h step size
     * @param b external weights of Butcher array
     * @return state vector
     */
    static <T extends CalculusFieldElement<T>> T[] applyExternalButcherWeights(final T[] y0, final T[][] yDotK,
                                                                               final T h, final double[] b) {
        final T[] y = y0.clone();
        final int stages = b.length;
        for (int j = 0; j < y0.length; ++j) {
            T sum = yDotK[0][j].multiply(b[0]);
            for (int l = 1; l < stages; ++l) {
                sum = sum.add(yDotK[l][j].multiply(b[l]));
            }
            y[j] = y[j].add(h.multiply(sum));
        }
        return y;
    }

}
