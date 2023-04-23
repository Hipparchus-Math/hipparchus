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


import java.lang.reflect.InvocationTargetException;

import org.hipparchus.Field;
import org.hipparchus.CalculusFieldElement;
import org.hipparchus.ode.EquationsMapper;
import org.hipparchus.ode.FieldEquationsMapper;
import org.hipparchus.ode.FieldExpandableODE;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.ode.FieldOrdinaryDifferentialEquation;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.sampling.AbstractFieldODEStateInterpolator;
import org.hipparchus.ode.sampling.FieldODEStateInterpolator;
import org.hipparchus.util.MathArrays;
import org.junit.Assert;

public abstract class RungeKuttaFieldStateInterpolatorAbstractTest extends FieldODEStateInterpolatorAbstractTest {

    protected abstract <T extends CalculusFieldElement<T>> RungeKuttaFieldStateInterpolator<T>
        createInterpolator(Field<T> field, boolean forward, T[][] yDotK,
                           FieldODEStateAndDerivative<T> globalPreviousState,
                           FieldODEStateAndDerivative<T> globalCurrentState,
                           FieldODEStateAndDerivative<T> softPreviousState,
                           FieldODEStateAndDerivative<T> softCurrentState,
                           FieldEquationsMapper<T> mapper);

    protected abstract <T extends CalculusFieldElement<T>> FieldButcherArrayProvider<T>
        createButcherArrayProvider(final Field<T> field);

    protected <T extends CalculusFieldElement<T>>
    RungeKuttaFieldStateInterpolator<T> setUpInterpolator(final Field<T> field,
                                                          final ReferenceFieldODE<T> eqn,
                                                          final double t0, final double[] y0,
                                                          final double t1) {

        // get the Butcher arrays from the field integrator
        FieldButcherArrayProvider<T> provider = createButcherArrayProvider(field);
        T[][] a = provider.getA();
        T[]   b = provider.getB();
        T[]   c = provider.getC();

        // store initial state
        T     t          = field.getZero().add(t0);
        T[]   fieldY     = MathArrays.buildArray(field, eqn.getDimension());
        T[][] fieldYDotK = MathArrays.buildArray(field, b.length, -1);
        for (int i = 0; i < y0.length; ++i) {
            fieldY[i] = field.getZero().add(y0[i]);
        }
        fieldYDotK[0] = eqn.computeDerivatives(t, fieldY);
        FieldODEStateAndDerivative<T> s0 = new FieldODEStateAndDerivative<T>(t, fieldY, fieldYDotK[0]);

        // perform one integration step, in order to get consistent derivatives
        T h = field.getZero().add(t1 - t0);
        for (int k = 0; k < a.length; ++k) {
            for (int i = 0; i < y0.length; ++i) {
                fieldY[i] = field.getZero().add(y0[i]);
                for (int s = 0; s <= k; ++s) {
                    fieldY[i] = fieldY[i].add(h.multiply(a[k][s].multiply(fieldYDotK[s][i])));
                }
            }
            fieldYDotK[k + 1] = eqn.computeDerivatives(h.multiply(c[k]).add(t0), fieldY);
        }

        // store state at step end
        t = field.getZero().add(t1);
        for (int i = 0; i < y0.length; ++i) {
            fieldY[i] = field.getZero().add(y0[i]);
            for (int s = 0; s < b.length; ++s) {
                fieldY[i] = fieldY[i].add(h.multiply(b[s].multiply(fieldYDotK[s][i])));
            }
        }
        FieldODEStateAndDerivative<T> s1 = new FieldODEStateAndDerivative<T>(t, fieldY,
                                                                             eqn.computeDerivatives(t, fieldY));

        return createInterpolator(field, t1 > t0, fieldYDotK, s0, s1, s0, s1,
                                  new FieldExpandableODE<T>(eqn).getMapper());

    }

    protected <T extends CalculusFieldElement<T>>
    RungeKuttaStateInterpolator convertInterpolator(final FieldODEStateInterpolator<T> fieldInterpolator,
                                                    final FieldOrdinaryDifferentialEquation<T> eqn) {

        RungeKuttaFieldStateInterpolator<T> rkFieldInterpolator =
                        (RungeKuttaFieldStateInterpolator<T>) fieldInterpolator;

        RungeKuttaStateInterpolator regularInterpolator = null;
        try {

            String interpolatorName = rkFieldInterpolator.getClass().getName();
            String integratorName = interpolatorName.replaceAll("Field", "");
            @SuppressWarnings("unchecked")
            Class<RungeKuttaStateInterpolator> clz = (Class<RungeKuttaStateInterpolator>) Class.forName(integratorName);

            java.lang.reflect.Field fYD = RungeKuttaFieldStateInterpolator.class.getDeclaredField("yDotK");
            fYD.setAccessible(true);
            @SuppressWarnings("unchecked")
            final double[][] yDotK = convertArray((T[][]) fYD.get(rkFieldInterpolator));

            java.lang.reflect.Field fMapper = AbstractFieldODEStateInterpolator.class.getDeclaredField("mapper");
            fMapper.setAccessible(true);
            @SuppressWarnings("unchecked")
            EquationsMapper regularMapper = convertMapper((FieldEquationsMapper<T>) fMapper.get(rkFieldInterpolator));

            java.lang.reflect.Constructor<RungeKuttaStateInterpolator> regularInterpolatorConstructor =
                            clz.getDeclaredConstructor(Boolean.TYPE,
                                                       double[][].class,
                                                       ODEStateAndDerivative.class,
                                                       ODEStateAndDerivative.class,
                                                       ODEStateAndDerivative.class,
                                                       ODEStateAndDerivative.class,
                                                       EquationsMapper.class);
            return regularInterpolatorConstructor.newInstance(rkFieldInterpolator.isForward(),
                                                              yDotK,
                                                              convertODEStateAndDerivative(rkFieldInterpolator.getGlobalPreviousState()),
                                                              convertODEStateAndDerivative(rkFieldInterpolator.getGlobalCurrentState()),
                                                              convertODEStateAndDerivative(rkFieldInterpolator.getPreviousState()),
                                                              convertODEStateAndDerivative(rkFieldInterpolator.getCurrentState()),
                                                              regularMapper);

        } catch (ClassNotFoundException | InstantiationException   | IllegalAccessException    |
                 NoSuchFieldException   | IllegalArgumentException | InvocationTargetException |
                 NoSuchMethodException  | SecurityException e) {
            Assert.fail(e.getLocalizedMessage());
        }

        return regularInterpolator;

    }

}
