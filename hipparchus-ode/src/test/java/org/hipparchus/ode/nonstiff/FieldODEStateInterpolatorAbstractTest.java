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
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.ode.FieldOrdinaryDifferentialEquation;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.sampling.AbstractFieldODEStateInterpolator;
import org.hipparchus.ode.sampling.FieldODEStateInterpolator;
import org.hipparchus.ode.sampling.ODEStateInterpolator;
import org.hipparchus.util.Binary64Field;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;
import org.junit.Assert;
import org.junit.Test;

public abstract class FieldODEStateInterpolatorAbstractTest {

    @Test
    public abstract void interpolationAtBounds();

    protected <T extends CalculusFieldElement<T>> void doInterpolationAtBounds(final Field<T> field, double epsilon) {

        FieldODEStateInterpolator<T> interpolator = setUpInterpolator(field,
                                                                      new SinCos<T>(field),
                                                                      0.0, new double[] { 0.0, 1.0 }, 0.125);

        Assert.assertEquals(0.0, interpolator.getPreviousState().getTime().getReal(), 1.0e-15);
        for (int i = 0; i < 2; ++i) {
            Assert.assertEquals(interpolator.getPreviousState().getPrimaryState()[i].getReal(),
                                interpolator.getInterpolatedState(interpolator.getPreviousState().getTime()).getPrimaryState()[i].getReal(),
                                epsilon);
        }
        Assert.assertEquals(0.125, interpolator.getCurrentState().getTime().getReal(), 1.0e-15);
        for (int i = 0; i < 2; ++i) {
            Assert.assertEquals(interpolator.getCurrentState().getPrimaryState()[i].getReal(),
                                interpolator.getInterpolatedState(interpolator.getCurrentState().getTime()).getPrimaryState()[i].getReal(),
                                epsilon);
        }

        Assert.assertEquals(false, interpolator.isPreviousStateInterpolated());
        Assert.assertEquals(false, interpolator.isCurrentStateInterpolated());
    }

    @Test
    public abstract void interpolationInside();

    protected <T extends CalculusFieldElement<T>> void doInterpolationInside(final Field<T> field,
                                                                         double epsilonSin, double epsilonCos) {

        ReferenceFieldODE<T> sinCos =  new SinCos<T>(field);
        FieldODEStateInterpolator<T> interpolator = setUpInterpolator(field, sinCos,
                                                                      0.0, new double[] { 0.0, 1.0 }, 0.0125);

        int n = 100;
        double maxErrorSin = 0;
        double maxErrorCos = 0;
        for (int i = 0; i <= n; ++i) {
            T t =     interpolator.getPreviousState().getTime().multiply(n - i).
                  add(interpolator.getCurrentState().getTime().multiply(i)).
                  divide(n);
            FieldODEStateAndDerivative<T> state = interpolator.getInterpolatedState(t);
            T[] ref = sinCos.theoreticalState(t);
            maxErrorSin = FastMath.max(maxErrorSin, state.getPrimaryState()[0].subtract(ref[0]).norm());
            maxErrorCos = FastMath.max(maxErrorCos, state.getPrimaryState()[1].subtract(ref[1]).norm());
        }
        Assert.assertEquals(0.0, maxErrorSin, epsilonSin);
        Assert.assertEquals(0.0, maxErrorCos, epsilonCos);

        Assert.assertEquals(false, interpolator.isPreviousStateInterpolated());
        Assert.assertEquals(false, interpolator.isCurrentStateInterpolated());
    }

    @Test
    public void restrictPrevious() {
        doRestrictPrevious(Binary64Field.getInstance(), 1e-15, 1e-15);
    }

    protected <T extends CalculusFieldElement<T>> void doRestrictPrevious(
            Field<T> field,
            double epsilon,
            double epsilonDot) {

        AbstractFieldODEStateInterpolator<T> original = setUpInterpolator(
                field, new SinCos<>(field), 0.0, new double[]{0.0, 1.0}, 0.125);

        Assert.assertEquals(false, original.isPreviousStateInterpolated());
        Assert.assertEquals(false, original.isCurrentStateInterpolated());

        AbstractFieldODEStateInterpolator<T> restricted = original.restrictStep(
                original.getInterpolatedState(field.getZero().add(1.0 / 32)),
                original.getCurrentState());

        Assert.assertSame(original.getPreviousState(),       original.getGlobalPreviousState());
        Assert.assertSame(original.getCurrentState(),        original.getGlobalCurrentState());
        Assert.assertSame(original.getGlobalPreviousState(), restricted.getGlobalPreviousState());
        Assert.assertSame(original.getGlobalCurrentState(),  restricted.getGlobalCurrentState());
        Assert.assertNotSame(restricted.getPreviousState(),  restricted.getGlobalPreviousState());
        Assert.assertSame(restricted.getCurrentState(),      restricted.getGlobalCurrentState());
        Assert.assertEquals(1.0 / 32, restricted.getPreviousState().getTime().getReal(), 1.0e-15);
        Assert.assertEquals(true, restricted.isPreviousStateInterpolated());
        Assert.assertEquals(false, restricted.isCurrentStateInterpolated());

        checkRestricted(original, restricted, epsilon, epsilonDot);

    }

    @Test
    public void restrictCurrent() {
        doRestrictCurrent(Binary64Field.getInstance(), 1e-15, 1e-15);
    }

    protected <T extends CalculusFieldElement<T>> void doRestrictCurrent(Field<T> field,
                                                                     double epsilon,
                                                                     double epsilonDot) {

        AbstractFieldODEStateInterpolator<T> original = setUpInterpolator(
                field, new SinCos<>(field), 0.0, new double[]{0.0, 1.0}, 0.125);

        Assert.assertEquals(false, original.isPreviousStateInterpolated());
        Assert.assertEquals(false, original.isCurrentStateInterpolated());

        AbstractFieldODEStateInterpolator<T> restricted = original.restrictStep(
                original.getPreviousState(),
                original.getInterpolatedState(field.getZero().add(3.0 / 32)));

        Assert.assertSame(original.getPreviousState(),       original.getGlobalPreviousState());
        Assert.assertSame(original.getCurrentState(),        original.getGlobalCurrentState());
        Assert.assertSame(original.getGlobalPreviousState(), restricted.getGlobalPreviousState());
        Assert.assertSame(original.getGlobalCurrentState(),  restricted.getGlobalCurrentState());
        Assert.assertSame(restricted.getPreviousState(),     restricted.getGlobalPreviousState());
        Assert.assertNotSame(restricted.getCurrentState(),   restricted.getGlobalCurrentState());
        Assert.assertEquals(3.0 / 32, restricted.getCurrentState().getTime().getReal(), 1.0e-15);
        Assert.assertEquals(false, restricted.isPreviousStateInterpolated());
        Assert.assertEquals(true, restricted.isCurrentStateInterpolated());

        checkRestricted(original, restricted, epsilon, epsilonDot);

    }

    @Test
    public void restrictBothEnds() {
        doRestrictBothEnds(Binary64Field.getInstance(), 1e-15, 1e-15);
    }

    protected <T extends CalculusFieldElement<T>> void doRestrictBothEnds(Field<T> field,
                                                                      double epsilon,
                                                                      double epsilonDot) {

        AbstractFieldODEStateInterpolator<T> original = setUpInterpolator(
                field, new SinCos<>(field), 0.0, new double[]{0.0, 1.0}, 0.125);

        Assert.assertEquals(false, original.isPreviousStateInterpolated());
        Assert.assertEquals(false, original.isCurrentStateInterpolated());

        AbstractFieldODEStateInterpolator<T> restricted = original.restrictStep(
                original.getInterpolatedState(field.getZero().add(1.0 / 32)),
                original.getInterpolatedState(field.getZero().add(3.0 / 32)));

        Assert.assertSame(original.getPreviousState(),       original.getGlobalPreviousState());
        Assert.assertSame(original.getCurrentState(),        original.getGlobalCurrentState());
        Assert.assertSame(original.getGlobalPreviousState(), restricted.getGlobalPreviousState());
        Assert.assertSame(original.getGlobalCurrentState(),  restricted.getGlobalCurrentState());
        Assert.assertNotSame(restricted.getPreviousState(),  restricted.getGlobalPreviousState());
        Assert.assertNotSame(restricted.getCurrentState(),   restricted.getGlobalCurrentState());
        Assert.assertEquals(1.0 / 32, restricted.getPreviousState().getTime().getReal(), 1.0e-15);
        Assert.assertEquals(3.0 / 32, restricted.getCurrentState().getTime().getReal(), 1.0e-15);
        Assert.assertEquals(true, restricted.isPreviousStateInterpolated());
        Assert.assertEquals(true, restricted.isCurrentStateInterpolated());

        checkRestricted(original, restricted, epsilon, epsilonDot);

    }

    @Test
    public void degenerateInterpolation() {
        doDegenerateInterpolation(Binary64Field.getInstance());
    }

    protected <T extends CalculusFieldElement<T>> void doDegenerateInterpolation(Field<T> field) {
        AbstractFieldODEStateInterpolator<T> interpolator = setUpInterpolator(
                field, new SinCos<>(field), 0.0, new double[] { 0.0, 1.0 }, 0.0);
        FieldODEStateAndDerivative<T> interpolatedState = interpolator.getInterpolatedState(field.getZero());
        Assert.assertEquals(0.0, interpolatedState.getTime().getReal(), 0.0);
        Assert.assertEquals(0.0, interpolatedState.getPrimaryState()[0].getReal(), 0.0);
        Assert.assertEquals(1.0, interpolatedState.getPrimaryState()[1].getReal(), 0.0);
        Assert.assertEquals(1.0, interpolatedState.getPrimaryDerivative()[0].getReal(), 0.0);
        Assert.assertEquals(0.0, interpolatedState.getPrimaryDerivative()[1].getReal(), 0.0);
    }

    private <T extends CalculusFieldElement<T>> void checkRestricted(
            AbstractFieldODEStateInterpolator<T> original,
            AbstractFieldODEStateInterpolator<T> restricted,
            double epsilon,
            double epsilonDot) {

        for (T t = restricted.getPreviousState().getTime();
             t.getReal() <= restricted.getCurrentState().getTime().getReal();
             t = t.add(1.0 / 256)) {
            FieldODEStateAndDerivative<T> originalInterpolated   = original.getInterpolatedState(t);
            FieldODEStateAndDerivative<T> restrictedInterpolated = restricted.getInterpolatedState(t);
            Assert.assertEquals(t.getReal(), originalInterpolated.getTime().getReal(), 1.0e-15);
            Assert.assertEquals(t.getReal(), restrictedInterpolated.getTime().getReal(), 1.0e-15);
            Assert.assertEquals(originalInterpolated.getPrimaryState()[0].getReal(),
                                restrictedInterpolated.getPrimaryState()[0].getReal(),
                                epsilon);
            Assert.assertEquals(originalInterpolated.getPrimaryState()[1].getReal(),
                                restrictedInterpolated.getPrimaryState()[1].getReal(),
                                epsilon);
            Assert.assertEquals(originalInterpolated.getPrimaryDerivative()[0].getReal(),
                                restrictedInterpolated.getPrimaryDerivative()[0].getReal(),
                                epsilonDot);
            Assert.assertEquals(originalInterpolated.getPrimaryDerivative()[1].getReal(),
                                restrictedInterpolated.getPrimaryDerivative()[1].getReal(),
                                epsilonDot);
        }

    }

    @Test
    public abstract void nonFieldInterpolatorConsistency();

    protected <T extends CalculusFieldElement<T>> void doNonFieldInterpolatorConsistency(final Field<T> field,
                                                                                     double epsilonSin, double epsilonCos,
                                                                                     double epsilonSinDot, double epsilonCosDot) {

        ReferenceFieldODE<T> eqn = new SinCos<T>(field);
        FieldODEStateInterpolator<T> fieldInterpolator =
                        setUpInterpolator(field, eqn, 0.0, new double[] { 0.0, 1.0 }, 0.125);
        ODEStateInterpolator regularInterpolator = convertInterpolator(fieldInterpolator, eqn);

        int n = 100;
        double maxErrorSin    = 0;
        double maxErrorCos    = 0;
        double maxErrorSinDot = 0;
        double maxErrorCosDot = 0;
        for (int i = 0; i <= n; ++i) {

            T t =     fieldInterpolator.getPreviousState().getTime().multiply(n - i).
                  add(fieldInterpolator.getCurrentState().getTime().multiply(i)).
                  divide(n);

            FieldODEStateAndDerivative<T> state = fieldInterpolator.getInterpolatedState(t);
            T[] fieldY    = state.getPrimaryState();
            T[] fieldYDot = state.getPrimaryDerivative();

            ODEStateAndDerivative regularState = regularInterpolator.getInterpolatedState(t.getReal());
            double[] regularY     = regularState.getPrimaryState();
            double[] regularYDot  = regularState.getPrimaryDerivative();

            maxErrorSin    = FastMath.max(maxErrorSin,    fieldY[0].subtract(regularY[0]).norm());
            maxErrorCos    = FastMath.max(maxErrorCos,    fieldY[1].subtract(regularY[1]).norm());
            maxErrorSinDot = FastMath.max(maxErrorSinDot, fieldYDot[0].subtract(regularYDot[0]).norm());
            maxErrorCosDot = FastMath.max(maxErrorCosDot, fieldYDot[1].subtract(regularYDot[1]).norm());

        }
        Assert.assertEquals(0.0, maxErrorSin,    epsilonSin);
        Assert.assertEquals(0.0, maxErrorCos,    epsilonCos);
        Assert.assertEquals(0.0, maxErrorSinDot, epsilonSinDot);
        Assert.assertEquals(0.0, maxErrorCosDot, epsilonCosDot);

    }

    public interface ReferenceFieldODE<T extends CalculusFieldElement<T>> extends FieldOrdinaryDifferentialEquation<T> {
        T[] theoreticalState(T t);
    }

    protected abstract <T extends CalculusFieldElement<T>>
    AbstractFieldODEStateInterpolator<T> setUpInterpolator(final Field<T> field,
                                                           final ReferenceFieldODE<T> eqn,
                                                           final double t0,
                                                           final double[] y0,
                                                           final double t1);

    protected abstract <T extends CalculusFieldElement<T>>
    ODEStateInterpolator convertInterpolator(final FieldODEStateInterpolator<T> fieldInterpolator,
                                             final FieldOrdinaryDifferentialEquation<T> eqn);

    protected <T extends CalculusFieldElement<T>>
    ODEStateAndDerivative convertODEStateAndDerivative(final FieldODEStateAndDerivative<T> s) {
        final double[][] secondaryStates;
        final double[][] secondaryDerivatives;
        if (s.getNumberOfSecondaryStates() == 0) {
            secondaryStates      = null;
            secondaryDerivatives = null;
        } else {
            secondaryStates      = new double[s.getNumberOfSecondaryStates()][];
            secondaryDerivatives = new double[s.getNumberOfSecondaryStates()][];
            for (int i = 0; i < secondaryStates.length; ++i) {
                secondaryStates[i]      = convertArray(s.getSecondaryState(i));
                secondaryDerivatives[i] = convertArray(s.getSecondaryDerivative(i));
            }
        }
        return new ODEStateAndDerivative(s.getTime().getReal(),
                                         convertArray(s.getPrimaryState()),
                                         convertArray(s.getPrimaryDerivative()),
                                         secondaryStates,
                                         secondaryDerivatives);
    }

    protected <T extends CalculusFieldElement<T>> double[][] convertArray(final T[][] fieldArray) {
        if (fieldArray == null) {
            return null;
        }
        double[][] array = new double[fieldArray.length][];
        for (int i = 0; i < array.length; ++i) {
            array[i] = convertArray(fieldArray[i]);
        }
        return array;
    }

    protected <T extends CalculusFieldElement<T>> double[] convertArray(final T[] fieldArray) {
        if (fieldArray == null) {
            return null;
        }
        double[] array = new double[fieldArray.length];
        for (int i = 0; i < array.length; ++i) {
            array[i] = fieldArray[i].getReal();
        }
        return array;
    }

    protected <T extends CalculusFieldElement<T>>
    EquationsMapper convertMapper(final FieldEquationsMapper<T> fieldmapper)
        throws NoSuchMethodException, SecurityException, NoSuchFieldException,
               IllegalArgumentException, IllegalAccessException,
               InstantiationException, InvocationTargetException {
        java.lang.reflect.Field fStart = FieldEquationsMapper.class.getDeclaredField("start");
        fStart.setAccessible(true);
        int[] start = (int[]) fStart.get(fieldmapper);

        java.lang.reflect.Constructor<EquationsMapper> regularMapperConstructor =
                        EquationsMapper.class.getDeclaredConstructor(EquationsMapper.class,
                                                                     Integer.TYPE);
        regularMapperConstructor.setAccessible(true);
        EquationsMapper regularMapper = regularMapperConstructor.newInstance(null, start[1]);
        for (int k = 1; k < fieldmapper.getNumberOfEquations(); ++k) {
            regularMapper = regularMapperConstructor.newInstance(regularMapper,
                                                                 start[k + 1] - start[k]);
        }
        return regularMapper;
    }

    private static class SinCos<T extends CalculusFieldElement<T>> implements ReferenceFieldODE<T> {
        private final Field<T> field;
        protected SinCos(final Field<T> field) {
            this.field = field;
        }
        public int getDimension() {
            return 2;
        }
        public void init(final T t0, final T[] y0, final T finalTime) {
        }
        public T[] computeDerivatives(final T t, final T[] y) {
            T[] yDot = MathArrays.buildArray(field, 2);
            yDot[0] = y[1];
            yDot[1] = y[0].negate();
            return yDot;
        }
        public T[] theoreticalState(final T t) {
            T[] state = MathArrays.buildArray(field, 2);
            state[0] = t.sin();
            state[1] = t.cos();
            return state;
        }
    }

}
