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


import java.lang.reflect.Array;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.analysis.differentiation.DSFactory;
import org.hipparchus.analysis.differentiation.DerivativeStructure;
import org.hipparchus.analysis.differentiation.Gradient;
import org.hipparchus.analysis.solvers.BracketedRealFieldUnivariateSolver;
import org.hipparchus.analysis.solvers.FieldBracketingNthOrderBrentSolver;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.FieldExpandableODE;
import org.hipparchus.ode.FieldODEIntegrator;
import org.hipparchus.ode.FieldODEState;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.ode.FieldOrdinaryDifferentialEquation;
import org.hipparchus.ode.FieldSecondaryODE;
import org.hipparchus.ode.LocalizedODEFormats;
import org.hipparchus.ode.TestFieldProblem1;
import org.hipparchus.ode.TestFieldProblem2;
import org.hipparchus.ode.TestFieldProblem3;
import org.hipparchus.ode.TestFieldProblem4;
import org.hipparchus.ode.TestFieldProblem5;
import org.hipparchus.ode.TestFieldProblem6;
import org.hipparchus.ode.TestFieldProblemAbstract;
import org.hipparchus.ode.TestFieldProblemHandler;
import org.hipparchus.ode.events.Action;
import org.hipparchus.ode.events.FieldAdaptableInterval;
import org.hipparchus.ode.events.FieldODEEventDetector;
import org.hipparchus.ode.events.FieldODEEventHandler;
import org.hipparchus.ode.sampling.FieldODEStateInterpolator;
import org.hipparchus.ode.sampling.FieldODEStepHandler;
import org.hipparchus.ode.sampling.StepInterpolatorTestUtils;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;
import org.junit.Assert;
import org.junit.Test;

public abstract class RungeKuttaFieldIntegratorAbstractTest {

    protected abstract <T extends CalculusFieldElement<T>> RungeKuttaFieldIntegrator<T>
        createIntegrator(Field<T> field, T step);

    @Test
    public abstract void testNonFieldIntegratorConsistency();

    protected <T extends CalculusFieldElement<T>> void doTestNonFieldIntegratorConsistency(final Field<T> field) {
        try {

            // get the Butcher arrays from the field integrator
            RungeKuttaFieldIntegrator<T> fieldIntegrator = createIntegrator(field, field.getZero().add(1));
            T[][] fieldA = fieldIntegrator.getA();
            T[]   fieldB = fieldIntegrator.getB();
            T[]   fieldC = fieldIntegrator.getC();

            String fieldName   = fieldIntegrator.getClass().getName();
            String regularName = fieldName.replaceAll("Field", "");

            // get the Butcher arrays from the regular integrator
            @SuppressWarnings("unchecked")
            Constructor<RungeKuttaIntegrator> constructor =
                (Constructor<RungeKuttaIntegrator>) Class.forName(regularName).getConstructor(Double.TYPE);
            final RungeKuttaIntegrator regularIntegrator =
                            constructor.newInstance(1.0);
            double[][] regularA = regularIntegrator.getA();
            double[]   regularB = regularIntegrator.getB();
            double[]   regularC = regularIntegrator.getC();

            Assert.assertEquals(regularA.length, fieldA.length);
            for (int i = 0; i < regularA.length; ++i) {
                checkArray(regularA[i], fieldA[i]);
            }
            checkArray(regularB, fieldB);
            checkArray(regularC, fieldC);

        } catch (ClassNotFoundException | IllegalAccessException | IllegalArgumentException  |
                 SecurityException      | NoSuchMethodException  | InvocationTargetException |
                 InstantiationException e) {
            Assert.fail(e.getLocalizedMessage());
        }
    }

    private <T extends CalculusFieldElement<T>> void checkArray(double[] regularArray, T[] fieldArray) {
        Assert.assertEquals(regularArray.length, fieldArray.length);
        for (int i = 0; i < regularArray.length; ++i) {
            if (regularArray[i] == 0) {
                Assert.assertTrue(0.0 == fieldArray[i].getReal());
            } else {
                Assert.assertEquals(regularArray[i], fieldArray[i].getReal(), FastMath.ulp(regularArray[i]));
            }
        }
    }

    @Test
    public abstract void testMissedEndEvent();

    protected <T extends CalculusFieldElement<T>> void doTestMissedEndEvent(final Field<T> field,
                                                                            final double epsilonT, final double epsilonY)
        throws MathIllegalArgumentException, MathIllegalStateException {
        final T   t0     = field.getZero().add(1878250320.0000029);
        final T   tEvent = field.getZero().add(1878250379.9999986);
        final T[] k      = MathArrays.buildArray(field, 3);
        k[0] = field.getZero().add(1.0e-4);
        k[1] = field.getZero().add(1.0e-5);
        k[2] = field.getZero().add(1.0e-6);
        FieldOrdinaryDifferentialEquation<T> ode = new FieldOrdinaryDifferentialEquation<T>() {

            public int getDimension() {
                return k.length;
            }

            public T[] computeDerivatives(T t, T[] y) {
                T[] yDot = MathArrays.buildArray(field, k.length);
                for (int i = 0; i < y.length; ++i) {
                    yDot[i] = k[i].multiply(y[i]);
                }
                return yDot;
            }
        };

        RungeKuttaFieldIntegrator<T> integrator = createIntegrator(field, field.getZero().add(60.0));

        T[] y0   = MathArrays.buildArray(field, k.length);
        for (int i = 0; i < y0.length; ++i) {
            y0[i] = field.getOne().add(i);
        }

        FieldODEStateAndDerivative<T> result = integrator.integrate(new FieldExpandableODE<T>(ode),
                                                                    new FieldODEState<T>(t0, y0),
                                                                    tEvent);
        Assert.assertEquals(tEvent.getReal(), result.getTime().getReal(), epsilonT);
        T[] y = result.getPrimaryState();
        for (int i = 0; i < y.length; ++i) {
            Assert.assertEquals(y0[i].multiply(k[i].multiply(result.getTime().subtract(t0)).exp()).getReal(),
                                y[i].getReal(),
                                epsilonY);
        }

        integrator.addEventDetector(new FieldODEEventDetector<T>() {
            public FieldAdaptableInterval<T> getMaxCheckInterval() {
                return s -> Double.POSITIVE_INFINITY;
            }
            public int getMaxIterationCount() {
                return 100;
            }
            public BracketedRealFieldUnivariateSolver<T> getSolver() {
                return new FieldBracketingNthOrderBrentSolver<T>(field.getZero(),
                                                                 field.getZero().newInstance(1.0e-20),
                                                                 field.getZero(),
                                                                 5);
            }
            public T g(FieldODEStateAndDerivative<T> state) {
                return state.getTime().subtract(tEvent);
            }
            public FieldODEEventHandler<T> getHandler() {
                return (state, detector, increasing) -> {
                    Assert.assertEquals(tEvent.getReal(), state.getTime().getReal(), epsilonT);
                    return Action.CONTINUE;
                };
            }
        });
        result = integrator.integrate(new FieldExpandableODE<T>(ode),
                                      new FieldODEState<T>(t0, y0),
                                      tEvent.add(120));
        Assert.assertEquals(tEvent.add(120).getReal(), result.getTime().getReal(), epsilonT);
        y = result.getPrimaryState();
        for (int i = 0; i < y.length; ++i) {
            Assert.assertEquals(y0[i].multiply(k[i].multiply(result.getTime().subtract(t0)).exp()).getReal(),
                                y[i].getReal(),
                                epsilonY);
        }

    }

    @Test
    public abstract void testSanityChecks();

    protected <T extends CalculusFieldElement<T>> void doTestSanityChecks(Field<T> field)
        throws MathIllegalArgumentException, MathIllegalStateException {
        RungeKuttaFieldIntegrator<T> integrator = createIntegrator(field, field.getZero().add(0.01));
        try  {
            TestFieldProblem1<T> pb = new TestFieldProblem1<T>(field);
            integrator.integrate(new FieldExpandableODE<T>(pb),
                                 new FieldODEState<T>(field.getZero(), MathArrays.buildArray(field, pb.getDimension() + 10)),
                                 field.getOne());
            Assert.fail("an exception should have been thrown");
        } catch(MathIllegalArgumentException ie) {
            Assert.assertEquals(LocalizedCoreFormats.DIMENSIONS_MISMATCH, ie.getSpecifier());
        }
        try  {
            TestFieldProblem1<T> pb = new TestFieldProblem1<T>(field);
            integrator.integrate(new FieldExpandableODE<T>(pb),
                                 new FieldODEState<T>(field.getZero(), MathArrays.buildArray(field, pb.getDimension())),
                                 field.getZero());
            Assert.fail("an exception should have been thrown");
        } catch(MathIllegalArgumentException ie) {
            Assert.assertEquals(LocalizedODEFormats.TOO_SMALL_INTEGRATION_INTERVAL, ie.getSpecifier());
        }
    }

    @Test
    public abstract void testDecreasingSteps();

    protected <T extends CalculusFieldElement<T>> void doTestDecreasingSteps(Field<T> field,
                                                                         final double safetyValueFactor,
                                                                         final double safetyTimeFactor,
                                                                         final double epsilonT)
        throws MathIllegalArgumentException, MathIllegalStateException {

        @SuppressWarnings("unchecked")
        TestFieldProblemAbstract<T>[] allProblems =
                        (TestFieldProblemAbstract<T>[]) Array.newInstance(TestFieldProblemAbstract.class, 6);
        allProblems[0] = new TestFieldProblem1<T>(field);
        allProblems[1] = new TestFieldProblem2<T>(field);
        allProblems[2] = new TestFieldProblem3<T>(field);
        allProblems[3] = new TestFieldProblem4<T>(field);
        allProblems[4] = new TestFieldProblem5<T>(field);
        allProblems[5] = new TestFieldProblem6<T>(field);
        for (TestFieldProblemAbstract<T> pb :  allProblems) {

            T previousValueError = null;
            T previousTimeError  = null;
            for (int i = 4; i < 10; ++i) {

                T step = pb.getFinalTime().subtract(pb.getInitialState().getTime()).multiply(FastMath.pow(2.0, -i));

                RungeKuttaFieldIntegrator<T> integ = createIntegrator(field, step);
                TestFieldProblemHandler<T> handler = new TestFieldProblemHandler<T>(pb, integ);
                integ.addStepHandler(handler);
                final double maxCheck = Double.POSITIVE_INFINITY;
                final T eventTol = step.multiply(1.0e-6);
                FieldODEEventDetector<T>[] functions = pb.getEventDetectors(maxCheck, eventTol, 1000);
                for (int l = 0; l < functions.length; ++l) {
                    integ.addEventDetector(functions[l]);
                }
                Assert.assertEquals(functions.length, integ.getEventDetectors().size());
                FieldODEStateAndDerivative<T> stop = integ.integrate(new FieldExpandableODE<T>(pb),
                                                                     pb.getInitialState(),
                                                                     pb.getFinalTime());
                if (functions.length == 0) {
                    Assert.assertEquals(pb.getFinalTime().getReal(), stop.getTime().getReal(), epsilonT);
                }

                T error = handler.getMaximalValueError();
                if (i > 4) {
                    Assert.assertTrue(error.subtract(previousValueError.abs().multiply(safetyValueFactor)).getReal() < 0);
                }
                previousValueError = error;

                T timeError = handler.getMaximalTimeError();
                if (i > 4) {
                    // can't expect time error to be less than event finding tolerance
                    T timeTol = max(eventTol, previousTimeError.abs().multiply(safetyTimeFactor));
                    Assert.assertTrue(timeError.subtract(timeTol).getReal() <= 0);
                }
                previousTimeError = timeError;

                integ.clearEventDetectors();
                Assert.assertEquals(0, integ.getEventDetectors().size());
            }

        }

    }

    /**
     * Get the larger of two numbers.
     *
     * @param a first number.
     * @param b second number.
     * @return the larger of a and b.
     */
    private <T extends CalculusFieldElement<T>> T max(T a, T b) {
        return a.getReal() > b.getReal() ? a : b;
    }

    @Test
    public abstract void testSmallStep();

    protected <T extends CalculusFieldElement<T>> void doTestSmallStep(Field<T> field,
                                                                   final double epsilonLast,
                                                                   final double epsilonMaxValue,
                                                                   final double epsilonMaxTime,
                                                                   final String name)
         throws MathIllegalArgumentException, MathIllegalStateException {

        TestFieldProblem1<T> pb = new TestFieldProblem1<T>(field);
        T step = pb.getFinalTime().subtract(pb.getInitialState().getTime()).multiply(0.001);

        RungeKuttaFieldIntegrator<T> integ = createIntegrator(field, step);
        TestFieldProblemHandler<T> handler = new TestFieldProblemHandler<T>(pb, integ);
        integ.addStepHandler(handler);
        integ.integrate(new FieldExpandableODE<T>(pb), pb.getInitialState(), pb.getFinalTime());

        Assert.assertEquals(0, handler.getLastError().getReal(),         epsilonLast);
        Assert.assertEquals(0, handler.getMaximalValueError().getReal(), epsilonMaxValue);
        Assert.assertEquals(0, handler.getMaximalTimeError().getReal(),  epsilonMaxTime);
        Assert.assertEquals(name, integ.getName());

    }

    @Test
    public abstract void testBigStep();

    protected <T extends CalculusFieldElement<T>> void doTestBigStep(Field<T> field,
                                                                 final double belowLast,
                                                                 final double belowMaxValue,
                                                                 final double epsilonMaxTime,
                                                                 final String name)
        throws MathIllegalArgumentException, MathIllegalStateException {

        TestFieldProblem1<T> pb = new TestFieldProblem1<T>(field);
        T step = pb.getFinalTime().subtract(pb.getInitialState().getTime()).multiply(0.2);

        RungeKuttaFieldIntegrator<T> integ = createIntegrator(field, step);
        TestFieldProblemHandler<T> handler = new TestFieldProblemHandler<T>(pb, integ);
        integ.addStepHandler(handler);
        integ.integrate(new FieldExpandableODE<T>(pb), pb.getInitialState(), pb.getFinalTime());

        Assert.assertTrue(handler.getLastError().getReal()         > belowLast);
        Assert.assertTrue(handler.getMaximalValueError().getReal() > belowMaxValue);
        Assert.assertEquals(0, handler.getMaximalTimeError().getReal(),  epsilonMaxTime);
        Assert.assertEquals(name, integ.getName());

    }

    @Test
    public abstract void testBackward();

    protected <T extends CalculusFieldElement<T>> void doTestBackward(Field<T> field,
                                                                  final double epsilonLast,
                                                                  final double epsilonMaxValue,
                                                                  final double epsilonMaxTime,
                                                                  final String name)
        throws MathIllegalArgumentException, MathIllegalStateException {

        TestFieldProblem5<T> pb = new TestFieldProblem5<T>(field);
        T step = pb.getFinalTime().subtract(pb.getInitialState().getTime()).multiply(0.001).abs();

        RungeKuttaFieldIntegrator<T> integ = createIntegrator(field, step);
        TestFieldProblemHandler<T> handler = new TestFieldProblemHandler<T>(pb, integ);
        integ.addStepHandler(handler);
        integ.integrate(new FieldExpandableODE<T>(pb), pb.getInitialState(), pb.getFinalTime());

        Assert.assertEquals(0, handler.getLastError().getReal(),         epsilonLast);
        Assert.assertEquals(0, handler.getMaximalValueError().getReal(), epsilonMaxValue);
        Assert.assertEquals(0, handler.getMaximalTimeError().getReal(),  epsilonMaxTime);
        Assert.assertEquals(name, integ.getName());

    }

    @Test
    public abstract void testKepler();

    protected <T extends CalculusFieldElement<T>> void doTestKepler(Field<T> field, double expectedMaxError, double epsilon)
        throws MathIllegalArgumentException, MathIllegalStateException {

        final TestFieldProblem3<T> pb  = new TestFieldProblem3<T>(field.getZero().add(0.9));
        T step = pb.getFinalTime().subtract(pb.getInitialState().getTime()).multiply(0.0003);

        RungeKuttaFieldIntegrator<T> integ = createIntegrator(field, step);
        integ.addStepHandler(new KeplerHandler<T>(pb, expectedMaxError, epsilon));
        final FieldExpandableODE<T> expandable = new FieldExpandableODE<T>(pb);
        Assert.assertSame(pb, expandable.getPrimary());
        integ.integrate(expandable, pb.getInitialState(), pb.getFinalTime());
    }

    private static class KeplerHandler<T extends CalculusFieldElement<T>> implements FieldODEStepHandler<T> {
        private T maxError;
        private final TestFieldProblem3<T> pb;
        private final double expectedMaxError;
        private final double epsilon;
        public KeplerHandler(TestFieldProblem3<T> pb, double expectedMaxError, double epsilon) {
            this.pb               = pb;
            this.expectedMaxError = expectedMaxError;
            this.epsilon          = epsilon;
            maxError = pb.getField().getZero();
        }
        public void init(FieldODEStateAndDerivative<T> state0, T t) {
            maxError = pb.getField().getZero();
        }
        public void handleStep(FieldODEStateInterpolator<T> interpolator) {

            FieldODEStateAndDerivative<T> current = interpolator.getCurrentState();
            T[] theoreticalY  = pb.computeTheoreticalState(current.getTime());
            T dx = current.getPrimaryState()[0].subtract(theoreticalY[0]);
            T dy = current.getPrimaryState()[1].subtract(theoreticalY[1]);
            T error = dx.multiply(dx).add(dy.multiply(dy));
            if (error.subtract(maxError).getReal() > 0) {
                maxError = error;
            }
        }
        public void finish(FieldODEStateAndDerivative<T> finalState) {
            Assert.assertEquals(expectedMaxError, maxError.getReal(), epsilon);
        }
    }

    @Test
    public abstract void testStepSize();

    protected <T extends CalculusFieldElement<T>> void doTestStepSize(final Field<T> field, final double epsilon)
        throws MathIllegalArgumentException, MathIllegalStateException {
        final T finalTime = field.getZero().add(5.0);
        final T step = field.getZero().add(1.23456);
        RungeKuttaFieldIntegrator<T> integ = createIntegrator(field, step);
        integ.addStepHandler(new FieldODEStepHandler<T>() {
            public void handleStep(FieldODEStateInterpolator<T> interpolator) {
                if (interpolator.getCurrentState().getTime().subtract(finalTime).getReal() < -0.001) {
                    Assert.assertEquals(step.getReal(),
                                        interpolator.getCurrentState().getTime().subtract(interpolator.getPreviousState().getTime()).getReal(),
                                        epsilon);
                }
            }
        });
        integ.integrate(new FieldExpandableODE<T>(new FieldOrdinaryDifferentialEquation<T>() {
            public T[] computeDerivatives(T t, T[] y) {
                T[] dot = MathArrays.buildArray(t.getField(), 1);
                dot[0] = t.getField().getOne();
                return dot;
            }
            public int getDimension() {
                return 1;
            }
        }), new FieldODEState<T>(field.getZero(), MathArrays.buildArray(field, 1)), finalTime);
    }

    @Test
    public abstract void testSingleStep();

    protected <T extends CalculusFieldElement<T>> void doTestSingleStep(final Field<T> field, final double epsilon) {

        final TestFieldProblem3<T> pb  = new TestFieldProblem3<T>(field.getZero().add(0.9));
        T h = pb.getFinalTime().subtract(pb.getInitialState().getTime()).multiply(0.0003);

        RungeKuttaFieldIntegrator<T> integ = createIntegrator(field, field.getZero().add(Double.NaN));
        T   t = pb.getInitialState().getTime();
        T[] y = pb.getInitialState().getPrimaryState();
        for (int i = 0; i < 100; ++i) {
            y = integ.singleStep(pb, t, y, t.add(h));
            t = t.add(h);
        }
        T[] yth = pb.computeTheoreticalState(t);
        T dx = y[0].subtract(yth[0]);
        T dy = y[1].subtract(yth[1]);
        T error = dx.multiply(dx).add(dy.multiply(dy));
        Assert.assertEquals(0.0, error.getReal(), epsilon);
    }

    @Test
    public abstract void testTooLargeFirstStep();

    protected <T extends CalculusFieldElement<T>> void doTestTooLargeFirstStep(final Field<T> field) {

        RungeKuttaFieldIntegrator<T> integ = createIntegrator(field, field.getZero().add(0.5));
        final T t0 = field.getZero();
        final T[] y0 = MathArrays.buildArray(field, 1);
        y0[0] = field.getOne();
        final T t   = field.getZero().add(0.001);
        FieldOrdinaryDifferentialEquation<T> equations = new FieldOrdinaryDifferentialEquation<T>() {

            public int getDimension() {
                return 1;
            }

            public T[] computeDerivatives(T t, T[] y) {
                Assert.assertTrue(t.getReal() >= FastMath.nextAfter(t0.getReal(), Double.NEGATIVE_INFINITY));
                Assert.assertTrue(t.getReal() <= FastMath.nextAfter(t.getReal(),   Double.POSITIVE_INFINITY));
                T[] yDot = MathArrays.buildArray(field, 1);
                yDot[0] = y[0].multiply(-100.0);
                return yDot;
            }

        };

        integ.integrate(new FieldExpandableODE<T>(equations), new FieldODEState<T>(t0, y0), t);

    }

    @Test
    public abstract void testUnstableDerivative();

    protected <T extends CalculusFieldElement<T>> void doTestUnstableDerivative(Field<T> field, double epsilon) {
      final StepFieldProblem<T> stepProblem = new StepFieldProblem<T>(field,
                                                                      s -> 999.0,
                                                                      field.getZero().newInstance(1.0e+12),
                                                                      1000000,
                                                                      field.getZero().newInstance(0.0),
                                                                      field.getZero().newInstance(1.0),
                                                                      field.getZero().newInstance(2.0)).
                                              withMaxCheck(field.getZero().newInstance(1.0)).
                                              withMaxIter(1000).
                                              withThreshold(field.getZero().newInstance(1.0e-12));
      Assert.assertEquals(1.0,     stepProblem.getMaxCheckInterval().currentInterval(null), 1.0e-15);
      Assert.assertEquals(1000,    stepProblem.getMaxIterationCount());
      Assert.assertEquals(1.0e-12, stepProblem.getSolver().getAbsoluteAccuracy().getReal(), 1.0e-25);
      Assert.assertNotNull(stepProblem.getHandler());
      RungeKuttaFieldIntegrator<T> integ = createIntegrator(field, field.getZero().add(0.3));
      integ.addEventDetector(stepProblem);
      FieldODEStateAndDerivative<T> result = integ.integrate(new FieldExpandableODE<T>(stepProblem),
                                                             new FieldODEState<T>(field.getZero(), MathArrays.buildArray(field, 1)),
                                                             field.getZero().add(10.0));
      Assert.assertEquals(8.0, result.getPrimaryState()[0].getReal(), epsilon);
    }

    @Test
    public abstract void testDerivativesConsistency();

    protected <T extends CalculusFieldElement<T>> void doTestDerivativesConsistency(final Field<T> field, double epsilon) {
        TestFieldProblem3<T> pb = new TestFieldProblem3<T>(field);
        T step = pb.getFinalTime().subtract(pb.getInitialState().getTime()).multiply(0.001);
        RungeKuttaFieldIntegrator<T> integ = createIntegrator(field, step);
        StepInterpolatorTestUtils.checkDerivativesConsistency(integ, pb, 1.0e-10);
    }

    @Test
    public abstract void testPartialDerivatives();

    protected void doTestPartialDerivatives(final double epsilonY,
                                            final double[] epsilonPartials) {

        // parameters indices
        final DSFactory factory = new DSFactory(5, 1);
        final int parOmega   = 0;
        final int parTO      = 1;
        final int parY00     = 2;
        final int parY01     = 3;
        final int parT       = 4;

        DerivativeStructure omega = factory.variable(parOmega, 1.3);
        DerivativeStructure t0    = factory.variable(parTO, 1.3);
        DerivativeStructure[] y0  = new DerivativeStructure[] {
            factory.variable(parY00, 3.0),
            factory.variable(parY01, 4.0)
        };
        DerivativeStructure t     = factory.variable(parT, 6.0);
        SinCos sinCos = new SinCos(omega);

        RungeKuttaFieldIntegrator<DerivativeStructure> integrator =
                        createIntegrator(omega.getField(), t.subtract(t0).multiply(0.001));
        FieldODEStateAndDerivative<DerivativeStructure> result =
                        integrator.integrate(new FieldExpandableODE<DerivativeStructure>(sinCos),
                                             new FieldODEState<DerivativeStructure>(t0, y0),
                                             t);

        // check values
        for (int i = 0; i < sinCos.getDimension(); ++i) {
            Assert.assertEquals(sinCos.theoreticalY(t.getReal())[i], result.getPrimaryState()[i].getValue(), epsilonY);
        }

        // check derivatives
        final double[][] derivatives = sinCos.getDerivatives(t.getReal());
        for (int i = 0; i < sinCos.getDimension(); ++i) {
            for (int parameter = 0; parameter < factory.getCompiler().getFreeParameters(); ++parameter) {
                Assert.assertEquals(derivatives[i][parameter],
                                    dYdP(result.getPrimaryState()[i], parameter),
                                    epsilonPartials[parameter]);
            }
        }

    }

    @Test
    public abstract void testSecondaryEquations();

    protected <T extends CalculusFieldElement<T>> void doTestSecondaryEquations(final Field<T> field,
                                                                                final double epsilonSinCos,
                                                                                final double epsilonLinear) {
        FieldOrdinaryDifferentialEquation<T> sinCos = new FieldOrdinaryDifferentialEquation<T>() {

            @Override
            public int getDimension() {
                return 2;
            }

            @Override
            public T[] computeDerivatives(T t, T[] y) {
                // here, we compute only half of the derivative
                // we will compute the full derivatives by multiplying
                // the main equation from within the additional equation
                // it is not the proper way, but it is intended to check
                // additional equations *can* change main equation
                T[] yDot = y.clone();
                yDot[0] = y[1].multiply( 0.5);
                yDot[1] = y[0].multiply(-0.5);
                return yDot;
            }

        };

        FieldSecondaryODE<T> linear = new FieldSecondaryODE<T>() {

            @Override
            public int getDimension() {
                return 1;
            }

            @Override
            public T[] computeDerivatives(T t, T[] primary, T[] primaryDot, T[] secondary) {
                for (int i = 0; i < primaryDot.length; ++i) {
                    // this secondary equation also changes the primary state derivative
                    // a proper example of this is for example optimal control when
                    // the secondary equations handle co-state, which changes control,
                    // and the control changes the primary state
                    primaryDot[i] = primaryDot[i].multiply(2);
                }
                T[] secondaryDot = secondary.clone();
                secondaryDot[0] = t.getField().getOne().negate();
                return secondaryDot;
            }

        };

        FieldExpandableODE<T> expandable = new FieldExpandableODE<>(sinCos);
        expandable.addSecondaryEquations(linear);

        FieldODEIntegrator<T> integrator = createIntegrator(field, field.getZero().add(0.001));
        final double[] max = new double[2];
        integrator.addStepHandler(new FieldODEStepHandler<T>() {
            @Override
            public void handleStep(FieldODEStateInterpolator<T> interpolator) {
                for (int i = 0; i <= 10; ++i) {
                    T tPrev = interpolator.getPreviousState().getTime();
                    T tCurr = interpolator.getCurrentState().getTime();
                    T t     = tPrev.multiply(10 - i).add(tCurr.multiply(i)).divide(10);
                    FieldODEStateAndDerivative<T> state = interpolator.getInterpolatedState(t);
                    Assert.assertEquals(2, state.getPrimaryStateDimension());
                    Assert.assertEquals(1, state.getNumberOfSecondaryStates());
                    Assert.assertEquals(2, state.getSecondaryStateDimension(0));
                    Assert.assertEquals(1, state.getSecondaryStateDimension(1));
                    Assert.assertEquals(3, state.getCompleteStateDimension());
                    max[0] = FastMath.max(max[0],
                                          t.sin().subtract(state.getPrimaryState()[0]).norm());
                    max[0] = FastMath.max(max[0],
                                          t.cos().subtract(state.getPrimaryState()[1]).norm());
                    max[1] = FastMath.max(max[1],
                                          field.getOne().subtract(t).subtract(state.getSecondaryState(1)[0]).norm());
                }
            }
        });

        T[] primary0 = MathArrays.buildArray(field, 2);
        primary0[0] = field.getZero();
        primary0[1] = field.getOne();
        T[][] secondary0 = MathArrays.buildArray(field, 1, 1);
        secondary0[0][0] = field.getOne();
        FieldODEState<T> initialState = new FieldODEState<T>(field.getZero(), primary0, secondary0);

        FieldODEStateAndDerivative<T> finalState =
                        integrator.integrate(expandable, initialState, field.getZero().add(10.0));
        Assert.assertEquals(10.0, finalState.getTime().getReal(), 1.0e-12);
        Assert.assertEquals(0, max[0], epsilonSinCos);
        Assert.assertEquals(0, max[1], epsilonLinear);

    }

    private double dYdP(final DerivativeStructure y, final int parameter) {
        int[] orders = new int[y.getFreeParameters()];
        orders[parameter] = 1;
        return y.getPartialDerivative(orders);
    }

    private static class SinCos implements FieldOrdinaryDifferentialEquation<DerivativeStructure> {

        private final DerivativeStructure omega;
        private       DerivativeStructure r;
        private       DerivativeStructure alpha;

        private double dRdY00;
        private double dRdY01;
        private double dAlphadOmega;
        private double dAlphadT0;
        private double dAlphadY00;
        private double dAlphadY01;

        protected SinCos(final DerivativeStructure omega) {
            this.omega = omega;
        }

        public int getDimension() {
            return 2;
        }

        public void init(final DerivativeStructure t0, final DerivativeStructure[] y0,
                         final DerivativeStructure finalTime) {

            // theoretical solution is y(t) = { r * sin(omega * t + alpha), r * cos(omega * t + alpha) }
            // so we retrieve alpha by identification from the initial state
            final DerivativeStructure r2 = y0[0].multiply(y0[0]).add(y0[1].multiply(y0[1]));

            this.r            = r2.sqrt();
            this.dRdY00       = y0[0].divide(r).getReal();
            this.dRdY01       = y0[1].divide(r).getReal();

            this.alpha        = y0[0].atan2(y0[1]).subtract(t0.multiply(omega));
            this.dAlphadOmega = -t0.getReal();
            this.dAlphadT0    = -omega.getReal();
            this.dAlphadY00   = y0[1].divide(r2).getReal();
            this.dAlphadY01   = y0[0].negate().divide(r2).getReal();

        }

        public DerivativeStructure[] computeDerivatives(final DerivativeStructure t, final DerivativeStructure[] y) {
            return new DerivativeStructure[] {
                omega.multiply(y[1]),
                omega.multiply(y[0]).negate()
            };
        }

        public double[] theoreticalY(final double t) {
            final double theta = omega.getReal() * t + alpha.getReal();
            return new double[] {
                r.getReal() * FastMath.sin(theta), r.getReal() * FastMath.cos(theta)
            };
        }

        public double[][] getDerivatives(final double t) {

            // intermediate angle and state
            final double theta        = omega.getReal() * t + alpha.getReal();
            final double sin          = FastMath.sin(theta);
            final double cos          = FastMath.cos(theta);
            final double y0           = r.getReal() * sin;
            final double y1           = r.getReal() * cos;

            // partial derivatives of the state first component
            final double dY0dOmega    =                y1 * (t + dAlphadOmega);
            final double dY0dT0       =                y1 * dAlphadT0;
            final double dY0dY00      = dRdY00 * sin + y1 * dAlphadY00;
            final double dY0dY01      = dRdY01 * sin + y1 * dAlphadY01;
            final double dY0dT        =                y1 * omega.getReal();

            // partial derivatives of the state second component
            final double dY1dOmega    =              - y0 * (t + dAlphadOmega);
            final double dY1dT0       =              - y0 * dAlphadT0;
            final double dY1dY00      = dRdY00 * cos - y0 * dAlphadY00;
            final double dY1dY01      = dRdY01 * cos - y0 * dAlphadY01;
            final double dY1dT        =              - y0 * omega.getReal();

            return new double[][] {
                { dY0dOmega, dY0dT0, dY0dY00, dY0dY01, dY0dT },
                { dY1dOmega, dY1dT0, dY1dY00, dY1dY01, dY1dT }
            };

        }

    }

    @Test
    public void testIssue250() {
        final Gradient defaultStep = Gradient.constant(3, 60.);
        RungeKuttaFieldIntegrator<Gradient> integrator = createIntegrator(defaultStep.getField(), defaultStep);
        Assert.assertEquals(defaultStep.getReal(), integrator.getDefaultStep().getReal(), 0.);
    }

}
