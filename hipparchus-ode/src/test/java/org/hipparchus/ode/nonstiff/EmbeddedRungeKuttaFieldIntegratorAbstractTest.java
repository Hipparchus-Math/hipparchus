/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.hipparchus.ode.nonstiff;


import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;

import org.hipparchus.Field;
import org.hipparchus.RealFieldElement;
import org.hipparchus.analysis.differentiation.DSFactory;
import org.hipparchus.analysis.differentiation.DerivativeStructure;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.ode.FieldExpandableODE;
import org.hipparchus.ode.FieldODEState;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.ode.FieldOrdinaryDifferentialEquation;
import org.hipparchus.ode.FieldSecondaryODE;
import org.hipparchus.ode.FieldODEIntegrator;
import org.hipparchus.ode.TestFieldProblem1;
import org.hipparchus.ode.TestFieldProblem3;
import org.hipparchus.ode.TestFieldProblem4;
import org.hipparchus.ode.TestFieldProblem5;
import org.hipparchus.ode.TestFieldProblemHandler;
import org.hipparchus.ode.events.Action;
import org.hipparchus.ode.events.FieldODEEventHandler;
import org.hipparchus.ode.sampling.FieldODEStepHandler;
import org.hipparchus.ode.sampling.FieldODEStateInterpolator;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;
import org.junit.Assert;
import org.junit.Test;

public abstract class EmbeddedRungeKuttaFieldIntegratorAbstractTest {

    protected abstract <T extends RealFieldElement<T>> EmbeddedRungeKuttaFieldIntegrator<T>
    createIntegrator(Field<T> field, final double minStep, final double maxStep,
                     final double scalAbsoluteTolerance, final double scalRelativeTolerance);

    protected abstract <T extends RealFieldElement<T>> EmbeddedRungeKuttaFieldIntegrator<T>
    createIntegrator(Field<T> field, final double minStep, final double maxStep,
                     final double[] vecAbsoluteTolerance, final double[] vecRelativeTolerance);

    @Test
    public abstract void testNonFieldIntegratorConsistency();

    protected <T extends RealFieldElement<T>> void doTestNonFieldIntegratorConsistency(final Field<T> field) {
        try {

            // get the Butcher arrays from the field integrator
            EmbeddedRungeKuttaFieldIntegrator<T> fieldIntegrator = createIntegrator(field, 0.001, 1.0, 1.0, 1.0);
            T[][] fieldA = fieldIntegrator.getA();
            T[]   fieldB = fieldIntegrator.getB();
            T[]   fieldC = fieldIntegrator.getC();

            String fieldName   = fieldIntegrator.getClass().getName();
            String regularName = fieldName.replaceAll("Field", "");

            // get the Butcher arrays from the regular integrator
            @SuppressWarnings("unchecked")
            Constructor<EmbeddedRungeKuttaIntegrator> constructor =
                (Constructor<EmbeddedRungeKuttaIntegrator>) Class.forName(regularName).getConstructor(Double.TYPE,
                                                                                                      Double.TYPE,
                                                                                                      Double.TYPE,
                                                                                                      Double.TYPE);
            final EmbeddedRungeKuttaIntegrator regularIntegrator =
                            constructor.newInstance(0.0, 1.0, 1.0e-10, 1.0e-10);
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

    private <T extends RealFieldElement<T>> void checkArray(double[] regularArray, T[] fieldArray) {
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
    public abstract void testForwardBackwardExceptions();

    protected <T extends RealFieldElement<T>> void doTestForwardBackwardExceptions(final Field<T> field) {
        FieldOrdinaryDifferentialEquation<T> equations = new FieldOrdinaryDifferentialEquation<T>() {

            public int getDimension() {
                return 1;
            }

            public void init(T t0, T[] y0, T t) {
            }

            public T[] computeDerivatives(T t, T[] y) {
                if (t.getReal() < -0.5) {
                    throw new LocalException();
                } else {
                    throw new RuntimeException("oops");
                }
            }
        };

        EmbeddedRungeKuttaFieldIntegrator<T> integrator = createIntegrator(field, 0.0, 1.0, 1.0e-10, 1.0e-10);

        try  {
            integrator.integrate(new FieldExpandableODE<T>(equations),
                                 new FieldODEState<T>(field.getOne().negate(),
                                                      MathArrays.buildArray(field, 1)),
                                 field.getZero());
            Assert.fail("an exception should have been thrown");
          } catch(LocalException de) {
            // expected behavior
          }

          try  {
              integrator.integrate(new FieldExpandableODE<T>(equations),
                                   new FieldODEState<T>(field.getZero(),
                                                        MathArrays.buildArray(field, 1)),
                                   field.getOne());
               Assert.fail("an exception should have been thrown");
          } catch(RuntimeException de) {
            // expected behavior
          }
    }

    protected static class LocalException extends RuntimeException {
        private static final long serialVersionUID = 20151208L;
    }

    @Test(expected=MathIllegalArgumentException.class)
    public abstract void testMinStep();

    protected <T extends RealFieldElement<T>> void doTestMinStep(final Field<T> field)
        throws MathIllegalArgumentException {

        TestFieldProblem1<T> pb = new TestFieldProblem1<T>(field);
        double minStep = pb.getFinalTime().subtract(pb.getInitialState().getTime()).multiply(0.1).getReal();
        double maxStep = pb.getFinalTime().subtract(pb.getInitialState().getTime()).getReal();
        double[] vecAbsoluteTolerance = { 1.0e-15, 1.0e-16 };
        double[] vecRelativeTolerance = { 1.0e-15, 1.0e-16 };

        FieldODEIntegrator<T> integ = createIntegrator(field, minStep, maxStep,
                                                              vecAbsoluteTolerance, vecRelativeTolerance);
        TestFieldProblemHandler<T> handler = new TestFieldProblemHandler<T>(pb, integ);
        integ.addStepHandler(handler);
        integ.integrate(new FieldExpandableODE<T>(pb), pb.getInitialState(), pb.getFinalTime());
        Assert.fail("an exception should have been thrown");

    }

    @Test
    public abstract void testIncreasingTolerance();

    protected <T extends RealFieldElement<T>> void doTestIncreasingTolerance(final Field<T> field,
                                                                             double factor,
                                                                             double epsilon) {

        int previousCalls = Integer.MAX_VALUE;
        for (int i = -12; i < -2; ++i) {
            TestFieldProblem1<T> pb = new TestFieldProblem1<T>(field);
            double minStep = 0;
            double maxStep = pb.getFinalTime().subtract(pb.getInitialState().getTime()).getReal();
            double scalAbsoluteTolerance = FastMath.pow(10.0, i);
            double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

            FieldODEIntegrator<T> integ = createIntegrator(field, minStep, maxStep,
                                                                  scalAbsoluteTolerance, scalRelativeTolerance);
            TestFieldProblemHandler<T> handler = new TestFieldProblemHandler<T>(pb, integ);
            integ.addStepHandler(handler);
            integ.integrate(new FieldExpandableODE<T>(pb), pb.getInitialState(), pb.getFinalTime());

            Assert.assertTrue(handler.getMaximalValueError().getReal() < (factor * scalAbsoluteTolerance));
            Assert.assertEquals(0, handler.getMaximalTimeError().getReal(), epsilon);

            int calls = pb.getCalls();
            Assert.assertEquals(integ.getEvaluations(), calls);
            Assert.assertTrue(calls <= previousCalls);
            previousCalls = calls;

        }

    }

    @Test
    public abstract void testEvents();

    protected <T extends RealFieldElement<T>> void doTestEvents(final Field<T> field,
                                                                final double epsilonMaxValue,
                                                                final String name) {

      TestFieldProblem4<T> pb = new TestFieldProblem4<T>(field);
      double minStep = 0;
      double maxStep = pb.getFinalTime().subtract(pb.getInitialState().getTime()).getReal();
      double scalAbsoluteTolerance = 1.0e-8;
      double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

      FieldODEIntegrator<T> integ = createIntegrator(field, minStep, maxStep,
                                                            scalAbsoluteTolerance, scalRelativeTolerance);
      TestFieldProblemHandler<T> handler = new TestFieldProblemHandler<T>(pb, integ);
      integ.addStepHandler(handler);
      FieldODEEventHandler<T>[] functions = pb.getEventsHandlers();
      double convergence = 1.0e-8 * maxStep;
      for (int l = 0; l < functions.length; ++l) {
          integ.addEventHandler(functions[l], Double.POSITIVE_INFINITY, convergence, 1000);
      }
      Assert.assertEquals(functions.length, integ.getEventHandlers().size());
      integ.integrate(new FieldExpandableODE<T>(pb), pb.getInitialState(), pb.getFinalTime());

      Assert.assertEquals(0, handler.getMaximalValueError().getReal(), epsilonMaxValue);
      Assert.assertEquals(0, handler.getMaximalTimeError().getReal(), convergence);
      Assert.assertEquals(12.0, handler.getLastTime().getReal(), convergence);
      Assert.assertEquals(name, integ.getName());
      integ.clearEventHandlers();
      Assert.assertEquals(0, integ.getEventHandlers().size());

    }

    @Test(expected=LocalException.class)
    public abstract void testEventsErrors();

    protected <T extends RealFieldElement<T>> void doTestEventsErrors(final Field<T> field)
        throws LocalException {
        final TestFieldProblem1<T> pb = new TestFieldProblem1<T>(field);
        double minStep = 0;
        double maxStep = pb.getFinalTime().subtract(pb.getInitialState().getTime()).getReal();
        double scalAbsoluteTolerance = 1.0e-8;
        double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

        FieldODEIntegrator<T> integ = createIntegrator(field, minStep, maxStep,
                                                              scalAbsoluteTolerance, scalRelativeTolerance);
        TestFieldProblemHandler<T> handler = new TestFieldProblemHandler<T>(pb, integ);
        integ.addStepHandler(handler);

        integ.addEventHandler(new FieldODEEventHandler<T>() {
          public void init(FieldODEStateAndDerivative<T> state0, T t) {
          }
          public Action eventOccurred(FieldODEStateAndDerivative<T> state, boolean increasing) {
            return Action.CONTINUE;
          }
          public T g(FieldODEStateAndDerivative<T> state) {
            T middle = pb.getInitialState().getTime().add(pb.getFinalTime()).multiply(0.5);
            T offset = state.getTime().subtract(middle);
            if (offset.getReal() > 0) {
              throw new LocalException();
            }
            return offset;
          }
          public FieldODEState<T> resetState(FieldODEStateAndDerivative<T> state) {
              return state;
          }
        }, Double.POSITIVE_INFINITY, 1.0e-8 * maxStep, 1000);

        integ.integrate(new FieldExpandableODE<T>(pb), pb.getInitialState(), pb.getFinalTime());

    }

    @Test
    public abstract void testEventsNoConvergence();

    protected <T extends RealFieldElement<T>> void doTestEventsNoConvergence(final Field<T> field){

        final TestFieldProblem1<T> pb = new TestFieldProblem1<T>(field);
        double minStep = 0;
        double maxStep = pb.getFinalTime().subtract(pb.getInitialState().getTime()).getReal();
        double scalAbsoluteTolerance = 1.0e-8;
        double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

        FieldODEIntegrator<T> integ = createIntegrator(field, minStep, maxStep,
                                                              scalAbsoluteTolerance, scalRelativeTolerance);
        TestFieldProblemHandler<T> handler = new TestFieldProblemHandler<T>(pb, integ);
        integ.addStepHandler(handler);

        integ.addEventHandler(new FieldODEEventHandler<T>() {
            public void init(FieldODEStateAndDerivative<T> state0, T t) {
            }
            public Action eventOccurred(FieldODEStateAndDerivative<T> state, boolean increasing) {
                return Action.CONTINUE;
            }
            public T g(FieldODEStateAndDerivative<T> state) {
                T middle = pb.getInitialState().getTime().add(pb.getFinalTime()).multiply(0.5);
                T offset = state.getTime().subtract(middle);
                return (offset.getReal() > 0) ? offset.add(0.5) : offset.subtract(0.5);
            }
            public FieldODEState<T> resetState(FieldODEStateAndDerivative<T> state) {
                return state;
            }
        }, Double.POSITIVE_INFINITY, 1.0e-8 * maxStep, 3);

        try {
            integ.integrate(new FieldExpandableODE<T>(pb), pb.getInitialState(), pb.getFinalTime());
            Assert.fail("an exception should have been thrown");
        } catch (MathIllegalStateException mcee) {
            // Expected.
        }

    }

    @Test
    public abstract void testSanityChecks();

    protected <T extends RealFieldElement<T>> void doTestSanityChecks(Field<T> field) {
        TestFieldProblem3<T> pb = new TestFieldProblem3<T>(field);
        try  {
            EmbeddedRungeKuttaFieldIntegrator<T> integrator = createIntegrator(field, 0,
                                                                               pb.getFinalTime().subtract(pb.getInitialState().getTime()).getReal(),
                                                                               new double[4], new double[4]);
            integrator.integrate(new FieldExpandableODE<T>(pb),
                                 new FieldODEState<T>(pb.getInitialState().getTime(),
                                                      MathArrays.buildArray(field, 6)),
                                 pb.getFinalTime());
            Assert.fail("an exception should have been thrown");
        } catch(MathIllegalArgumentException ie) {
        }
        try  {
            EmbeddedRungeKuttaFieldIntegrator<T> integrator =
                            createIntegrator(field, 0,
                                             pb.getFinalTime().subtract(pb.getInitialState().getTime()).getReal(),
                                             new double[2], new double[4]);
            integrator.integrate(new FieldExpandableODE<T>(pb), pb.getInitialState(), pb.getFinalTime());
            Assert.fail("an exception should have been thrown");
        } catch(MathIllegalArgumentException ie) {
        }
        try  {
            EmbeddedRungeKuttaFieldIntegrator<T> integrator =
                            createIntegrator(field, 0,
                                             pb.getFinalTime().subtract(pb.getInitialState().getTime()).getReal(),
                                             new double[4], new double[4]);
            integrator.integrate(new FieldExpandableODE<T>(pb), pb.getInitialState(), pb.getInitialState().getTime());
            Assert.fail("an exception should have been thrown");
        } catch(MathIllegalArgumentException ie) {
        }
    }

    @Test
    public abstract void testBackward();

    protected <T extends RealFieldElement<T>> void doTestBackward(Field<T> field,
                                                                  final double epsilonLast,
                                                                  final double epsilonMaxValue,
                                                                  final double epsilonMaxTime,
                                                                  final String name)
        throws MathIllegalArgumentException, MathIllegalStateException {

        TestFieldProblem5<T> pb = new TestFieldProblem5<T>(field);
        double minStep = 0;
        double maxStep = pb.getFinalTime().subtract(pb.getInitialState().getTime()).abs().getReal();
        double scalAbsoluteTolerance = 1.0e-8;
        double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

        EmbeddedRungeKuttaFieldIntegrator<T> integ = createIntegrator(field, minStep, maxStep,
                                                                      scalAbsoluteTolerance,
                                                                      scalRelativeTolerance);
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

    protected <T extends RealFieldElement<T>> void doTestKepler(Field<T> field, double epsilon) {

        final TestFieldProblem3<T> pb  = new TestFieldProblem3<T>(field.getZero().add(0.9));
        double minStep = 0;
        double maxStep = pb.getFinalTime().subtract(pb.getInitialState().getTime()).getReal();
        double[] vecAbsoluteTolerance = { 1.0e-8, 1.0e-8, 1.0e-10, 1.0e-10 };
        double[] vecRelativeTolerance = { 1.0e-10, 1.0e-10, 1.0e-8, 1.0e-8 };

        FieldODEIntegrator<T> integ = createIntegrator(field, minStep, maxStep,
                                                              vecAbsoluteTolerance, vecRelativeTolerance);
        integ.addStepHandler(new KeplerHandler<T>(pb, epsilon));
        integ.integrate(new FieldExpandableODE<T>(pb), pb.getInitialState(), pb.getFinalTime());
    }

    private static class KeplerHandler<T extends RealFieldElement<T>> implements FieldODEStepHandler<T> {
        private T maxError;
        private final TestFieldProblem3<T> pb;
        private final double epsilon;
        public KeplerHandler(TestFieldProblem3<T> pb, double epsilon) {
            this.pb      = pb;
            this.epsilon = epsilon;
            maxError = pb.getField().getZero();
        }
        public void init(FieldODEStateAndDerivative<T> state0, T t) {
            maxError = pb.getField().getZero();
        }
        public void handleStep(FieldODEStateInterpolator<T> interpolator, boolean isLast)
                        throws MathIllegalStateException {

            FieldODEStateAndDerivative<T> current = interpolator.getCurrentState();
            T[] theoreticalY  = pb.computeTheoreticalState(current.getTime());
            T dx = current.getPrimaryState()[0].subtract(theoreticalY[0]);
            T dy = current.getPrimaryState()[1].subtract(theoreticalY[1]);
            T error = dx.multiply(dx).add(dy.multiply(dy));
            if (error.subtract(maxError).getReal() > 0) {
                maxError = error;
            }
            if (isLast) {
                Assert.assertEquals(0.0, maxError.getReal(), epsilon);
            }
        }
    }

    @Test
    public abstract void testSecondaryEquations();

    protected <T extends RealFieldElement<T>> void doTestSecondaryEquations(final Field<T> field,
                                                                            final double epsilonSinCos,
                                                                            final double epsilonLinear) {
        FieldOrdinaryDifferentialEquation<T> sinCos = new FieldOrdinaryDifferentialEquation<T>() {

            @Override
            public int getDimension() {
                return 2;
            }

            @Override
            public T[] computeDerivatives(T t, T[] y) {
                T[] yDot = y.clone();
                yDot[0] = y[1];
                yDot[1] = y[0].negate();
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
                T[] secondaryDot = secondary.clone();
                secondaryDot[0] = t.getField().getOne().negate();
                return secondaryDot;
            }

        };

        FieldExpandableODE<T> expandable = new FieldExpandableODE<>(sinCos);
        expandable.addSecondaryEquations(linear);

        FieldODEIntegrator<T> integrator = createIntegrator(field, 0.001, 1.0, 1.0e-12, 1.0e-12);
        final double[] max = new double[2];
        integrator.addStepHandler(new FieldODEStepHandler<T>() {
            @Override
            public void handleStep(FieldODEStateInterpolator<T> interpolator, boolean isLast) {
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
                                          t.sin().subtract(state.getPrimaryState()[0]).abs().getReal());
                    max[0] = FastMath.max(max[0],
                                          t.cos().subtract(state.getPrimaryState()[1]).abs().getReal());
                    max[1] = FastMath.max(max[1],
                                          field.getOne().subtract(t).subtract(state.getSecondaryState(1)[0]).abs().getReal());
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

        EmbeddedRungeKuttaFieldIntegrator<DerivativeStructure> integrator =
                        createIntegrator(omega.getField(),
                                         t.subtract(t0).multiply(0.001).getReal(), t.subtract(t0).getReal(),
                                         1.0e-12, 1.0e-12);
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
                Assert.assertEquals(derivatives[i][parameter], dYdP(result.getPrimaryState()[i], parameter), epsilonPartials[parameter]);
            }
        }

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

}
