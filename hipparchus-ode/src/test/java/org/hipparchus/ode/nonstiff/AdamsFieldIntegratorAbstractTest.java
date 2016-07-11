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


import org.hipparchus.Field;
import org.hipparchus.RealFieldElement;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.ode.AbstractFieldIntegrator;
import org.hipparchus.ode.FieldExpandableODE;
import org.hipparchus.ode.FieldODEState;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.ode.FieldOrdinaryDifferentialEquation;
import org.hipparchus.ode.FieldSecondaryODE;
import org.hipparchus.ode.FieldODEIntegrator;
import org.hipparchus.ode.MultistepFieldIntegrator;
import org.hipparchus.ode.TestFieldProblem1;
import org.hipparchus.ode.TestFieldProblem5;
import org.hipparchus.ode.TestFieldProblem6;
import org.hipparchus.ode.TestFieldProblemAbstract;
import org.hipparchus.ode.TestFieldProblemHandler;
import org.hipparchus.ode.sampling.FieldODEStepHandler;
import org.hipparchus.ode.sampling.FieldODEStateInterpolator;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;
import org.junit.Assert;
import org.junit.Test;

public abstract class AdamsFieldIntegratorAbstractTest {

    protected abstract <T extends RealFieldElement<T>> AdamsFieldIntegrator<T>
    createIntegrator(Field<T> field, final int nSteps, final double minStep, final double maxStep,
                     final double scalAbsoluteTolerance, final double scalRelativeTolerance);

    protected abstract <T extends RealFieldElement<T>> AdamsFieldIntegrator<T>
    createIntegrator(Field<T> field, final int nSteps, final double minStep, final double maxStep,
                     final double[] vecAbsoluteTolerance, final double[] vecRelativeTolerance);

    @Test(expected=MathIllegalArgumentException.class)
    public abstract void testMinStep();

    protected <T extends RealFieldElement<T>> void doDimensionCheck(final Field<T> field) {
        TestFieldProblem1<T> pb = new TestFieldProblem1<T>(field);

        double minStep = pb.getFinalTime().subtract(pb.getInitialState().getTime()).multiply(0.1).getReal();
        double maxStep = pb.getFinalTime().subtract(pb.getInitialState().getTime()).getReal();
        double[] vecAbsoluteTolerance = { 1.0e-15, 1.0e-16 };
        double[] vecRelativeTolerance = { 1.0e-15, 1.0e-16 };

        FieldODEIntegrator<T> integ = createIntegrator(field, 4, minStep, maxStep,
                                                              vecAbsoluteTolerance,
                                                              vecRelativeTolerance);
        TestFieldProblemHandler<T> handler = new TestFieldProblemHandler<T>(pb, integ);
        integ.addStepHandler(handler);
        integ.integrate(new FieldExpandableODE<T>(pb), pb.getInitialState(), pb.getFinalTime());

    }

    @Test
    public abstract void testIncreasingTolerance();

    protected <T extends RealFieldElement<T>> void doTestIncreasingTolerance(final Field<T> field,
                                                                             double ratioMin, double ratioMax) {

        int previousCalls = Integer.MAX_VALUE;
        for (int i = -12; i < -2; ++i) {
            TestFieldProblem1<T> pb = new TestFieldProblem1<T>(field);
            double minStep = 0;
            double maxStep = pb.getFinalTime().subtract(pb.getInitialState().getTime()).getReal();
            double scalAbsoluteTolerance = FastMath.pow(10.0, i);
            double scalRelativeTolerance = 0.01 * scalAbsoluteTolerance;

            FieldODEIntegrator<T> integ = createIntegrator(field, 4, minStep, maxStep,
                                                                  scalAbsoluteTolerance,
                                                                  scalRelativeTolerance);
            TestFieldProblemHandler<T> handler = new TestFieldProblemHandler<T>(pb, integ);
            integ.addStepHandler(handler);
            integ.integrate(new FieldExpandableODE<T>(pb), pb.getInitialState(), pb.getFinalTime());

            Assert.assertTrue(handler.getMaximalValueError().getReal() > ratioMin * scalAbsoluteTolerance);
            Assert.assertTrue(handler.getMaximalValueError().getReal() < ratioMax * scalAbsoluteTolerance);

            int calls = pb.getCalls();
            Assert.assertEquals(integ.getEvaluations(), calls);
            Assert.assertTrue(calls <= previousCalls);
            previousCalls = calls;

        }

    }

    @Test(expected = MathIllegalStateException.class)
    public abstract void exceedMaxEvaluations();

    protected <T extends RealFieldElement<T>> void doExceedMaxEvaluations(final Field<T> field, final int max) {

        TestFieldProblem1<T> pb  = new TestFieldProblem1<T>(field);
        double range = pb.getFinalTime().subtract(pb.getInitialState().getTime()).getReal();

        FieldODEIntegrator<T> integ = createIntegrator(field, 2, 0, range, 1.0e-12, 1.0e-12);
        TestFieldProblemHandler<T> handler = new TestFieldProblemHandler<T>(pb, integ);
        integ.addStepHandler(handler);
        integ.setMaxEvaluations(max);
        integ.integrate(new FieldExpandableODE<T>(pb), pb.getInitialState(), pb.getFinalTime());

    }

    @Test
    public abstract void backward();

    protected <T extends RealFieldElement<T>> void doBackward(final Field<T> field,
                                                              final double epsilonLast,
                                                              final double epsilonMaxValue,
                                                              final double epsilonMaxTime,
                                                              final String name) {

        TestFieldProblem5<T> pb = new TestFieldProblem5<T>(field);
        double range = pb.getFinalTime().subtract(pb.getInitialState().getTime()).getReal();

        AdamsFieldIntegrator<T> integ = createIntegrator(field, 4, 0, range, 1.0e-12, 1.0e-12);
        TestFieldProblemHandler<T> handler = new TestFieldProblemHandler<T>(pb, integ);
        integ.addStepHandler(handler);
        integ.integrate(new FieldExpandableODE<T>(pb), pb.getInitialState(), pb.getFinalTime());

        Assert.assertEquals(0.0, handler.getLastError().getReal(), epsilonLast);
        Assert.assertEquals(0.0, handler.getMaximalValueError().getReal(), epsilonMaxValue);
        Assert.assertEquals(0, handler.getMaximalTimeError().getReal(), epsilonMaxTime);
        Assert.assertEquals(name, integ.getName());
    }

    @Test
    public abstract void polynomial();

    protected <T extends RealFieldElement<T>> void doPolynomial(final Field<T> field,
                                                                final int nLimit,
                                                                final double epsilonBad,
                                                                final double epsilonGood) {
        TestFieldProblem6<T> pb = new TestFieldProblem6<T>(field);
        double range = pb.getFinalTime().subtract(pb.getInitialState().getTime()).abs().getReal();

        for (int nSteps = 2; nSteps < 8; ++nSteps) {
            AdamsFieldIntegrator<T> integ = createIntegrator(field, nSteps, 1.0e-6 * range, 0.1 * range, 1.0e-4, 1.0e-4);
            integ.setStarterIntegrator(new PerfectStarter<T>(pb, nSteps));
            TestFieldProblemHandler<T> handler = new TestFieldProblemHandler<T>(pb, integ);
            integ.addStepHandler(handler);
            integ.integrate(new FieldExpandableODE<T>(pb), pb.getInitialState(), pb.getFinalTime());
            if (nSteps < nLimit) {
                Assert.assertTrue(handler.getMaximalValueError().getReal() > epsilonBad);
            } else {
                Assert.assertTrue(handler.getMaximalValueError().getReal() < epsilonGood);
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

        FieldODEIntegrator<T> integrator = createIntegrator(field, 6, 0.001, 1.0, 1.0e-12, 1.0e-12);
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

    @Test(expected=MathIllegalStateException.class)
    public abstract void testStartFailure();

        protected <T extends RealFieldElement<T>> void doTestStartFailure(final Field<T> field) {
            TestFieldProblem1<T> pb = new TestFieldProblem1<T>(field);
        double minStep = pb.getFinalTime().subtract(pb.getInitialState().getTime()).multiply(0.0001).getReal();
        double maxStep = pb.getFinalTime().subtract(pb.getInitialState().getTime()).getReal();
        double scalAbsoluteTolerance = 1.0e-6;
        double scalRelativeTolerance = 1.0e-7;

        MultistepFieldIntegrator<T> integ = createIntegrator(field, 6, minStep, maxStep,
                                                             scalAbsoluteTolerance,
                                                             scalRelativeTolerance);
        integ.setStarterIntegrator(new DormandPrince853FieldIntegrator<T>(field, maxStep * 0.5, maxStep, 0.1, 0.1));
        TestFieldProblemHandler<T> handler = new TestFieldProblemHandler<T>(pb, integ);
        integ.addStepHandler(handler);
        integ.integrate(new FieldExpandableODE<T>(pb), pb.getInitialState(), pb.getFinalTime());

    }

    private static class PerfectStarter<T extends RealFieldElement<T>> extends AbstractFieldIntegrator<T> {

        private final PerfectInterpolator<T> interpolator;
        private final int nbSteps;

        public PerfectStarter(final TestFieldProblemAbstract<T> problem, final int nbSteps) {
            super(problem.getField(), "perfect-starter");
            this.interpolator = new PerfectInterpolator<T>(problem);
            this.nbSteps      = nbSteps;
        }

        public FieldODEStateAndDerivative<T> integrate(FieldExpandableODE<T> equations,
                                                       FieldODEState<T> initialState, T finalTime) {
            T tStart = initialState.getTime().add(finalTime.subtract(initialState.getTime()).multiply(0.01));
            getEvaluationsCounter().increment(nbSteps);
            interpolator.setCurrentTime(initialState.getTime());
            for (int i = 0; i < nbSteps; ++i) {
                T tK = initialState.getTime().multiply(nbSteps - 1 - (i + 1)).add(tStart.multiply(i + 1)).divide(nbSteps - 1);
                interpolator.setPreviousTime(interpolator.getCurrentTime());
                interpolator.setCurrentTime(tK);
                for (FieldODEStepHandler<T> handler : getStepHandlers()) {
                    handler.handleStep(interpolator, i == nbSteps - 1);
                }
            }
            return interpolator.getInterpolatedState(tStart);
        }

    }

    private static class PerfectInterpolator<T extends RealFieldElement<T>> implements FieldODEStateInterpolator<T> {
        private final TestFieldProblemAbstract<T> problem;
        private T previousTime;
        private T currentTime;

        public PerfectInterpolator(final TestFieldProblemAbstract<T> problem) {
            this.problem = problem;
        }

        public void setPreviousTime(T previousTime) {
            this.previousTime = previousTime;
        }

        public void setCurrentTime(T currentTime) {
            this.currentTime = currentTime;
        }

        public T getCurrentTime() {
            return currentTime;
        }

        public boolean isForward() {
            return problem.getFinalTime().subtract(problem.getInitialState().getTime()).getReal() >= 0;
        }

        public FieldODEStateAndDerivative<T> getPreviousState() {
            return getInterpolatedState(previousTime);
        }

        @Override
        public boolean isPreviousStateInterpolated() {
            return false;
        }

        public FieldODEStateAndDerivative<T> getCurrentState() {
            return getInterpolatedState(currentTime);
        }

        @Override
        public boolean isCurrentStateInterpolated() {
            return false;
        }

        public FieldODEStateAndDerivative<T> getInterpolatedState(T time) {
            T[] y    = problem.computeTheoreticalState(time);
            T[] yDot = problem.computeDerivatives(time, y);
            return new FieldODEStateAndDerivative<T>(time, y, yDot);
        }

    }

}
