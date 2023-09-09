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

package org.hipparchus.ode;

import java.lang.reflect.Array;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.Field;
import org.hipparchus.analysis.solvers.BracketedRealFieldUnivariateSolver;
import org.hipparchus.analysis.solvers.FieldBracketingNthOrderBrentSolver;
import org.hipparchus.ode.events.Action;
import org.hipparchus.ode.events.FieldAdaptableInterval;
import org.hipparchus.ode.events.FieldODEEventDetector;
import org.hipparchus.ode.events.FieldODEEventHandler;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;

/**
 * This class is used in the junit tests for the ODE integrators.

 * <p>This specific problem is the following differential equation :
 * <pre>
 *    x'' = -x
 * </pre>
 * And when x decreases down to 0, the state should be changed as follows :
 * <pre>
 *   x' → -x'
 * </pre>
 * The theoretical solution of this problem is x = |sin(t+a)|
 * </p>

 * @param <T> the type of the field elements
 */
public class TestFieldProblem4<T extends CalculusFieldElement<T>>
    extends TestFieldProblemAbstract<T> {

    private static final double OFFSET = 1.2;

    /** Time offset. */
    private T a;

    /** Simple constructor.
     * @param field field to which elements belong
     */
    public TestFieldProblem4(Field<T> field) {
        super(convert(field, 0.0),
              createY0(field),
              convert(field, 15),
              convert(field, 1.0, 0.0));
        a = convert(field, OFFSET);
    }

    private static <T extends CalculusFieldElement<T>> T[] createY0(final Field<T> field) {
        final T a = convert(field, OFFSET);
        T[] y0 = MathArrays.buildArray(field, 2);
        y0[0] = a.sin();
        y0[1] = a.cos();
        return y0;
    }

    @Override
    public FieldODEEventDetector<T>[] getEventDetectors(final double maxCheck, final T threshold, final int maxIter) {
        @SuppressWarnings("unchecked")
        FieldODEEventDetector<T>[] handlers =
                        (FieldODEEventDetector<T>[]) Array.newInstance(FieldODEEventDetector.class, 2);
        handlers[0] = new Bounce<>(maxCheck, threshold, maxIter);
        handlers[1] = new Stop<>(maxCheck, threshold, maxIter);
        return handlers;
    }

    /**
     * Get the theoretical events times.
     * @return theoretical events times
     */
    @Override
    public T[] getTheoreticalEventsTimes() {
        T[] array = MathArrays.buildArray(getField(), 5);
        array[0] = a.negate().add(1 * FastMath.PI);
        array[1] = a.negate().add(2 * FastMath.PI);
        array[2] = a.negate().add(3 * FastMath.PI);
        array[3] = a.negate().add(4 * FastMath.PI);
        array[4] = convert(a.getField(), 120.0);
        return array;
    }

    @Override
    public T[] doComputeDerivatives(T t, T[] y) {
        final T[] yDot = MathArrays.buildArray(getField(), getDimension());
        yDot[0] = y[1];
        yDot[1] = y[0].negate();
        return yDot;
    }

    @Override
    public T[] computeTheoreticalState(T t) {
        T sin = t.add(a).sin();
        T cos = t.add(a).cos();
        final T[] y = MathArrays.buildArray(getField(), getDimension());
        y[0] = sin.abs();
        y[1] = (sin.getReal() >= 0) ? cos : cos.negate();
        return y;
    }

    private static abstract class BaseDetector<T extends CalculusFieldElement<T>>
        implements FieldODEEventDetector<T>, FieldODEEventHandler<T> {

        final double                                maxCheck;
        final int                                   maxIter;
        final BracketedRealFieldUnivariateSolver<T> solver;

        protected BaseDetector(final double maxCheck, final T threshold, final int maxIter) {
            this.maxCheck  = maxCheck;
            this.maxIter   = maxIter;
            this.solver    = new FieldBracketingNthOrderBrentSolver<>(threshold.getField().getZero(),
                                                                      threshold,
                                                                      threshold.getField().getZero(),
                                                                      5);
        }

        public FieldAdaptableInterval<T> getMaxCheckInterval() {
            return s -> maxCheck;
        }

        public int getMaxIterationCount() {
            return maxIter;
        }

        public BracketedRealFieldUnivariateSolver<T> getSolver() {
            return solver;
        }

        public FieldODEEventHandler<T> getHandler() {
            return this;
        }

    }

    private static class Bounce<T extends CalculusFieldElement<T>> extends BaseDetector<T> {

        private int sign;

        public Bounce(final double maxCheck, final T threshold, final int maxIter) {
            super(maxCheck, threshold, maxIter);
            sign = +1;
        }

        public void init(FieldODEStateAndDerivative<T> state0, T t) {
        }

        public T g(FieldODEStateAndDerivative<T> state) {
            return state.getPrimaryState()[0].multiply(sign);
        }

        public Action eventOccurred(FieldODEStateAndDerivative<T> state,
                                    FieldODEEventDetector<T> detector,
                                    boolean increasing) {
            // this sign change is needed because the state will be reset soon
            sign = -sign;
            return Action.RESET_STATE;
        }

        public FieldODEState<T> resetState(FieldODEEventDetector<T> detector,
                                           FieldODEStateAndDerivative<T> state) {
            T[] y = state.getPrimaryState();
            y[0] = y[0].negate();
            y[1] = y[1].negate();
            return new FieldODEState<T>(state.getTime(), y);
        }

    }

    private static class Stop<T extends CalculusFieldElement<T>> extends BaseDetector<T> {

        public Stop(final double maxCheck, final T threshold, final int maxIter) {
            super(maxCheck, threshold, maxIter);
        }

        public void init(FieldODEStateAndDerivative<T> state0, T t) {
        }

        public T g(FieldODEStateAndDerivative<T> state) {
            return state.getTime().subtract(12.0);
        }

        public Action eventOccurred(FieldODEStateAndDerivative<T> state,
                                    FieldODEEventDetector<T> detector,
                                    boolean increasing) {
            return Action.STOP;
        }

        public FieldODEState<T> resetState(FieldODEEventDetector<T> detector,
                                           FieldODEStateAndDerivative<T> state) {
            return state;
        }

    }

}
