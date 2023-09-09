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
import org.hipparchus.analysis.solvers.BracketedRealFieldUnivariateSolver;
import org.hipparchus.analysis.solvers.FieldBracketingNthOrderBrentSolver;
import org.hipparchus.ode.FieldODEStateAndDerivative;
import org.hipparchus.ode.FieldOrdinaryDifferentialEquation;
import org.hipparchus.ode.events.AbstractFieldODEDetector;
import org.hipparchus.ode.events.Action;
import org.hipparchus.ode.events.FieldAdaptableInterval;
import org.hipparchus.ode.events.FieldODEEventDetector;
import org.hipparchus.ode.events.FieldODEEventHandler;
import org.hipparchus.util.MathArrays;


public class StepFieldProblem<T extends CalculusFieldElement<T>>
    extends AbstractFieldODEDetector<StepFieldProblem<T>, T>
    implements FieldOrdinaryDifferentialEquation<T> {

    private Field<T> field;
    private T        rateBefore;
    private T        rateAfter;
    private T        rate;
    private T        switchTime;

    public StepFieldProblem(Field<T> field,
                            final FieldAdaptableInterval<T> maxCheck, final T threshold, final int maxIter,
                            T rateBefore, T rateAfter, T switchTime) {
        this(field, maxCheck, maxIter,
             new FieldBracketingNthOrderBrentSolver<>(field.getZero(), threshold, field.getZero(), 5),
             new LocalHandler<>(),
             rateBefore, rateAfter, switchTime);
    }

    private StepFieldProblem(Field<T> field,
                             final FieldAdaptableInterval<T> maxCheck, final int maxIter,
                             final BracketedRealFieldUnivariateSolver<T> solver,
                             final FieldODEEventHandler<T> handler,
                             final T rateBefore, final T rateAfter,
                             final T switchTime) {
        super(maxCheck, maxIter, solver, handler);
        this.field      = field;
        this.rateBefore = rateBefore;
        this.rateAfter  = rateAfter;
        this.switchTime = switchTime;
        setRate(rateBefore);
    }

    protected StepFieldProblem<T> create(FieldAdaptableInterval<T> newMaxCheck, int newMaxIter,
                                         BracketedRealFieldUnivariateSolver<T> newSolver,
                                         FieldODEEventHandler<T> newHandler) {
        return new StepFieldProblem<>(field, newMaxCheck, newMaxIter, newSolver, newHandler,
                                      rateBefore, rateAfter, switchTime);
    }

    public T[] computeDerivatives(T t, T[] y) {
        T[] yDot = MathArrays.buildArray(field, 1);
        yDot[0] = rate;
        return yDot;
    }

    public int getDimension() {
        return 1;
    }

    public void setRate(T rate) {
        this.rate = rate;
    }

    public void init(T t0, T[] y0, T t) {
    }

    public void init(FieldODEStateAndDerivative<T> state0, T t) {
    }

    private static class LocalHandler<T extends CalculusFieldElement<T>> implements FieldODEEventHandler<T> {
        public Action eventOccurred(FieldODEStateAndDerivative<T> state, FieldODEEventDetector<T> detector, boolean increasing) {
            final StepFieldProblem<T> sp = (StepFieldProblem<T>) detector;
            sp.setRate(sp.rateAfter);
            return Action.RESET_DERIVATIVES;
        }

    }

    public T g(FieldODEStateAndDerivative<T> state) {
        return state.getTime().subtract(switchTime);
    }

}
