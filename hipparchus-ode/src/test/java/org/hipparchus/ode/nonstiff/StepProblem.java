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

import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.solvers.BracketedUnivariateSolver;
import org.hipparchus.analysis.solvers.BracketingNthOrderBrentSolver;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.ode.events.AbstractODEDetector;
import org.hipparchus.ode.events.Action;
import org.hipparchus.ode.events.AdaptableInterval;
import org.hipparchus.ode.events.ODEEventDetector;
import org.hipparchus.ode.events.ODEEventHandler;


public class StepProblem extends AbstractODEDetector<StepProblem> implements OrdinaryDifferentialEquation {

    private double rateBefore;
    private double rateAfter;
    private double rate;
    private double switchTime;

    public StepProblem(final AdaptableInterval maxCheck, final double threshold, final int maxIter,
                       double rateBefore, double rateAfter, double switchTime) {
        this(maxCheck, maxIter, new BracketingNthOrderBrentSolver(0, threshold, 0, 5),
             new LocalHandler(), rateBefore, rateAfter, switchTime);
    }

    private StepProblem(final AdaptableInterval maxCheck, final int maxIter,
                        final BracketedUnivariateSolver<UnivariateFunction> solver,
                        final ODEEventHandler handler,
                        final double rateBefore, final double rateAfter,
                        final double switchTime) {
        super(maxCheck, maxIter, solver, handler);
        this.rateBefore = rateBefore;
        this.rateAfter  = rateAfter;
        this.switchTime = switchTime;
        setRate(rateBefore);
    }

    protected StepProblem create(AdaptableInterval newMaxCheck, int newMaxIter,
                                 BracketedUnivariateSolver<UnivariateFunction> newSolver,
                                 ODEEventHandler newHandler) {
        return new StepProblem(newMaxCheck, newMaxIter, newSolver, newHandler,
                               rateBefore, rateAfter, switchTime);
    }

    public double[] computeDerivatives(double t, double[] y) {
        return new double[] {
            rate
        };
    }

    public int getDimension() {
        return 1;
    }

    public void setRate(double rate) {
        this.rate = rate;
    }

    public void init(double t0, double[] y0, double t) {
    }

    private static class LocalHandler implements ODEEventHandler {
        public Action eventOccurred(ODEStateAndDerivative state, ODEEventDetector detector, boolean increasing) {
            final StepProblem sp = (StepProblem) detector;
            sp.setRate(sp.rateAfter);
            return Action.RESET_DERIVATIVES;
        }
    }

    public double g(ODEStateAndDerivative s) {
        return s.getTime() - switchTime;
    }

}
