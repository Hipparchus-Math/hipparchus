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

import org.hipparchus.analysis.UnivariateFunction;
import org.hipparchus.analysis.solvers.BracketedUnivariateSolver;
import org.hipparchus.analysis.solvers.BracketingNthOrderBrentSolver;
import org.hipparchus.ode.events.Action;
import org.hipparchus.ode.events.AdaptableInterval;
import org.hipparchus.ode.events.ODEEventDetector;
import org.hipparchus.ode.events.ODEEventHandler;
import org.hipparchus.util.FastMath;

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

 */
public class TestProblem4 extends TestProblemAbstract {

    private static final double OFFSET = 1.2;

    /** Time offset. */
    private double a;

    /** Simple constructor. */
    public TestProblem4() {
        super(0.0, new double[] { FastMath.sin(OFFSET), FastMath.cos(OFFSET) },
              15, new double[] { 1.0, 0.0 });
        a = OFFSET;
    }

    @Override
    public ODEEventDetector[] getEventDetectors(final double maxCheck, final double threshold, final int maxIter) {
        return new ODEEventDetector[] {
            new Bounce(maxCheck, threshold, maxIter),
            new Stop(maxCheck, threshold, maxIter)
        };
    }

    /**
     * Get the theoretical events times.
     * @return theoretical events times
     */
    @Override
    public double[] getTheoreticalEventsTimes() {
        return new double[] {
                             1 * FastMath.PI - a,
                             2 * FastMath.PI - a,
                             3 * FastMath.PI - a,
                             4 * FastMath.PI - a,
                             12.0
        };
    }

    @Override
    public double[] doComputeDerivatives(double t, double[] y) {
        return new double[] {
          y[1], -y[0]
        };
    }

    @Override
    public double[] computeTheoreticalState(double t) {
        double sin = FastMath.sin(t + a);
        double cos = FastMath.cos(t + a);
        return new double[] {
            FastMath.abs(sin),
            (sin >= 0) ? cos : -cos
        };
    }

    private static abstract class BaseDetector implements ODEEventDetector, ODEEventHandler {

        final double                                        maxCheck;
        final int                                           maxIter;
        final BracketedUnivariateSolver<UnivariateFunction> solver;

        protected BaseDetector(final double maxCheck, final double threshold, final int maxIter) {
            this.maxCheck  = maxCheck;
            this.maxIter   = maxIter;
            this.solver    = new BracketingNthOrderBrentSolver(0, threshold, 0, 5);
        }

        public AdaptableInterval getMaxCheckInterval() {
            return s -> maxCheck;
        }

        public int getMaxIterationCount() {
            return maxIter;
        }

        public BracketedUnivariateSolver<UnivariateFunction> getSolver() {
            return solver;
        }

        public ODEEventHandler getHandler() {
            return this;
        }

    }

    private static class Bounce extends BaseDetector {

        private int sign;

        public Bounce(final double maxCheck, final double threshold, final int maxIter) {
            super(maxCheck, threshold, maxIter);
            sign = +1;
        }

        public double g(ODEStateAndDerivative s) {
            return sign * s.getPrimaryState()[0];
        }

        public Action eventOccurred(ODEStateAndDerivative s, ODEEventDetector detector, boolean increasing) {
            // this sign change is needed because the state will be reset soon
            sign = -sign;
            return Action.RESET_STATE;
        }

        public ODEStateAndDerivative resetState(ODEEventDetector detector, ODEStateAndDerivative s) {
            final double[] y    = s.getPrimaryState();
            final double[] yDot = s.getPrimaryDerivative();
            y[0]    = -y[0];
            y[1]    = -y[1];
            yDot[0] = -yDot[0];
            yDot[1] = -yDot[1];
            return new ODEStateAndDerivative(s.getTime(), y, yDot);
        }

    }

    private static class Stop extends BaseDetector {

        public Stop(final double maxCheck, final double threshold, final int maxIter) {
            super(maxCheck, threshold, maxIter);
        }

        public double g(ODEStateAndDerivative s) {
            return s.getTime() - 12.0;
        }

        public Action eventOccurred(ODEStateAndDerivative s, ODEEventDetector detector, boolean increasing) {
            return Action.STOP;
        }

    }

}
