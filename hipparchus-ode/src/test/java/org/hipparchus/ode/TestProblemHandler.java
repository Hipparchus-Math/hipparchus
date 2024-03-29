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

import org.hipparchus.ode.sampling.ODEStateInterpolator;
import org.hipparchus.ode.sampling.ODEStepHandler;
import org.hipparchus.util.FastMath;

/**
 * This class is used to handle steps for the test problems
 * integrated during the junit tests for the ODE integrators.
 */
public class TestProblemHandler implements ODEStepHandler {

    /** Associated problem. */
    private TestProblemAbstract problem;

    /** Maximal errors encountered during the integration. */
    private double maxValueError;
    private double maxTimeError;

    /** Error at the end of the integration. */
    private double lastError;

    /** Time at the end of integration. */
    private double lastTime;

    /** ODE solver used. */
    private ODEIntegrator integrator;

    /** Expected start for step. */
    private double expectedStepStart;

    /**
     * Simple constructor.
     * @param problem problem for which steps should be handled
     * @param integrator ODE solver used
     */
    public TestProblemHandler(TestProblemAbstract problem, ODEIntegrator integrator) {
        this.problem = problem;
        this.integrator = integrator;
        maxValueError = 0;
        maxTimeError  = 0;
        lastError     = 0;
        expectedStepStart = Double.NaN;
    }

    public void init(ODEStateAndDerivative s0, double t) {
        maxValueError = 0;
        maxTimeError  = 0;
        lastError     = 0;
        expectedStepStart = Double.NaN;
    }

    public void handleStep(ODEStateInterpolator interpolator) {

        double start = interpolator.getPreviousState().getTime();
        if (FastMath.abs((start - problem.getInitialTime()) / integrator.getCurrentSignedStepsize()) > 0.001) {
            // multistep integrators do not handle the first steps themselves
            // so we have to make sure the integrator we look at has really started its work
            if (!Double.isNaN(expectedStepStart)) {
                // the step should either start at the end of the integrator step
                // or at an event if the step is split into several substeps
                double stepError = FastMath.max(maxTimeError, FastMath.abs(start - expectedStepStart));
                for (double eventTime : problem.getTheoreticalEventsTimes()) {
                    stepError = FastMath.min(stepError, FastMath.abs(start - eventTime));
                }
                maxTimeError = FastMath.max(maxTimeError, stepError);
            }
            expectedStepStart = interpolator.getCurrentState().getTime();
        }


        double pT = interpolator.getPreviousState().getTime();
        double cT = interpolator.getCurrentState().getTime();
        double[] errorScale = problem.getErrorScale();

        // walk through the step
        for (int k = 0; k <= 20; ++k) {

            double time = pT + (k * (cT - pT)) / 20;
            ODEStateAndDerivative interpolated = interpolator.getInterpolatedState(time);
            double[] interpolatedY = interpolated.getPrimaryState();
            double[] theoreticalY  = problem.computeTheoreticalState(interpolated.getTime());

            // update the errors
            for (int i = 0; i < interpolatedY.length; ++i) {
                double error = errorScale[i] * FastMath.abs(interpolatedY[i] - theoreticalY[i]);
                maxValueError = FastMath.max(error, maxValueError);
            }
        }

    }

    public void finish(ODEStateAndDerivative finalState) {
        double[] theoreticalY  = problem.computeTheoreticalState(finalState.getTime());
        for (int i = 0; i < finalState.getCompleteState().length; ++i) {
            double error = FastMath.abs(finalState.getCompleteState()[i] - theoreticalY[i]);
            lastError = FastMath.max(error, lastError);
        }
        lastTime = finalState.getTime();
    }

    /**
     * Get the maximal value error encountered during integration.
     * @return maximal value error
     */
    public double getMaximalValueError() {
        return maxValueError;
    }

    /**
     * Get the maximal time error encountered during integration.
     * @return maximal time error
     */
    public double getMaximalTimeError() {
        return maxTimeError;
    }


    public int getCalls() {
        return problem.getCalls();
    }

    /**
     * Get the error at the end of the integration.
     * @return error at the end of the integration
     */
    public double getLastError() {
        return lastError;
    }

    /**
     * Get the time at the end of the integration.
     * @return time at the end of the integration.
     */
    public double getLastTime() {
        return lastTime;
    }

}
