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


import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.ode.ExpandableODE;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.sampling.AbstractODEStateInterpolator;

public class AdamsStateInterpolatorTest extends ODEStateInterpolatorAbstractTest {

    @Override
    protected AbstractODEStateInterpolator
        setUpInterpolator(ReferenceODE eqn, double t0, double[] y0, double t1) {
        final int        nSteps   = 12;
        final double     h        = (t1 - t0) / (nSteps - 1);
        final int        nbPoints = (nSteps + 3) / 2;
        final double[]   t        = new double[nbPoints];
        final double[][] y        = new double[nbPoints][];
        final double[][] yDot     = new double[nbPoints][];
        for (int i = 0; i < nbPoints; ++i) {
            t[i]    = t0 + i * h;
            y[i]    = eqn.theoreticalState(t[i]);
            yDot[i] = eqn.computeDerivatives(t[i], y[i]);
        }
        AdamsNordsieckTransformer transformer = AdamsNordsieckTransformer.getInstance(nSteps);
        Array2DRowRealMatrix      nordsieck   = transformer.initializeHighOrderDerivatives(h, t, y, yDot);

        double[] scaled = new double[eqn.getDimension()];
        for (int i = 0; i < scaled.length; ++i) {
            scaled[i] = h * yDot[0][i];
        }
        double   tCurrent    = t1;
        double[] yCurrent    = eqn.theoreticalState(tCurrent);
        double[] yDotCurrent = eqn.computeDerivatives(tCurrent, yCurrent);

        ODEStateAndDerivative previous = new ODEStateAndDerivative(t[0], y[0], yDot[0]);
        ODEStateAndDerivative current  = new ODEStateAndDerivative(tCurrent, yCurrent, yDotCurrent);
        return new AdamsStateInterpolator(h, previous, scaled, nordsieck, t1 >= t0,
                                          previous, current, new ExpandableODE(eqn).getMapper());

    }

    @Override
    public void interpolationAtBounds() {
        // as the Adams step interpolator is based on a Taylor expansion since previous step,
        // the maximum error is at step end and not in the middle of the step
        // as for Runge-Kutta integrators
        doInterpolationAtBounds(1.4e-6);
    }

    @Override
    public void interpolationInside() {
        doInterpolationInside(3.3e-10, 1.4e-6);
    }

    @Override
    public void restrictPrevious() {
        doRestrictPrevious(1.0e-15, 1.0e-15);
    }

    @Override
    public void restrictCurrent() {
        doRestrictCurrent(1.0e-15, 1.0e-15);
    }

    @Override
    public void restrictBothEnds() {
        doRestrictBothEnds(1.0e-15, 1.0e-15);
    }

    @Override
    public void degenerateInterpolation() {
        doDegenerateInterpolation();
    }

}
