/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */

package org.hipparchus.migration.ode.sampling;

import org.hipparchus.ode.ODEIntegrator;
import org.junit.Assert;


@Deprecated
public class DeprecatedStepHandler implements StepHandler {

    private final ODEIntegrator integrator;
    private boolean initCalled   = false;
    private boolean lastStepSeen = false;

    public DeprecatedStepHandler(final ODEIntegrator integrator) {
        this.integrator = integrator;
    }

    public void init(double t0, double[] y0, double t) {
        initCalled = true;
    }

    public void handleStep(MigrationStepInterpolator interpolator, boolean isLast) {
        Assert.assertEquals(interpolator.getPreviousTime(),
                            integrator.getCurrentStepStart(),
                            1.0e-10);
        if (isLast) {
            lastStepSeen = true;
        }
    }

    public boolean isInitCalled() {
        return initCalled;
    }

    public boolean isLastStepSeen() {
        return lastStepSeen;
    }

}
