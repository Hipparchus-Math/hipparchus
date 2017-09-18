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

package org.hipparchus.migration.ode.events;

@Deprecated
public class DeprecatedEventHandler implements EventHandler {

    private boolean initCalled = false;
    private boolean resetCalled = false;

    public void init(double t0, double[] y0, double t) {
        initCalled = true;
    }

    public double g(double t, double[] y) {
        return (t - 2.0) * (t - 4.0);
    }

    public Action eventOccurred(double t, double[] y, boolean increasing) {
        return t < 3 ? Action.RESET_STATE : Action.STOP;
    }

    public void resetState(double t, double[] y) {
        resetCalled = true;
    }

    public boolean isInitCalled() {
        return initCalled;
    }

    public boolean isResetCalled() {
        return resetCalled;
    }

}
