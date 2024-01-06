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
package org.hipparchus.optim.nonlinear.vector.constrained;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.hipparchus.optim.OptimizationData;

public class InequalityConstraintSet implements OptimizationData {

    /** Constraints. */
    private final List<TwiceDifferentiableFunction> constraints;

    /**
     * Construct an inequality constraint set from a Collection of convex constraint functions
     * @param constraints the Collection of constraint functions to apply
     */
    public InequalityConstraintSet(Collection<TwiceDifferentiableFunction> constraints) {
        this.constraints = new ArrayList<>(constraints);
    }

    /**
     * Construct an inequality constraint set from an argument list (or array) of convex
     * constraint functions
     * @param constraints the list of constraint functions to apply
     */
    public InequalityConstraintSet(TwiceDifferentiableFunction... constraints) {
        this.constraints = new ArrayList<>(constraints.length);
        for (TwiceDifferentiableFunction f: constraints) {
            this.constraints.add(f);
        }
    }
}
