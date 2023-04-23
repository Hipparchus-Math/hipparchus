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
package org.hipparchus.optim;

import java.util.List;

/** Multiplexer for {@link ConvergenceChecker}, checking <em>all</em> the checkers converged.
 * <p>
 * The checkers are checked in the order of the initial list and the check loop
 * is interrupted as soon as one checker fails to converge (that is the remaining
 * checkers may <em>not</em> be called in first iterations.
 * </p>
 * @param <P> type of the evaluation
 * @since 2.1
 */
public class ConvergenceCheckerAndMultiplexer<P> implements ConvergenceChecker<P> {

    /** Underlying checkers. */
    private final List<ConvergenceChecker<P>> checkers;

    /** Simple constructor.
     * @param checkers checkers to use, convergence is reached when
     * <em>all</em> checkers have converged
     */
    public ConvergenceCheckerAndMultiplexer(final List<ConvergenceChecker<P>> checkers) {
        this.checkers = checkers;
    }

    /** {@inheritDoc} */
    @Override
    public boolean converged(final int iteration, final P previous, final P current) {
        for (final ConvergenceChecker<P> checker : checkers) {
            if (!checker.converged(iteration, previous, current)) {
                return false;
            }
        }
        return true;
    }

}
