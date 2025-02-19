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
package org.hipparchus.optim;

/**
 * Base class for all convergence checker implementations.
 *
 * @param <P> Type of (point, value) pair.
 *
 */
public abstract class AbstractConvergenceChecker<P>
    implements ConvergenceChecker<P> {
    /**
     * Relative tolerance threshold.
     */
    private final double relativeThreshold;
    /**
     * Absolute tolerance threshold.
     */
    private final double absoluteThreshold;

    /**
     * Build an instance with a specified thresholds.
     *
     * @param relativeThreshold relative tolerance threshold
     * @param absoluteThreshold absolute tolerance threshold
     */
    protected AbstractConvergenceChecker(final double relativeThreshold,
                                         final double absoluteThreshold) {
        this.relativeThreshold = relativeThreshold;
        this.absoluteThreshold = absoluteThreshold;
    }

    /** Get relative threshold.
     * @return the relative threshold.
     */
    public double getRelativeThreshold() {
        return relativeThreshold;
    }

    /** Get absolute threshold.
     * @return the absolute threshold.
     */
    public double getAbsoluteThreshold() {
        return absoluteThreshold;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public abstract boolean converged(int iteration,
                                      P previous,
                                      P current);
}
