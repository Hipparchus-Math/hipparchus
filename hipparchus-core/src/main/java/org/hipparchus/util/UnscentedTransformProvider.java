/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.util;

import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;

/**
 * Provider for unscented transform.
 */
public interface UnscentedTransformProvider {

    /**
     * Perform the unscented transform from a state and its covariance.
     * @param state state
     * @param covariance covariance associated with the state
     * @return an array of realvector containing the sigma points
     */
    RealVector[] unscentedTransform(RealVector state, RealMatrix covariance);

    /**
     * Computes mean from samples.
     * @param samples
     * @return mean
     */
    RealVector getMean(RealVector[] samples);

    /**
     * Computes covariance from state and samples.
     * @param samples
     * @param state
     * @return covariance matrix
     */
    RealMatrix getCovariance(RealVector[] samples, RealVector state);

    /**
     * Computes cross covariance from two states and two sets of samples.
     * @param firstSamplesSet first samples set
     * @param firstState first state
     * @param secondSamplesSet seconnd samples set
     * @param secondState second state
     * @return cross covariance matrix
     */
    RealMatrix getCrossCovariance(RealVector[] firstSamplesSet, RealVector firstState, RealVector[] secondSamplesSet, RealVector secondState);


}
