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
package org.hipparchus.filtering.kalman;

import org.hipparchus.linear.RealMatrix;

/** Interface representing a Kalman estimate.
 * @since 4.0
 */
public interface KalmanEstimate {

    /** Get the current predicted state.
     * @return current predicted state
     */
    ProcessEstimate getPredicted();

    /** Get the current corrected state.
     * @return current corrected state
     */
    ProcessEstimate getCorrected();

    /** Get the cross-covariance between the previous state and the prediction.
     * Not required for forward filtering, but required for the smoother.
     * @return cross-covariance
     */
    RealMatrix getStateCrossCovariance();

}
