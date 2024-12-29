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

/** Observer for Kalman filter recursions.
 * <p>
 * This interface is intended to be implemented by users to monitor
 * the progress of the Kalman filter estimator during estimation.
 * </p>
 */
public interface KalmanObserver {

    /** Callback for initialisation of observer.
     * @param estimate estimate calculated by a Kalman filter
     */
    default void init(KalmanEstimate estimate) {}

    /** Notification callback after each Kalman filter measurement update.
     * @param estimate estimate calculated by a Kalman filter
     */
    void updatePerformed(KalmanEstimate estimate);
}
