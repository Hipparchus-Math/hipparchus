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
package org.hipparchus.migration.linear;


/**
 * Exception to be thrown when a positive definite matrix is expected.
 *
 * @deprecated as of 1.0, this exception is replaced by {@link org.hipparchus.exception.MathIllegalArgumentException}
 */
@Deprecated
public class NonPositiveDefiniteMatrixException
    extends org.hipparchus.migration.exception.NumberIsTooSmallException {
    /** Serializable version Id. */
    private static final long serialVersionUID = 1641613838113738061L;
    /** Index (diagonal element). */
    private final int index;
    /** Threshold. */
    private final double threshold;

    /**
     * Construct an exception.
     *
     * @param wrong Value that fails the positivity check.
     * @param index Row (and column) index.
     * @param threshold Absolute positivity threshold.
     */
    public NonPositiveDefiniteMatrixException(double wrong,
                                              int index,
                                              double threshold) {
        super(wrong, threshold, false);
        this.index = index;
        this.threshold = threshold;

    }

    /** Get row index.
     * @return the row index
     */
    public int getRow() {
        return index;
    }

    /** Get column index.
     * @return the column index
     */
    public int getColumn() {
        return index;
    }

    /** Get absolute positivity threshold.
     * @return the absolute positivity threshold
     */
    public double getThreshold() {
        return threshold;
    }

}
