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
package org.hipparchus.migration.exception;

import org.hipparchus.exception.Localizable;
import org.hipparchus.exception.LocalizedCoreFormats;

/**
 * Exception to be thrown when two dimensions differ.
 *
 * @deprecated as of 1.0, this exception is replaced by {@link org.hipparchus.exception.MathIllegalArgumentException}
 */
@Deprecated
public class DimensionMismatchException extends MathIllegalNumberException {
    /** Serializable version Id. */
    private static final long serialVersionUID = -8415396756375798143L;
    /** Correct dimension. */
    private final int dimension;

    /**
     * Construct an exception from the mismatched dimensions.
     *
     * @param specific Specific context information pattern.
     * @param wrong Wrong dimension.
     * @param expected Expected dimension.
     */
    public DimensionMismatchException(Localizable specific,
                                      int wrong,
                                      int expected) {
        super(specific, Integer.valueOf(wrong), Integer.valueOf(expected));
        dimension = expected;
    }

    /**
     * Construct an exception from the mismatched dimensions.
     *
     * @param wrong Wrong dimension.
     * @param expected Expected dimension.
     */
    public DimensionMismatchException(int wrong,
                                      int expected) {
        this(LocalizedCoreFormats.DIMENSIONS_MISMATCH, wrong, expected);
    }

    /** Get expected dimension.
     * @return the expected dimension.
     */
    public int getDimension() {
        return dimension;
    }
}
