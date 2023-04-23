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

import org.hipparchus.exception.MathIllegalArgumentException;

/**
 * Exception to be thrown when a symmetric, definite positive
 * {@link org.hipparchus.linear.RealLinearOperator} is expected.
 * Since the coefficients of the matrix are not accessible, the most
 * general definition is used to check that {@code A} is not positive
 * definite, i.e.  there exists {@code x} such that {@code x' A x <= 0}.
 * In the terminology of this exception, {@code A} is the "offending"
 * linear operator and {@code x} the "offending" vector.
 *
 * @deprecated as of 1.0, this exception is replaced by {@link MathIllegalArgumentException}
 */
@Deprecated
public class NonPositiveDefiniteOperatorException
    extends MathIllegalArgumentException {
    /** Serializable version Id. */
    private static final long serialVersionUID = 917034489420549847L;

    /** Creates a new instance of this class. */
    public NonPositiveDefiniteOperatorException() {
        super(org.hipparchus.migration.exception.util.LocalizedFormats.NON_POSITIVE_DEFINITE_OPERATOR);
    }
}
