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
 * An exception to be thrown when the condition number of a
 * {@link org.hipparchus.linear.RealLinearOperator} is too high.
 *
 * @deprecated as of 1.0, this exception is replaced by {@link MathIllegalArgumentException}
 */
@Deprecated
public class IllConditionedOperatorException
    extends MathIllegalArgumentException {
    /** Serializable version Id. */
    private static final long serialVersionUID = -7883263944530490135L;

    /**
     * Creates a new instance of this class.
     *
     * @param cond An estimate of the condition number of the offending linear
     * operator.
     */
    public IllConditionedOperatorException(final double cond) {
        super(org.hipparchus.migration.exception.util.LocalizedFormats.ILL_CONDITIONED_OPERATOR, cond);
    }
}
