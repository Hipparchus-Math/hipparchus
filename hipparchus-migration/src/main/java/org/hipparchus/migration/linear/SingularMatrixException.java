/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
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
package org.hipparchus.migration.linear;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;

/**
 * Exception to be thrown when a non-singular matrix is expected.
 *
 * @deprecated as of 1.0, this exception is replaced by {@link MathIllegalStateException}
 */
@Deprecated
public class SingularMatrixException extends MathIllegalArgumentException {
    /** Serializable version Id. */
    private static final long serialVersionUID = -4206514844735401070L;

    /**
     * Construct an exception.
     */
    public SingularMatrixException() {
        super(org.hipparchus.migration.exception.util.LocalizedFormats.SINGULAR_MATRIX);
    }
}
