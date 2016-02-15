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
package org.hipparchus.exception;

import org.hipparchus.exception.util.LocalizedFormats;

/**
 * Exception to be thrown when the maximal number of iterations is exceeded.
 *
 * @since 3.1
 */
public class TooManyIterationsException extends MaxCountExceededException {
    /** Serializable version Id. */
    private static final long serialVersionUID = 20121211L;

    /**
     * Construct the exception.
     *
     * @param max Maximum number of evaluations.
     */
    public TooManyIterationsException(Number max) {
        super(max);
        getContext().addMessage(LocalizedFormats.ITERATIONS);
    }
}
