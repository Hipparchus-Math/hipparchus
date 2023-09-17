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
package org.hipparchus.migration.ode;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.migration.exception.util.LocalizedFormats;

/**
 * Exception to be thrown when a parameter is unknown.
 *
 * @deprecated as of 1.0, this exception is replaced by {@link MathIllegalArgumentException}
 */
@Deprecated
public class UnknownParameterException extends MathIllegalArgumentException {

    /** Serializable version Id. */
    private static final long serialVersionUID = 20120902L;

    /** Parameter name. */
    private final String name;

    /**
     * Construct an exception from the unknown parameter.
     *
     * @param name parameter name.
     */
    public UnknownParameterException(final String name) {
        super(LocalizedFormats.UNKNOWN_PARAMETER, name);
        this.name = name;
    }

    /** Get name of the unknown parameter.
     * @return the name of the unknown parameter
     */
    public String getName() {
        return name;
    }

}
