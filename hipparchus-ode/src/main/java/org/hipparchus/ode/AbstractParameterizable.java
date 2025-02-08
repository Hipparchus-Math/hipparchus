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
package org.hipparchus.ode;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import org.hipparchus.exception.MathIllegalArgumentException;

/** This abstract class provides boilerplate parameters list.
 *
 */

public abstract class AbstractParameterizable implements Parameterizable {

   /** List of the parameters names. */
    private final List<String> parametersNames;

    /** Simple constructor.
     * @param names names of the supported parameters
     */
    protected AbstractParameterizable(final String ... names) {
        parametersNames = new ArrayList<>();
        Collections.addAll(parametersNames, names);
    }

    /** Simple constructor.
     * @param names names of the supported parameters
     */
    protected AbstractParameterizable(final Collection<String> names) {
        parametersNames = new ArrayList<>();
        parametersNames.addAll(names);
    }

    /** {@inheritDoc} */
    @Override
    public List<String> getParametersNames() {
        return parametersNames;
    }

    /** {@inheritDoc} */
    @Override
    public boolean isSupported(final String name) {
        for (final String supportedName : parametersNames) {
            if (supportedName.equals(name)) {
                return true;
            }
        }
        return false;
    }

    /** Check if a parameter is supported and throw an IllegalArgumentException if not.
     * @param name name of the parameter to check
     * @exception MathIllegalArgumentException if the parameter is not supported
     * @see #isSupported(String)
     */
    public void complainIfNotSupported(final String name)
        throws MathIllegalArgumentException {
        if (!isSupported(name)) {
            throw new MathIllegalArgumentException(LocalizedODEFormats.UNKNOWN_PARAMETER, name);
        }
    }

}
