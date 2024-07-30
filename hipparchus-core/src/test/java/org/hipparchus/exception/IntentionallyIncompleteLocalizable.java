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
package org.hipparchus.exception;

import java.util.Locale;

/**
 * Intentionally incomplete localizable, for testing purposes
 * @author Luc Maisonobe
 */
enum IntentionallyIncompleteLocalizable implements Localizable
{
    MESSAGE_WITHOUT_ARGUMENT("message without argument"),

    MESSAGE_WITH_ONE_ARGUMENT("message with one argument {0}"),

    MESSAGE_WITH_TWO_ARGUMENTS("message with two arguments {0}, {1}");

    /**
     * Simple constructor.
     * @param sourceFormat source English format to use when no
     *                     localized version is available
     */
    IntentionallyIncompleteLocalizable(final String sourceFormat)
    {
        this.sourceFormat = sourceFormat;
    }

    /** Source English format. */
    private final String sourceFormat;

    /** {@inheritDoc} */
    @Override
    public String getSourceString()
    {
        return sourceFormat;
    }

    /** {@inheritDoc} */
    @Override
    public String getLocalizedString(final Locale locale)
    {
        return getLocalizedString("assets/" + IntentionallyIncompleteLocalizable.class.getName().replaceAll("\\.", "/"),
                                  name(), locale);
    }

}
