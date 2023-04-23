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
package org.hipparchus.ode;

import java.util.Locale;
import java.util.MissingResourceException;
import java.util.ResourceBundle;

import org.hipparchus.exception.Localizable;
import org.hipparchus.exception.UTF8Control;

/**
 * Enumeration for localized messages formats used in exceptions messages.
 * <p>
 * The constants in this enumeration represent the available
 * formats as localized strings. These formats are intended to be
 * localized using simple properties files, using the constant
 * name as the key and the property value as the message format.
 * The source English format is provided in the constants themselves
 * to serve both as a reminder for developers to understand the parameters
 * needed by each format, as a basis for translators to create
 * localized properties files, and as a default format if some
 * translation is missing.
 * </p>
 */
public enum LocalizedODEFormats implements Localizable {

    /** HOLE_BETWEEN_MODELS_TIME_RANGES. */
    HOLE_BETWEEN_MODELS_TIME_RANGES("{0} wide hole between models time ranges"),

    /** INTEGRATION_METHOD_NEEDS_AT_LEAST_TWO_PREVIOUS_POINTS. */
    INTEGRATION_METHOD_NEEDS_AT_LEAST_TWO_PREVIOUS_POINTS("multistep method needs at least {0} previous steps, got {1}"),

    /** MINIMAL_STEPSIZE_REACHED_DURING_INTEGRATION. */
    MINIMAL_STEPSIZE_REACHED_DURING_INTEGRATION("minimal step size ({1,number,0.00E00}) reached, integration needs {0,number,0.00E00}"),

    /** MULTISTEP_STARTER_STOPPED_EARLY. */
    MULTISTEP_STARTER_STOPPED_EARLY("multistep integrator starter stopped early, maybe too large step size"),

    /** PROPAGATION_DIRECTION_MISMATCH. */
    PROPAGATION_DIRECTION_MISMATCH("propagation direction mismatch"),

    /** TOO_SMALL_INTEGRATION_INTERVAL. */
    TOO_SMALL_INTEGRATION_INTERVAL("too small integration interval: length = {0}"),

    /** UNKNOWN_PARAMETER. */
    UNKNOWN_PARAMETER("unknown parameter {0}"),

    /** UNMATCHED_ODE_IN_EXPANDED_SET. */
    UNMATCHED_ODE_IN_EXPANDED_SET("ode does not match the main ode set in the extended set"),

    /** NAN_APPEARING_DURING_INTEGRATION. */
    NAN_APPEARING_DURING_INTEGRATION("NaN appears during integration near time {0}"),

    /** FIND_ROOT. */
    FIND_ROOT("{0} failed to find root between {1} (g={2,number,0.0##############E0}) and {3} (g={4,number,0.0##############E0})\nLast iteration at {5} (g={6,number,0.0##############E0})");

    /** Source English format. */
    private final String sourceFormat;

    /** Simple constructor.
     * @param sourceFormat source English format to use when no
     * localized version is available
     */
    LocalizedODEFormats(final String sourceFormat) {
        this.sourceFormat = sourceFormat;
    }

    /** {@inheritDoc} */
    @Override
    public String getSourceString() {
        return sourceFormat;
    }

    /** {@inheritDoc} */
    @Override
    public String getLocalizedString(final Locale locale) {
        try {
            final String path = LocalizedODEFormats.class.getName().replaceAll("\\.", "/");
            ResourceBundle bundle =
                    ResourceBundle.getBundle("assets/" + path, locale, new UTF8Control());
            if (bundle.getLocale().getLanguage().equals(locale.getLanguage())) {
                final String translated = bundle.getString(name());
                if ((translated != null) &&
                    (translated.length() > 0) &&
                    (!translated.toLowerCase(locale).contains("missing translation"))) {
                    // the value of the resource is the translated format
                    return translated;
                }
            }

        } catch (MissingResourceException mre) { // NOPMD
            // do nothing here
        }

        // either the locale is not supported or the resource is unknown
        // don't translate and fall back to using the source format
        return sourceFormat;

    }

}
