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
package org.hipparchus.exception;

import java.io.Serializable;
import java.util.Locale;
import java.util.MissingResourceException;
import java.util.ResourceBundle;

/**
 * Interface for localizable strings.
 *
 */
public interface Localizable extends Serializable {
    /**
     * Gets the source (non-localized) string.
     *
     * @return the source string.
     */
    String getSourceString();

    /**
     * Gets the localized string.
     *
     * @param locale locale into which to get the string.
     * @return the localized string or the source string if no
     * localized version is available.
     */
    String getLocalizedString(Locale locale);

    /**
     * Gets the localized string.
     *
     * @param baseName base name of the resource bundle
     * @param key key of the item in the bundle
     * @param locale locale into which to get the string.
     * @return the localized string or the source string if no
     * localized version is available.
     */
    default String getLocalizedString(final String baseName, final String key, final Locale locale) {

        try {
            ResourceBundle bundle;
            try {
                bundle = ResourceBundle.getBundle(baseName, locale, new UTF8Control());
            } catch (UnsupportedOperationException uoe) {
                // fix for Java 9+ on module path
                // (see issue https://github.com/Hipparchus-Math/hipparchus/issues/392)
                bundle = ResourceBundle.getBundle(baseName, locale);
            }

            if (bundle.getLocale().getLanguage().equals(locale.getLanguage())) {
                final String translated = bundle.getString(key);
                if (!(translated.isEmpty() || translated.toLowerCase(locale).contains("missing translation"))) {
                    // the value of the resource is the translated format
                    return translated;
                }
            }
        } catch (MissingResourceException mre) { // NOPMD
            // do nothing here
        }

        // either the locale is not supported or the resource is not translated, or
        // it is unknown: don't translate and fall back to using the source format
        return getSourceString();

    }

}
