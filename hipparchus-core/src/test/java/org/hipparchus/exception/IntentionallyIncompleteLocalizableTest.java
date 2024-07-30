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

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.Locale;

/**
 * Intentionally incomplete localizable, for testing purposes
 * @author Luc Maisonobe
 */
public class IntentionallyIncompleteLocalizableTest
{
    @Test
    public void testFixedLocale() {
        Assertions.assertEquals("message without argument",
                                IntentionallyIncompleteLocalizable.MESSAGE_WITHOUT_ARGUMENT.getLocalizedString(Locale.ENGLISH));
        Assertions.assertEquals("message with one argument {0}",
                                IntentionallyIncompleteLocalizable.MESSAGE_WITH_ONE_ARGUMENT.getLocalizedString(Locale.ENGLISH));
        Assertions.assertEquals("message with two arguments {0}, {1}",
                                IntentionallyIncompleteLocalizable.MESSAGE_WITH_TWO_ARGUMENTS.getLocalizedString(Locale.ENGLISH));
    }

    @Test
    public void testUnsupportedLocale() {
        Assertions.assertEquals("message without argument",
                                IntentionallyIncompleteLocalizable.MESSAGE_WITHOUT_ARGUMENT.getLocalizedString(Locale.TRADITIONAL_CHINESE));
        Assertions.assertEquals("message with one argument {0}",
                                IntentionallyIncompleteLocalizable.MESSAGE_WITH_ONE_ARGUMENT.getLocalizedString(Locale.TRADITIONAL_CHINESE));
        Assertions.assertEquals("message with two arguments {0}, {1}",
                                IntentionallyIncompleteLocalizable.MESSAGE_WITH_TWO_ARGUMENTS.getLocalizedString(Locale.TRADITIONAL_CHINESE));
    }

    @Test
    public void testIncompleteTranslations() {
        Assertions.assertEquals("message sans argument",
                                IntentionallyIncompleteLocalizable.MESSAGE_WITHOUT_ARGUMENT.getLocalizedString(Locale.FRENCH));
        Assertions.assertEquals("message with one argument {0}",
                                IntentionallyIncompleteLocalizable.MESSAGE_WITH_ONE_ARGUMENT.getLocalizedString(Locale.FRENCH));
        Assertions.assertEquals("message with two arguments {0}, {1}",
                                IntentionallyIncompleteLocalizable.MESSAGE_WITH_TWO_ARGUMENTS.getLocalizedString(Locale.FRENCH));
    }

}
