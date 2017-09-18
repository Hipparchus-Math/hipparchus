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

/*
 * This is not the original file distributed by the Apache Software Foundation
 * It has been modified by the Hipparchus project
 */
package org.hipparchus.exception;


import java.lang.reflect.InvocationTargetException;
import java.text.MessageFormat;
import java.util.Enumeration;
import java.util.Locale;
import java.util.ResourceBundle;

import org.junit.Assert;
import org.junit.Test;

public abstract class LocalizedFormatsAbstractTest {

    protected abstract Class<? extends Enum<?>> getFormatsClass();
    protected abstract int getExpectedNumber();

    private Localizable[] getValues() {
        Localizable[] localizable = null;
        try {
            Object[] a = (Object[]) getFormatsClass().getMethod("values").invoke(null);
            localizable = new Localizable[a.length];
            for (int i = 0; i < a.length; ++i) {
                localizable[i] = (Localizable) a[i];
            }
        } catch (NoSuchMethodException | SecurityException | IllegalAccessException |
                        IllegalArgumentException | InvocationTargetException e) {
            Assert.fail(e.getLocalizedMessage());
        }
        return localizable;
    }

    private Localizable valueOf(String s) {
        Localizable localizable = null;
        try {
            localizable = (Localizable) getFormatsClass().getMethod("valueOf", String.class).invoke(null, s);
        } catch (IllegalArgumentException iae) {
            localizable = null;
        } catch (NoSuchMethodException | SecurityException | IllegalAccessException | InvocationTargetException e) {
            Assert.fail(s + " <-> " + e.getLocalizedMessage());
        }
        return localizable;
    }

    @Test
    public void testMessageNumber() {
        Assert.assertEquals(getExpectedNumber(), getValues().length);
    }

    @Test
    public void testAllKeysPresentInPropertiesFiles() {
        Class<? extends Enum<?>> c = getFormatsClass();
        final String path = c.getName().replaceAll("\\.", "/");
        for (final String language : new String[] { "fr" } ) {
            ResourceBundle bundle =
                ResourceBundle.getBundle("assets/" + path, new Locale(language), c.getClassLoader());
            for (Localizable message : getValues()) {
                final String messageKey = message.toString();
                boolean keyPresent = false;
                for (final Enumeration<String> keys = bundle.getKeys(); keys.hasMoreElements();) {
                    keyPresent |= messageKey.equals(keys.nextElement());
                }
                Assert.assertTrue("missing key \"" + message.toString() + "\" for language " + language,
                                  keyPresent);
            }
            Assert.assertEquals(language, bundle.getLocale().getLanguage());
        }

    }

    @Test
    public void testAllPropertiesCorrespondToKeys() {
        Class<? extends Enum<?>> c = getFormatsClass();
        final String path = c.getName().replaceAll("\\.", "/");
        for (final String language : new String[] { "fr" } ) {
            ResourceBundle bundle =
                ResourceBundle.getBundle("assets/" + path, new Locale(language), c.getClassLoader());
            for (final Enumeration<String> keys = bundle.getKeys(); keys.hasMoreElements();) {
                final String propertyKey = keys.nextElement();
                try {
                    Assert.assertNotNull(valueOf(propertyKey));
                } catch (IllegalArgumentException iae) {
                    Assert.fail("unknown key \"" + propertyKey + "\" in language " + language);
                }
            }
            Assert.assertEquals(language, bundle.getLocale().getLanguage());
        }

    }

    @Test
    public void testNoMissingFrenchTranslation() {
        for (Localizable message : getValues()) {
            String translated = message.getLocalizedString(Locale.FRENCH);
            Assert.assertFalse(message.toString(), translated.toLowerCase().contains("missing translation"));
        }
    }

    @Test
    public void testNoOpEnglishTranslation() {
        for (Localizable message : getValues()) {
            String translated = message.getLocalizedString(Locale.ENGLISH);
            Assert.assertEquals(message.getSourceString(), translated);
        }
    }

    @Test
    public void testVariablePartsConsistency() {
        for (final String language : new String[] { "fr" } ) {
            Locale locale = new Locale(language);
            for (Localizable message : getValues()) {
                MessageFormat source     = new MessageFormat(message.getSourceString());
                MessageFormat translated = new MessageFormat(message.getLocalizedString(locale));
                Assert.assertEquals(message.toString() + " (" + language + ")",
                                    source.getFormatsByArgumentIndex().length,
                                    translated.getFormatsByArgumentIndex().length);
            }
        }
    }

}
