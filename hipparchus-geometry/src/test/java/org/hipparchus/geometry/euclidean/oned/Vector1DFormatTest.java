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

package org.hipparchus.geometry.euclidean.oned;

import org.hipparchus.exception.MathIllegalStateException;
import org.junit.jupiter.api.Test;

import java.text.NumberFormat;
import java.text.ParsePosition;
import java.util.Locale;

import static org.junit.jupiter.api.Assertions.*;


public class Vector1DFormatTest extends Vector1DFormatAbstractTest {

    @Override
    protected char getDecimalCharacter() {
        return '.';
    }

    @Override
    protected Locale getLocale() {
        return Locale.US;
    }

    @Test
    void testGetAvailableLocales() {
        Locale[] locales = Vector1DFormat.getAvailableLocales();
        assertNotNull(locales);
        assertEquals(NumberFormat.getAvailableLocales().length, locales.length);
    }

    @Test
    void testParse() {
        Vector1DFormat format = Vector1DFormat.getVector1DFormat();
        Vector1D expected = new Vector1D(1.0);
        assertEquals(expected, format.parse("{1}"));
        assertThrows(MathIllegalStateException.class, () -> format.parse("invalid"));
    }

    @Test
    void testParseCoordinatesFailure() {
        // This targets the parseCoordinates method in VectorFormat via a concrete subclass
        Vector1DFormat format = Vector1DFormat.getVector1DFormat();
        
        // Missing prefix
        ParsePosition pos = new ParsePosition(0);
        assertNull(format.parse("1}", pos));
        
        // Invalid number
        pos = new ParsePosition(0);
        assertNull(format.parse("{abc}", pos));
        assertEquals(0, pos.getIndex()); // Should reset to initial index
        
        // Missing suffix
        pos = new ParsePosition(0);
        assertNull(format.parse("{1", pos));
    }
}
