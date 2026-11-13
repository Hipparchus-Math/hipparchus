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

package org.hipparchus.complex;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.util.CompositeFormat;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.Locale;


public class ComplexFormatTest extends ComplexFormatAbstractTest {
    @Override
    protected char getDecimalCharacter() {
        return '.';
    }

    @Override
    protected Locale getLocale() {
        return Locale.US;
    }

    @Test
    void testIssue459() {

        checkWrongParsing(",,7+,,,2i");
        checkWrongParsing(",8+,,3i");
        checkWrongParsing(",7");
        checkWrongParsing(";7");
        checkWrongParsing("#7");

        Assertions.assertTrue(CompositeFormat.getDefaultNumberFormat(Locale.US).isGroupingUsed());
        Assertions.assertFalse(ComplexFormat.getComplexFormat().getRealFormat().isGroupingUsed());
        Assertions.assertFalse(ComplexFormat.getComplexFormat().getImaginaryFormat().isGroupingUsed());

    }

    private void checkWrongParsing(final String input) {
        try {
            final ComplexFormat format = new ComplexFormat(CompositeFormat.getDefaultNumberFormat(Locale.US));
            format.parse(input);
            Assertions.fail("an exception should have been thrown");
        } catch (MathIllegalStateException mise) {
            Assertions.assertEquals(LocalizedCoreFormats.CANNOT_PARSE_AS_TYPE, mise.getSpecifier());
        }
    }

}
