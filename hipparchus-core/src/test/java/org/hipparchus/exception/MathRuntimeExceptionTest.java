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

import org.hamcrest.CoreMatchers;
import org.hipparchus.util.FastMath;
import org.junit.jupiter.api.Test;

import java.io.PrintWriter;
import java.io.StringWriter;

import static org.hamcrest.MatcherAssert.assertThat;

/**
 * Unit tests for {@link MathRuntimeException}.
 *
 * @author Evan Ward
 */
class MathRuntimeExceptionTest {

    /**
     * Check that a helpful message and stack trace is still generated even when there is
     * an error while building the exception message.
     */
    @Test
    void testGetMessageError() {
        // setup
        Object bad = new Object() {
            @Override
            public String toString() {
                throw new RuntimeException("toString failed");
            }
        };
        StringWriter buffer = new StringWriter();
        PrintWriter writer = new PrintWriter(buffer);

        // action
        try {
            throw new MathRuntimeException(
                    LocalizedCoreFormats.URL_CONTAINS_NO_DATA,
                    bad);
        } catch (MathRuntimeException e) {
            // verify
            e.printStackTrace(writer);
            String message = buffer.toString();
            // check original reason is preserved
            assertThat(message,
                    CoreMatchers.containsString("contains no data"));
            assertThat(message,
                    CoreMatchers.containsString("MathRuntimeException"));
            // check exception during formatting is preserved
            assertThat(message,
                    CoreMatchers.containsString("toString failed"));
        }
    }

    /** Check the bracketing exception message uses full precision. */
    @Test
    void testGetMessageDecimalFormat() {
        // setup
        double a = FastMath.nextUp(1.0), b = FastMath.nextDown(1.0);
        double fa = -Double.MIN_NORMAL, fb = -12.345678901234567e-10;

        // action
        String message = new MathRuntimeException(
                LocalizedCoreFormats.NOT_BRACKETING_INTERVAL, a, b, fa, fb)
                .getMessage();

        // verify
        String expected = "interval does not bracket a root: " +
                "f(1.0000000000000002E0) = -22.250738585072014E-309, " +
                "f(999.9999999999999E-3) = -1.2345678901234566E-9";
        assertThat( message, CoreMatchers.is( expected) );
    }

}
