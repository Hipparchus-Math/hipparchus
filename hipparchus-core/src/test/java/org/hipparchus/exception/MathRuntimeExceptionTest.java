package org.hipparchus.exception;

import org.hamcrest.CoreMatchers;
import org.hamcrest.MatcherAssert;
import org.junit.Test;

import java.io.PrintWriter;
import java.io.StringWriter;

/**
 * Unit tests for {@link MathRuntimeException}.
 *
 * @author Evan Ward
 */
public class MathRuntimeExceptionTest {

    /**
     * Check that a helpful message and stack trace is still generated even when there is
     * an error while building the exception message.
     */
    @Test
    public void testGetMessageError() {
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
                    LocalizedCoreFormats.NOT_BRACKETING_INTERVAL,
                    bad);
        } catch (MathRuntimeException e) {
            // verify
            e.printStackTrace(writer);
            String message = buffer.toString();
            System.out.println(message);
            // check original reason is preserved
            MatcherAssert.assertThat(message,
                    CoreMatchers.containsString("interval does not bracket a root"));
            MatcherAssert.assertThat(message,
                    CoreMatchers.containsString("MathRuntimeException"));
            // check exception during formatting is preserved
            MatcherAssert.assertThat(message,
                    CoreMatchers.containsString("toString failed"));
        }
    }

}
