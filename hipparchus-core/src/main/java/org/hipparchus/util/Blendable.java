package org.hipparchus.util;

import org.hipparchus.exception.MathIllegalArgumentException;

/**
 * Interface representing classes that can blend with other instances of themselves using a given smoothstep function.
 *
 * @param <B> blendable class
 */
public interface Blendable<B> {

    /**
     * Blend arithmetically this instance with another one.
     *
     * @param other other instance to blend arithmetically with
     * @param blendingValue value from smoothstep function B(x). It is expected to be between [0:1] and will
     *         throw an exception otherwise.
     * @return this * (1 - B(x)) + other * B(x)
     * @throws MathIllegalArgumentException if blending value is not within [0:1]
     */
    B blendArithmeticallyWith(B other, double blendingValue) throws MathIllegalArgumentException;
}
