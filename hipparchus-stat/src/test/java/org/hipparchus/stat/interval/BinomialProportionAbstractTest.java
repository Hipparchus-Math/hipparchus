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
package org.hipparchus.stat.interval;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertThrows;

/**
 * Test cases for the BinomialProportion implementations.
 */
public abstract class BinomialProportionAbstractTest {

    @FunctionalInterface
    public interface BinomialProportionMethod {
        ConfidenceInterval calculate(int trials, double probability, double confidenceLevel);
    }

    protected BinomialProportionMethod testMethod;

    private final int trials                  = 500;
    private final double probabilityOfSuccess = 0.1;
    private final double confidenceLevel      = 0.9;

    protected abstract BinomialProportionMethod getBinomialProportionMethod();

    /**
     * Returns the confidence interval for the given statistic with the following values:
     *
     * <ul>
     *  <li>trials: 500</li>
     *  <li>probabilityOfSuccess: 0.1</li>
     *  <li>confidenceLevel: 0.9</li>
     * </ul>
     * @return the Confidence Interval for the given values
     */
    protected ConfidenceInterval createStandardTestInterval() {
        return testMethod.calculate(trials, probabilityOfSuccess, confidenceLevel);
    }

    @BeforeEach
    public void setUp() {
        testMethod = getBinomialProportionMethod();
    }

    @Test
    public void testZeroConfidencelevel() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            testMethod.calculate(trials, probabilityOfSuccess, 0d);
        });
    }

    @Test
    public void testOneConfidencelevel() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            testMethod.calculate(trials, probabilityOfSuccess, 1d);
        });
    }

    @Test
    public void testZeroTrials() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            testMethod.calculate(0, 0, confidenceLevel);
        });
    }

    @Test
    public void testNegativeSuccesses() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            testMethod.calculate(trials, -1, confidenceLevel);
        });
    }

    @Test
    public void testSuccessesExceedingTrials() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            testMethod.calculate(trials, trials + 1, confidenceLevel);
        });
    }
}
