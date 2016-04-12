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
package org.hipparchus.stat.descriptive;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.util.MathUtils;

/**
 * Extends the definition of {@link UnivariateStatistic} for statistics
 * based on quantiles.
 */
public interface QuantiledUnivariateStatistic extends UnivariateStatistic {

    /**
     * Returns the value of the quantile field (determines what percentile is
     * computed when evaluate() is called with no quantile argument).
     *
     * @return quantile set while construction or {@link #setQuantile(double)}
     */
    double getQuantile();

    /**
     * Sets the value of the quantile field (determines what percentile is
     * computed when evaluate() is called with no quantile argument).
     *
     * @param p a value between 0 < p <= 100
     * @throws MathIllegalArgumentException  if p is not greater than 0 and less
     * than or equal to 100
     */
    void setQuantile(final double p) throws MathIllegalArgumentException;

    /**
     * Returns an estimate of the <code>p</code>th percentile of the values
     * in the <code>values</code> array.
     * <p>
     * <ul>
     * <li>Returns <code>Double.NaN</code> if <code>values</code> has length
     * <code>0</code></li>
     * <li>Returns (for any value of <code>p</code>) <code>values[0]</code>
     *  if <code>values</code> has length <code>1</code></li>
     * <li>Throws <code>MathIllegalArgumentException</code> if <code>values</code>
     * is null or p is not a valid quantile value (p must be greater than 0
     * and less than or equal to 100) </li>
     * </ul>
     * <p>
     * The default implementation delegates to
     * <code>evaluate(double[], int, int, double)</code> in the natural way.
     *
     * @param values input array of values
     * @param p the percentile value to compute
     * @return the percentile value or Double.NaN if the array is empty
     * @throws MathIllegalArgumentException if <code>values</code> is null or p is invalid
     */
    default double evaluate(final double[] values, final double p)
        throws MathIllegalArgumentException {
        MathUtils.checkNotNull(values, LocalizedCoreFormats.INPUT_ARRAY);
        return evaluate(values, 0, values.length, p);
    }

    /**
     * Returns an estimate of the <code>p</code>th percentile of the values
     * in the <code>values</code> array, starting with the element in (0-based)
     * position <code>begin</code> in the array and including <code>length</code>
     * values.
     * <p>
     * Calls to this method do not modify the internal <code>quantile</code>
     * state of this statistic.
     * <p>
     * <ul>
     * <li>Returns <code>Double.NaN</code> if <code>length = 0</code></li>
     * <li>Returns (for any value of <code>p</code>) <code>values[begin]</code>
     *  if <code>length = 1 </code></li>
     * <li>Throws <code>MathIllegalArgumentException</code> if <code>values</code>
     *  is null , <code>begin</code> or <code>length</code> is invalid, or
     * <code>p</code> is not a valid quantile value (p must be greater than 0
     * and less than or equal to 100)</li>
     * </ul>
     *
     * @param values array of input values
     * @param p  the percentile to compute
     * @param begin  the first (0-based) element to include in the computation
     * @param length  the number of array elements to include
     * @return  the percentile value
     * @throws MathIllegalArgumentException if the parameters are not valid or the
     * input array is null
     */
    double evaluate(final double[] values, final int begin,
                    final int length, final double p);

    /**
     * Returns a copy of the statistic with the same internal state.
     *
     * @return a copy of the statistic
     */
    @Override
    QuantiledUnivariateStatistic copy();

}
