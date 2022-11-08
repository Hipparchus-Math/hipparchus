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
package org.hipparchus.stat.correlation;

import java.util.Arrays;

import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.linear.BlockRealMatrix;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathArrays;

/**
 * Implementation of Kendall's Tau-b rank correlation.
 * <p>
 * A pair of observations (<var>x<sub>1</sub></var>, <var>y<sub>1</sub></var>) and
 * (<var>x<sub>2</sub></var>, <var>y<sub>2</sub></var>) are considered <i>concordant</i> if
 * <var>x<sub>1</sub></var> &lt; <var>x<sub>2</sub></var> and <var>y<sub>1</sub></var> &lt; <var>y<sub>2</sub></var>
 * or <var>x<sub>2</sub></var> &lt; <var>x<sub>1</sub></var> and <var>y<sub>2</sub></var> &lt; <var>y<sub>1</sub></var>.
 * The pair is <i>discordant</i> if <var><var>x<sub>1</sub></var></var> &lt; <var><var>x<sub>2</sub></var></var> and
 * <var><var>y<sub>2</sub></var></var> &lt; <var><var>y<sub>1</sub></var></var> or <var><var>x<sub>2</sub></var></var> &lt; <var><var>x<sub>1</sub></var></var> and
 * <var><var>y<sub>1</sub></var></var> &lt; <var><var>y<sub>2</sub></var></var>.  If either <var><var>x<sub>1</sub></var></var> = <var><var>x<sub>2</sub></var></var>
 * or <var><var>y<sub>1</sub></var></var> = <var><var>y<sub>2</sub></var></var>, the pair is neither concordant nor
 * discordant.
 * <p>
 * Kendall's Tau-b is defined as:
 * <pre>
 * <var><var>tau<sub>b</sub></var></var> = (<var><var>n<sub>c</sub></var></var> - <var><var>n<sub>d</sub></var></var>) / sqrt((<var>n<sub>0</sub></var> - <var>n<sub>1</sub></var>) * (<var>n<sub>0</sub></var> - <var>n<sub>2</sub></var>))
 * </pre>
 * <p>
 * where:
 * <ul>
 *     <li><var>n<sub>0</sub></var> = n * (n - 1) / 2</li>
 *     <li><var>n<sub>c</sub></var> = Number of concordant pairs</li>
 *     <li><var>n<sub>d</sub></var> = Number of discordant pairs</li>
 *     <li><var>n<sub>1</sub></var> = sum of <var>t<sub>i</sub></var> * (<var>t<sub>i</sub></var> - 1) / 2 for all i</li>
 *     <li><var>n<sub>2</sub></var> = sum of <var>u<sub>j</sub></var> * (<var>u<sub>j</sub></var> - 1) / 2 for all j</li>
 *     <li><var>t<sub>i</sub></var> = Number of tied values in the i<sup>th</sup> group of ties in x</li>
 *     <li><var>u<sub>j</sub></var> = Number of tied values in the j<sup>th</sup> group of ties in y</li>
 * </ul>
 * <p>
 * This implementation uses the O(n log n) algorithm described in
 * William R. Knight's 1966 paper "A Computer Method for Calculating
 * Kendall's Tau with Ungrouped Data" in the Journal of the American
 * Statistical Association.
 *
 * @see <a href="http://en.wikipedia.org/wiki/Kendall_tau_rank_correlation_coefficient">
 * Kendall tau rank correlation coefficient (Wikipedia)</a>
 * @see <a href="http://www.jstor.org/stable/2282833">A Computer
 * Method for Calculating Kendall's Tau with Ungrouped Data</a>
 */
public class KendallsCorrelation {

    /** correlation matrix */
    private final RealMatrix correlationMatrix;

    /**
     * Create a KendallsCorrelation instance without data.
     */
    public KendallsCorrelation() {
        correlationMatrix = null;
    }

    /**
     * Create a KendallsCorrelation from a rectangular array
     * whose columns represent values of variables to be correlated.
     *
     * @param data rectangular array with columns representing variables
     * @throws IllegalArgumentException if the input data array is not
     * rectangular with at least two rows and two columns.
     */
    public KendallsCorrelation(double[][] data) {
        this(MatrixUtils.createRealMatrix(data));
    }

    /**
     * Create a KendallsCorrelation from a RealMatrix whose columns
     * represent variables to be correlated.
     *
     * @param matrix matrix with columns representing variables to correlate
     */
    public KendallsCorrelation(RealMatrix matrix) {
        correlationMatrix = computeCorrelationMatrix(matrix);
    }

    /**
     * Returns the correlation matrix.
     *
     * @return correlation matrix
     */
    public RealMatrix getCorrelationMatrix() {
        return correlationMatrix;
    }

    /**
     * Computes the Kendall's Tau rank correlation matrix for the columns of
     * the input matrix.
     *
     * @param matrix matrix with columns representing variables to correlate
     * @return correlation matrix
     */
    public RealMatrix computeCorrelationMatrix(final RealMatrix matrix) {
        int nVars = matrix.getColumnDimension();
        RealMatrix outMatrix = new BlockRealMatrix(nVars, nVars);
        for (int i = 0; i < nVars; i++) {
            for (int j = 0; j < i; j++) {
                double corr = correlation(matrix.getColumn(i), matrix.getColumn(j));
                outMatrix.setEntry(i, j, corr);
                outMatrix.setEntry(j, i, corr);
            }
            outMatrix.setEntry(i, i, 1d);
        }
        return outMatrix;
    }

    /**
     * Computes the Kendall's Tau rank correlation matrix for the columns of
     * the input rectangular array.  The columns of the array represent values
     * of variables to be correlated.
     *
     * @param matrix matrix with columns representing variables to correlate
     * @return correlation matrix
     */
    public RealMatrix computeCorrelationMatrix(final double[][] matrix) {
       return computeCorrelationMatrix(new BlockRealMatrix(matrix));
    }

    /**
     * Computes the Kendall's Tau rank correlation coefficient between the two arrays.
     *
     * @param xArray first data array
     * @param yArray second data array
     * @return Returns Kendall's Tau rank correlation coefficient for the two arrays
     * @throws MathIllegalArgumentException if the arrays lengths do not match
     */
    public double correlation(final double[] xArray, final double[] yArray)
            throws MathIllegalArgumentException {

        MathArrays.checkEqualLength(xArray, yArray);

        final int n = xArray.length;
        final long numPairs = sum(n - 1);

        DoublePair[] pairs = new DoublePair[n];
        for (int i = 0; i < n; i++) {
            pairs[i] = new DoublePair(xArray[i], yArray[i]);
        }

        Arrays.sort(pairs, (p1, p2) -> {
            int compareKey = Double.compare(p1.getFirst(), p2.getFirst());
            return compareKey != 0 ? compareKey : Double.compare(p1.getSecond(), p2.getSecond());
        });

        long tiedXPairs = 0;
        long tiedXYPairs = 0;
        long consecutiveXTies = 1;
        long consecutiveXYTies = 1;
        DoublePair prev = pairs[0];
        for (int i = 1; i < n; i++) {
            final DoublePair curr = pairs[i];
            if (Double.compare(curr.getFirst(), prev.getFirst()) == 0) {
                consecutiveXTies++;
                if (Double.compare(curr.getSecond(), prev.getSecond()) == 0) {
                    consecutiveXYTies++;
                } else {
                    tiedXYPairs += sum(consecutiveXYTies - 1);
                    consecutiveXYTies = 1;
                }
            } else {
                tiedXPairs += sum(consecutiveXTies - 1);
                consecutiveXTies = 1;
                tiedXYPairs += sum(consecutiveXYTies - 1);
                consecutiveXYTies = 1;
            }
            prev = curr;
        }
        tiedXPairs += sum(consecutiveXTies - 1);
        tiedXYPairs += sum(consecutiveXYTies - 1);

        long swaps = 0;
        DoublePair[] pairsDestination = new DoublePair[n];
        for (int segmentSize = 1; segmentSize < n; segmentSize <<= 1) {
            for (int offset = 0; offset < n; offset += 2 * segmentSize) {
                int i = offset;
                final int iEnd = FastMath.min(i + segmentSize, n);
                int j = iEnd;
                final int jEnd = FastMath.min(j + segmentSize, n);

                int copyLocation = offset;
                while (i < iEnd || j < jEnd) {
                    if (i < iEnd) {
                        if (j < jEnd) {
                            if (Double.compare(pairs[i].getSecond(), pairs[j].getSecond()) <= 0) {
                                pairsDestination[copyLocation] = pairs[i];
                                i++;
                            } else {
                                pairsDestination[copyLocation] = pairs[j];
                                j++;
                                swaps += iEnd - i;
                            }
                        } else {
                            pairsDestination[copyLocation] = pairs[i];
                            i++;
                        }
                    } else {
                        pairsDestination[copyLocation] = pairs[j];
                        j++;
                    }
                    copyLocation++;
                }
            }
            final DoublePair[] pairsTemp = pairs;
            pairs = pairsDestination;
            pairsDestination = pairsTemp;
        }

        long tiedYPairs = 0;
        long consecutiveYTies = 1;
        prev = pairs[0];
        for (int i = 1; i < n; i++) {
            final DoublePair curr = pairs[i];
            if (Double.compare(curr.getSecond(), prev.getSecond()) == 0) {
                consecutiveYTies++;
            } else {
                tiedYPairs += sum(consecutiveYTies - 1);
                consecutiveYTies = 1;
            }
            prev = curr;
        }
        tiedYPairs += sum(consecutiveYTies - 1);

        final long concordantMinusDiscordant = numPairs - tiedXPairs - tiedYPairs + tiedXYPairs - 2 * swaps;
        final double nonTiedPairsMultiplied = (numPairs - tiedXPairs) * (double) (numPairs - tiedYPairs);
        return concordantMinusDiscordant / FastMath.sqrt(nonTiedPairsMultiplied);
    }

    /**
     * Returns the sum of the number from 1 .. n according to Gauss' summation formula:
     * \[ \sum\limits_{k=1}^n k = \frac{n(n + 1)}{2} \]
     *
     * @param n the summation end
     * @return the sum of the number from 1 to n
     */
    private static long sum(long n) {
        return n * (n + 1) / 2l;
    }

    /**
     * Helper data structure holding a (double, double) pair.
     */
    private static class DoublePair {
        /** The first value */
        private final double first;
        /** The second value */
        private final double second;

        /**
         * @param first first value.
         * @param second second value.
         */
        DoublePair(double first, double second) {
            this.first = first;
            this.second = second;
        }

        /** @return the first value. */
        public double getFirst() {
            return first;
        }

        /** @return the second value. */
        public double getSecond() {
            return second;
        }

    }

}
