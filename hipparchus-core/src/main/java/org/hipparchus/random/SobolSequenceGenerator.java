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
package org.hipparchus.random;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.Charset;
import java.util.Arrays;
import java.util.NoSuchElementException;
import java.util.StringTokenizer;

import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.exception.MathRuntimeException;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;

/**
 * Implementation of a Sobol sequence.
 * <p>
 * A Sobol sequence is a low-discrepancy sequence with the property that for all values of N,
 * its subsequence (x1, ... xN) has a low discrepancy. It can be used to generate pseudo-random
 * points in a space S, which are equi-distributed.
 * <p>
 * The implementation already comes with support for up to 21201 dimensions with direction numbers
 * calculated from <a href="http://web.maths.unsw.edu.au/~fkuo/sobol/">Stephen Joe and Frances Kuo</a>.
 * <p>
 * The generator supports two modes:
 * <ul>
 *   <li>sequential generation of points: {@link #nextVector()}</li>
 *   <li>random access to the i-th point in the sequence: {@link #skipTo(int)}</li>
 * </ul>
 *
 * @see <a href="http://en.wikipedia.org/wiki/Sobol_sequence">Sobol sequence (Wikipedia)</a>
 * @see <a href="http://web.maths.unsw.edu.au/~fkuo/sobol/">Sobol sequence direction numbers</a>
 *
 */
public class SobolSequenceGenerator implements RandomVectorGenerator {

    /** The number of bits to use. */
    private static final int BITS = 52;

    /** The scaling factor. */
    private static final double SCALE = FastMath.pow(2, BITS);

    /** The maximum supported space dimension. */
    private static final int MAX_DIMENSION = 21201;

    /** The resource containing the direction numbers. */
    private static final String RESOURCE_NAME = "/assets/org/hipparchus/random/new-joe-kuo-6.21201";

    /** Character set for file input. */
    private static final String FILE_CHARSET = "US-ASCII";

    /** Space dimension. */
    private final int dimension;

    /** The current index in the sequence. */
    private int count;

    /** The direction vector for each component. */
    private final long[][] direction;

    /** The current state. */
    private final long[] x;

    /**
     * Construct a new Sobol sequence generator for the given space dimension.
     *
     * @param dimension the space dimension
     * @throws MathIllegalArgumentException if the space dimension is outside the allowed range of [1, 1000]
     */
    public SobolSequenceGenerator(final int dimension) throws MathIllegalArgumentException {
        MathUtils.checkRangeInclusive(dimension, 1, MAX_DIMENSION);

        // initialize the other dimensions with direction numbers from a resource
        try (InputStream is = getClass().getResourceAsStream(RESOURCE_NAME)) {
            if (is == null) {
                throw MathRuntimeException.createInternalError();
            }

            this.dimension = dimension;

            // init data structures
            direction = new long[dimension][BITS + 1];
            x = new long[dimension];

            initFromStream(is);
        } catch (IOException | MathIllegalStateException e) {
            // the internal resource file could not be parsed -> should not happen
            throw MathRuntimeException.createInternalError(e);
        }
    }

    /**
     * Construct a new Sobol sequence generator for the given space dimension with
     * direction vectors loaded from the given stream.
     * <p>
     * The expected format is identical to the files available from
     * <a href="http://web.maths.unsw.edu.au/~fkuo/sobol/">Stephen Joe and Frances Kuo</a>.
     * The first line will be ignored as it is assumed to contain only the column headers.
     * The columns are:
     * <ul>
     *  <li>d: the dimension</li>
     *  <li>s: the degree of the primitive polynomial</li>
     *  <li>a: the number representing the coefficients</li>
     *  <li>m: the list of initial direction numbers</li>
     * </ul>
     * Example:
     * <pre>
     * d       s       a       m_i
     * 2       1       0       1
     * 3       2       1       1 3
     * </pre>
     * <p>
     * The input stream <i>must</i> be an ASCII text containing one valid direction vector per line.
     *
     * @param dimension the space dimension
     * @param is the stream to read the direction vectors from
     * @throws MathIllegalArgumentException if the space dimension is &lt; 1
     * @throws MathIllegalArgumentException if the space dimension is outside the range [1, max], where
     *   max refers to the maximum dimension found in the input stream
     * @throws MathIllegalStateException if the content in the stream could not be parsed successfully
     * @throws IOException if an error occurs while reading from the input stream
     */
    public SobolSequenceGenerator(final int dimension, final InputStream is)
            throws MathIllegalArgumentException, MathIllegalStateException, IOException {

        if (dimension < 1) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.NUMBER_TOO_SMALL,
                                                   dimension, 1);
        }

        this.dimension = dimension;

        // init data structures
        direction = new long[dimension][BITS + 1];
        x = new long[dimension];

        // initialize the other dimensions with direction numbers from the stream
        int lastDimension = initFromStream(is);
        MathUtils.checkRangeInclusive(dimension, 1, lastDimension);
    }

    /**
     * Load the direction vector for each dimension from the given stream.
     * <p>
     * The input stream <i>must</i> be an ASCII text containing one
     * valid direction vector per line.
     *
     * @param is the input stream to read the direction vector from
     * @return the last dimension that has been read from the input stream
     * @throws IOException if the stream could not be read
     * @throws MathIllegalStateException if the content could not be parsed successfully
     */
    private int initFromStream(final InputStream is) throws MathIllegalStateException, IOException {

        // special case: dimension 1 -> use unit initialization
        for (int i = 1; i <= BITS; i++) {
            direction[0][i] = 1L << (BITS - i);
        }

        final Charset charset = Charset.forName(FILE_CHARSET);
        int dim = -1;

        try (BufferedReader reader = new BufferedReader(new InputStreamReader(is, charset))) {
            // ignore first line
            reader.readLine();

            int lineNumber = 2;
            int index = 1;
            for (String line = reader.readLine(); line != null; line = reader.readLine()) {
                StringTokenizer st = new StringTokenizer(line, " ");
                try {
                    dim = Integer.parseInt(st.nextToken());
                    if (dim >= 2 && dim <= dimension) { // we have found the right dimension
                        final int s = Integer.parseInt(st.nextToken());
                        final int a = Integer.parseInt(st.nextToken());
                        final int[] m = new int[s + 1];
                        for (int i = 1; i <= s; i++) {
                            m[i] = Integer.parseInt(st.nextToken());
                        }
                        initDirectionVector(index++, a, m);
                    }

                    if (dim > dimension) {
                        return dim;
                    }
                } catch (NoSuchElementException|NumberFormatException e) {
                    throw new MathIllegalStateException(e, LocalizedCoreFormats.CANNOT_PARSE,
                                                        line, lineNumber);
                }
                lineNumber++;
            }
        }

        return dim;
    }

    /**
     * Calculate the direction numbers from the given polynomial.
     *
     * @param d the dimension, zero-based
     * @param a the coefficients of the primitive polynomial
     * @param m the initial direction numbers
     */
    private void initDirectionVector(final int d, final int a, final int[] m) {
        final int s = m.length - 1;
        for (int i = 1; i <= s; i++) {
            direction[d][i] = ((long) m[i]) << (BITS - i);
        }
        for (int i = s + 1; i <= BITS; i++) {
            direction[d][i] = direction[d][i - s] ^ (direction[d][i - s] >> s);
            for (int k = 1; k <= s - 1; k++) {
                direction[d][i] ^= ((a >> (s - 1 - k)) & 1) * direction[d][i - k];
            }
        }
    }

    /** {@inheritDoc} */
    @Override
    public double[] nextVector() {
        final double[] v = new double[dimension];
        if (count == 0) {
            count++;
            return v;
        }

        // find the index c of the rightmost 0
        int c = 1;
        int value = count - 1;
        while ((value & 1) == 1) {
            value >>= 1;
            c++;
        }

        for (int i = 0; i < dimension; i++) {
            x[i] ^= direction[i][c];
            v[i] = x[i] / SCALE;
        }
        count++;
        return v;
    }

    /**
     * Skip to the i-th point in the Sobol sequence.
     * <p>
     * This operation can be performed in O(1).
     *
     * @param index the index in the sequence to skip to
     * @return the i-th point in the Sobol sequence
     * @throws MathIllegalArgumentException if index &lt; 0
     */
    public double[] skipTo(final int index) throws MathIllegalArgumentException {
        if (index == 0) {
            // reset x vector
            Arrays.fill(x, 0);
        } else {
            final int i = index - 1;
            final long grayCode = i ^ (i >> 1); // compute the gray code of i = i XOR floor(i / 2)
            for (int j = 0; j < dimension; j++) {
                long result = 0;
                for (int k = 1; k <= BITS; k++) {
                    final long shift = grayCode >> (k - 1);
                    if (shift == 0) {
                        // stop, as all remaining bits will be zero
                        break;
                    }
                    // the k-th bit of i
                    final long ik = shift & 1;
                    result ^= ik * direction[j][k];
                }
                x[j] = result;
            }
        }
        count = index;
        return nextVector();
    }

    /**
     * Returns the index i of the next point in the Sobol sequence that will be returned
     * by calling {@link #nextVector()}.
     *
     * @return the index of the next point
     */
    public int getNextIndex() {
        return count;
    }

}
