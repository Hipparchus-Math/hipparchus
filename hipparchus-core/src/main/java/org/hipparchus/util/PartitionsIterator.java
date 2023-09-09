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
package org.hipparchus.util;

import java.lang.reflect.Array;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Queue;

/** Iterator for generating partitions.
 * <p>
 * This class implements the iterative algorithm described in
 * <a href="https://academic.oup.com/comjnl/article/32/3/281/331557">Short Note:
 * A Fast Iterative Algorithm for Generating Set Partitions</a>
 * by B. Djokić, M. Miyakawa, S. Sekiguchi, I. Semba, and I. Stojmenović
 * (The Computer Journal, Volume 32, Issue 3, 1989, Pages 281–282,
 * <a href="https://doi.org/10.1093/comjnl/32.3.281">https://doi.org/10.1093/comjnl/32.3.281</a>
 * </p>
 * @param <T> type of the elements
 * @since 2.2
 */
class PartitionsIterator<T> implements Iterator<List<T>[]> {

    /** List to partition. */
    private final List<T> list;

    /** Number of elements to partition. */
    private final int   n;

    /** Mapping from elements indices to parts indices. */
    private final int[] partIndex;

    /** Backtracking array. */
    private final int[] backTrack;

    /** Current part index. */
    private int   r;

    /** Current backtrack index. */
    private int   j;

    /** Pending parts already generated. */
    private final Queue<List<T>[]> pending;

    /** Indicator for exhausted partitions. */
    private boolean exhausted;

    /** Simple constructor.
     * @param list list to partition
     */
    PartitionsIterator(final List<T> list) {

        this.list      = list;
        this.n         = list.size();
        this.partIndex = new int[list.size()];
        this.backTrack = new int[list.size() - 1];
        this.r         = 0;
        this.j         = 0;
        this.pending   = new ArrayDeque<>(n);

        // generate a first set of partitions
        generate();

    }

    /** Generate one set of partitions.
     */
    private void generate() {

        // put elements in the first part
        while (r < n - 2) {
            partIndex[++r] = 0;
            backTrack[++j] = r;
        }

        // generate partitions
        for (int i = 0; i < n - j; ++i) {

            // fill-up final element
            partIndex[n - 1] = i;

            // count the number of parts in this partition
            int max = 0;
            for (final int index : partIndex) {
                max = FastMath.max(max, index);
            }

            // prepare storage
            @SuppressWarnings("unchecked")
            final List<T>[] partition = (List<T>[]) Array.newInstance(List.class, max + 1);
            for (int k = 0; k < partition.length; ++k) {
                partition[k] = new ArrayList<>(n);
            }

            // distribute elements in the parts
            for (int k = 0; k < partIndex.length; ++k) {
                partition[partIndex[k]].add(list.get(k));
            }

            // add the generated partition to the pending queue
            pending.add(partition);

        }

        // backtrack to generate next partition
        r = backTrack[j];
        partIndex[r]++;
        if (partIndex[r] > r - j) {
            --j;
        }

        // keep track of end of generation
        exhausted = r == 0;

    }

    /** {@inheritDoc} */
    @Override
    public boolean hasNext() {
        return !(exhausted && pending.isEmpty());
    }

    /** {@inheritDoc} */
    @Override
    public List<T>[] next() {

        if (pending.isEmpty()) {
            // we need to generate more partitions
            if (exhausted) {
                throw new NoSuchElementException();
            }
            generate();
        }

        return pending.remove();

    }

}
