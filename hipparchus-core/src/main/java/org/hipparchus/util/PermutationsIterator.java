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

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

/** Iterator for generating permutations.
 * <p>
 * This class implements the Steinhaus–Johnson–Trotter algorithm
 * with Even's speedup
 * <a href="https://en.wikipedia.org/wiki/Steinhaus%E2%80%93Johnson%E2%80%93Trotter_algorithm">Steinhaus–Johnson–Trotter algorithm</a>
 * </p>
 * @param <T> type of the elements
 * @since 2.2
 */
class PermutationsIterator<T> implements Iterator<List<T>> {

    /** Current permuted list. */
    private final List<T> permuted;

    /** Value markers. */
    private final int[] value;

    /** Directions markers. */
    private final int[] direction;

    /** Indicator for exhausted partitions. */
    private boolean exhausted;

    /** Simple constructor.
     * @param list list to permute (will not be touched)
     */
    PermutationsIterator(final List<T> list) {

        this.permuted  = new ArrayList<>(list);

        this.value     = new int[list.size()];
        this.direction = new int[list.size()];
        for (int i = 0; i < value.length; ++i) {
            value[i]     = i;
            direction[i] = i == 0 ? 0 : -1;
        }

        exhausted = false;

    }

    /** {@inheritDoc} */
    @Override
    public boolean hasNext() {
        return !exhausted;
    }

    /** {@inheritDoc} */
    @Override
    public List<T> next() {

        if (exhausted) {
            throw new NoSuchElementException();
        }

        // the value that will be returned at the end
        final List<T> current = new ArrayList<>(permuted);

        // select element to swap for next permutation
        int selectedIndex     = -1;
        int selectedValue     = -1;
        int selectedDirection = -1;
        for (int i = 0; i < value.length; ++i) {
            if (direction[i] != 0 && value[i] > selectedValue) {
                selectedIndex     = i;
                selectedValue     = value[i];
                selectedDirection = direction[i];
            }
        }
        if (selectedIndex < 0) {
            exhausted = true;
        } else {

            // prepare next permutation

            // swap selected and peer elements
            final int selectedPeer   = selectedIndex + selectedDirection;
            final T tmp              = permuted.get(selectedIndex);
            permuted.set(selectedIndex, permuted.get(selectedPeer));
            permuted.set(selectedPeer, tmp);
            value[selectedIndex]     = value[selectedPeer];
            value[selectedPeer]      = selectedValue;
            direction[selectedIndex] = direction[selectedPeer];
            if (selectedPeer == 0 || selectedPeer == permuted.size() - 1 ||
                value[selectedPeer + selectedDirection] > selectedValue) {
                // we cannot move anymore
                direction[selectedPeer] = 0;
            } else {
                // we will continue moving in the same direction
                direction[selectedPeer] = selectedDirection;
            }

            // enable motion again for greater elements, towards selected one
            for (int i = 0; i < selectedIndex; ++i) {
                if (value[i] > selectedValue) {
                    direction[i] = +1;
                }
            }
            for (int i = selectedIndex + 1; i < value.length; ++i) {
                if (value[i] > selectedValue) {
                    direction[i] = -1;
                }
            }

        }

        return current;

    }

}
