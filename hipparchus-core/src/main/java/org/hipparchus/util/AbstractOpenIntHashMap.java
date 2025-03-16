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

import java.util.Arrays;
import java.util.ConcurrentModificationException;
import java.util.NoSuchElementException;

/** Base class for open addressed map from int.
 * @since 3.1
 */
public abstract class AbstractOpenIntHashMap {

    /** Default starting size.
     * <p>This must be a power of two for bit mask to work properly. </p>
     */
    protected static final int DEFAULT_EXPECTED_SIZE = 16;

    /** Multiplier for size growth when map fills up.
     * <p>This must be a power of two for bit mask to work properly. </p>
     */
    protected static final int RESIZE_MULTIPLIER = 2;

    /** Status indicator for free table entries. */
    private static final byte FREE    = 0;

    /** Status indicator for full table entries. */
    private static final byte FULL    = 1;

    /** Status indicator for removed table entries. */
    private static final byte REMOVED = 2;

    /** Load factor for the map. */
    private static final float LOAD_FACTOR = 0.5f;

    /** Number of bits to perturb the index when probing for collision resolution. */
    private static final int PERTURB_SHIFT = 5;

    /** Keys table. */
    private int[] keys;

    /** States table. */
    private byte[] states;

    /** Current size of the map. */
    private int size;

    /** Bit mask for hash values. */
    private int mask;

    /** Modifications count. */
    private transient int count;

    /** Build an empty map with default size.
     */
    protected AbstractOpenIntHashMap() {
        this(DEFAULT_EXPECTED_SIZE);
    }

    /**
     * Build an empty map with specified size.
     * @param expectedSize expected number of elements in the map
     */
    protected AbstractOpenIntHashMap(final int expectedSize) {
        final int capacity = computeCapacity(expectedSize);
        keys   = new int[capacity];
        states = new byte[capacity];
        mask   = capacity - 1;
        resetCount();
    }

    /**
     * Copy constructor.
     * @param source map to copy
     */
    protected AbstractOpenIntHashMap(final AbstractOpenIntHashMap source) {
        final int length = source.keys.length;
        keys = new int[length];
        System.arraycopy(source.keys, 0, keys, 0, length);
        states = new byte[length];
        System.arraycopy(source.states, 0, states, 0, length);
        size  = source.size;
        mask  = source.mask;
        count = source.count;
    }

    /** Get capacity.
     * @return capacity
     * @since 3.1
     */
    protected int getCapacity() {
        return keys.length;
    }

    /** Get the number of elements stored in the map.
     * @return number of elements stored in the map
     */
    public int getSize() {
        return size;
    }

    /**
     * Compute the capacity needed for a given size.
     * @param expectedSize expected size of the map
     * @return capacity to use for the specified size
     */
    private static int computeCapacity(final int expectedSize) {
        if (expectedSize == 0) {
            return 1;
        }
        final int capacity   = (int) FastMath.ceil(expectedSize / LOAD_FACTOR);
        final int powerOfTwo = Integer.highestOneBit(capacity);
        if (powerOfTwo == capacity) {
            return capacity;
        }
        return nextPowerOfTwo(capacity);
    }

    /** Reset count.
     * @since 3.1
     */
    protected final void resetCount() {
        count = 0;
    }

    /**
     * Find the smallest power of two greater than the input value
     * @param i input value
     * @return smallest power of two greater than the input value
     */
    private static int nextPowerOfTwo(final int i) {
        return Integer.highestOneBit(i) << 1;
    }

    /**
     * Check if a value is associated with a key.
     * @param key key to check
     * @return true if a value is associated with key
     */
    public boolean containsKey(final int key) {

        final int hash  = hashOf(key);
        int index = hash & mask;
        if (containsKey(key, index)) {
            return true;
        }

        if (states[index] == FREE) {
            return false;
        }

        int j = index;
        for (int perturb = perturb(hash); states[index] != FREE; perturb >>= PERTURB_SHIFT) {
            j = probe(perturb, j);
            index = j & mask;
            if (containsKey(key, index)) {
                return true;
            }
        }

        return false;

    }

    /**
     * Perturb the hash for starting probing.
     * @param hash initial hash
     * @return perturbed hash
     */
    private static int perturb(final int hash) {
        return hash & 0x7fffffff;
    }

    /**
     * Find the index at which a key should be inserted
     * @param keys keys table
     * @param states states table
     * @param key key to lookup
     * @param mask bit mask for hash values
     * @return index at which key should be inserted
     */
    private static int findInsertionIndex(final int[] keys, final byte[] states,
                                          final int key, final int mask) {
        final int hash = hashOf(key);
        int index = hash & mask;
        if (states[index] == FREE) {
            return index;
        } else if (states[index] == FULL && keys[index] == key) {
            return changeIndexSign(index);
        }

        int perturb = perturb(hash);
        int j = index;
        if (states[index] == FULL) {
            while (true) {
                j = probe(perturb, j);
                index = j & mask;
                perturb >>= PERTURB_SHIFT;

                if (states[index] != FULL || keys[index] == key) {
                    break;
                }
            }
        }

        if (states[index] == FREE) {
            return index;
        } else if (states[index] == FULL) {
            // due to the loop exit condition,
            // if (states[index] == FULL) then keys[index] == key
            return changeIndexSign(index);
        }

        final int firstRemoved = index;
        while (true) {
            j = probe(perturb, j);
            index = j & mask;

            if (states[index] == FREE) {
                return firstRemoved;
            } else if (states[index] == FULL && keys[index] == key) {
                return changeIndexSign(index);
            }

            perturb >>= PERTURB_SHIFT;

        }

    }

    /**
     * Compute next probe for collision resolution
     * @param perturb perturbed hash
     * @param j previous probe
     * @return next probe
     */
    private static int probe(final int perturb, final int j) {
        return (j << 2) + j + perturb + 1;
    }

    /**
     * Change the index sign
     * @param index initial index
     * @return changed index
     */
    private static int changeIndexSign(final int index) {
        return -index - 1;
    }

    /**
     * Get the number of elements stored in the map.
     * @return number of elements stored in the map
     */
    public int size() {
        return size;
    }

    /**
     * Check if the tables contain an element associated with specified key
     * at specified index.
     * @param key key to check
     * @param index index to check
     * @return true if an element is associated with key at index
     */
    public boolean containsKey(final int key, final int index) {
        return (key != 0 || states[index] == FULL) && keys[index] == key;
    }

    /** Locate the index of value associated with the given key
     * @param key key associated with the data
     * @return index of value associated with the given key or negative
     * if key not present
     */
    protected int locate(final int key) {

        final int hash  = hashOf(key);
        int index = hash & mask;
        if (containsKey(key, index)) {
            return index;
        }

        if (states[index] == FREE) {
            return -1;
        }

        int j = index;
        for (int perturb = perturb(hash); states[index] != FREE; perturb >>= PERTURB_SHIFT) {
            j = probe(perturb, j);
            index = j & mask;
            if (containsKey(key, index)) {
                return index;
            }
        }

        return -1;

    }

    /** Remove an element at specified index.
     * @param index index of the element to remove
     */
    protected void doRemove(int index) {
        keys[index]   = 0;
        states[index] = REMOVED;
        --size;
        ++count;
    }

    /** Put a value associated with a key in the map.
     * @param key key to which value is associated
     * @return holder to manage insertion
     */
    protected InsertionHolder put(final int key) {
        int     oldIndex   = findInsertionIndex(keys, states, key, mask);
        int     newIndex   = oldIndex;
        boolean existing   = false;
        boolean newMapping = true;
        if (oldIndex < 0) {
            oldIndex   = changeIndexSign(oldIndex);
            existing   = true;
            newMapping = false;
        }
        keys[oldIndex] = key;
        states[oldIndex] = FULL;
        if (newMapping) {
            ++size;
            if (shouldGrowTable()) {
                newIndex = growTable(oldIndex);
            }
            ++count;
        }
        return new InsertionHolder(existing ? oldIndex : newIndex, existing);

    }

    /** Grow the tables.
     * @param oldIndex index the entry being inserted should have used
     * @return index the entry being inserted should really use
     */
    protected abstract int growTable(int oldIndex);

    /** Grow the tables.
     * @param oldIndex index the entry being inserted should have used
     * @param valueCopier copier for existing values
     * @return index the entry being inserted should really use
     */
    protected int doGrowTable(final int oldIndex, final ValueCopier valueCopier) {

        int newIndex = oldIndex;
        final int    oldLength = states.length;
        final int[]  oldKeys   = keys;
        final byte[] oldStates = states;

        final int    newLength = RESIZE_MULTIPLIER * oldLength;
        final int[]  newKeys   = new int[newLength];
        final byte[] newStates = new byte[newLength];
        final int newMask = newLength - 1;
        for (int srcIndex = 0; srcIndex < oldLength; ++srcIndex) {
            if (oldStates[srcIndex] == FULL) {
                final int key   = oldKeys[srcIndex];
                final int dstIndex = findInsertionIndex(newKeys, newStates, key, newMask);
                newKeys[dstIndex]  = key;
                valueCopier.copyValue(srcIndex, dstIndex);
                if (srcIndex == oldIndex) {
                    newIndex = dstIndex;
                }
                newStates[dstIndex] = FULL;
            }
        }

        mask   = newMask;
        keys   = newKeys;
        states = newStates;

        return newIndex;

    }

    /**
     * Check if tables should grow due to increased size.
     * @return true if  tables should grow
     */
    private boolean shouldGrowTable() {
        return size > (mask + 1) * LOAD_FACTOR;
    }

    /**
     * Compute the hash value of a key
     * @param key key to hash
     * @return hash value of the key
     */
    private static int hashOf(final int key) {
        final int h = key ^ ((key >>> 20) ^ (key >>> 12));
        return h ^ (h >>> 7) ^ (h >>> 4);
    }

    /** Check if keys are equals.
     * @param other other map
     * @return true if keys are equals
     */
    protected boolean equalKeys(final AbstractOpenIntHashMap other) {
        return Arrays.equals(keys, other.keys);
    }

    /** Check if states are equals.
     * @param other other map
     * @return true if states are equals
     */
    protected boolean equalStates(final AbstractOpenIntHashMap other) {
        return Arrays.equals(states, other.states);
    }

    /** Compute partial hashcode on keys and states.
     * @return partial hashcode on keys and states
     */
    protected int keysStatesHashCode() {
        return  53 * Arrays.hashCode(keys) + 31 * Arrays.hashCode(states);
    }

    /** Iterator class for the map. */
    protected class BaseIterator {

        /** Reference modification count. */
        private final int referenceCount;

        /** Index of current element. */
        private int current;

        /** Index of next element. */
        private int next;

        /**
         * Simple constructor.
         */
        protected BaseIterator() {

            // preserve the modification count of the map to detect concurrent modifications later
            referenceCount = count;

            // initialize current index
            next = -1;
            try {
                advance();
            } catch (NoSuchElementException nsee) { // NOPMD
                // ignored
            }

        }

        /**
         * Check if there is a next element in the map.
         * @return true if there is a next element
         */
        public boolean hasNext() {
            return next >= 0;
        }

        /** Get index of current entry.
         * @return key of current entry
         * @since 3.1
         */
        protected int getCurrent() {
            return current;
        }

        /**
         * Get the key of current entry.
         * @return key of current entry
         * @exception ConcurrentModificationException if the map is modified during iteration
         * @exception NoSuchElementException if there is no element left in the map
         */
        public int key() throws ConcurrentModificationException, NoSuchElementException {
            return keys[getCurrent()];
        }

        /**
         * Advance iterator one step further.
         * @exception ConcurrentModificationException if the map is modified during iteration
         * @exception NoSuchElementException if there is no element left in the map
         */
        public void advance()
            throws ConcurrentModificationException, NoSuchElementException {

            if (referenceCount != count) {
                throw new ConcurrentModificationException();
            }

            // advance on step
            current = next;

            // prepare next step
            try {
                while (states[++next] != FULL) { // NOPMD
                    // nothing to do
                }
            } catch (ArrayIndexOutOfBoundsException e) {
                next = -2;
                if (current < 0) {
                    throw new NoSuchElementException(); // NOPMD
                }
            }

        }

    }

    /** Holder for handling values insertion.
     * @since 3.1
     */
    protected static class InsertionHolder {

        /** Index at which new value should be put. */
        private final int index;

        /** Indicator for value already present before insertion. */
        private final boolean existing;

        /** Simple constructor.
         * @param index index at which new value should be put
         * @param existing indicator for value already present before insertion
         */
        InsertionHolder(final int index, final boolean existing) {
            this.index    = index;
            this.existing = existing;
        }

        /** Get index at which new value should be put.
         * @return index at which new value should be put
         */
        public int getIndex() {
            return index;
        }

        /** Get indicator for value already present before insertion.
         * @return indicator for value already present before insertion
         */
        public boolean isExisting() {
            return existing;
        }
    }

    /** Interface for copying values.
     * @since 3.1
     */
    @FunctionalInterface
    protected interface ValueCopier {
        /** Copy a value.
         * @param src source index
         * @param dest destination index
         */
        void copyValue(int src, int dest);
    }

}
