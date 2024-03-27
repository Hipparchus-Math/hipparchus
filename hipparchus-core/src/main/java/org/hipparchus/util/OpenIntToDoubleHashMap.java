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

package org.hipparchus.util;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.Serializable;
import java.util.Arrays;
import java.util.ConcurrentModificationException;
import java.util.NoSuchElementException;

/**
 * Open addressed map from int to double.
 * <p>This class provides a dedicated map from integers to doubles with a
 * much smaller memory overhead than standard <code>java.util.Map</code>.</p>
 * <p>This class is not synchronized. The specialized iterators returned by
 * {@link #iterator()} are fail-fast: they throw a
 * <code>ConcurrentModificationException</code> when they detect the map has been
 * modified during iteration.</p>
 */
public class OpenIntToDoubleHashMap extends AbstractOpenIntHashMap implements Serializable {

    /** Serializable version identifier */
    private static final long serialVersionUID = 20240326L;

    /** Values table. */
    private double[] values;

    /** Return value for missing entries. */
    private final double missingEntries;

    /**
     * Build an empty map with default size and using NaN for missing entries.
     */
    public OpenIntToDoubleHashMap() {
        this(DEFAULT_EXPECTED_SIZE, Double.NaN);
    }

    /**
     * Build an empty map with default size
     * @param missingEntries value to return when a missing entry is fetched
     */
    public OpenIntToDoubleHashMap(final double missingEntries) {
        this(DEFAULT_EXPECTED_SIZE, missingEntries);
    }

    /**
     * Build an empty map with specified size and using NaN for missing entries.
     * @param expectedSize expected number of elements in the map
     */
    public OpenIntToDoubleHashMap(final int expectedSize) {
        this(expectedSize, Double.NaN);
    }

    /**
     * Build an empty map with specified size.
     * @param expectedSize expected number of elements in the map
     * @param missingEntries value to return when a missing entry is fetched
     */
    public OpenIntToDoubleHashMap(final int expectedSize,
                                  final double missingEntries) {
        super(expectedSize);
        values = new double[getCapacity()];
        this.missingEntries = missingEntries;
    }

    /**
     * Copy constructor.
     * @param source map to copy
     */
    public OpenIntToDoubleHashMap(final OpenIntToDoubleHashMap source) {
        super(source);
        values = new double[getCapacity()];
        System.arraycopy(source.values, 0, values, 0, getCapacity());
        missingEntries = source.missingEntries;
    }

    /**
     * Get the stored value associated with the given key
     * @param key key associated with the data
     * @return data associated with the key
     */
    public double get(final int key) {
        final int index = locate(key);
        return index < 0 ? missingEntries : values[index];
    }

    /**
     * Get an iterator over map elements.
     * <p>The specialized iterators returned are fail-fast: they throw a
     * <code>ConcurrentModificationException</code> when they detect the map
     * has been modified during iteration.</p>
     * @return iterator over the map elements
     */
    public Iterator iterator() {
        return new Iterator();
    }

    /** {@inheritDoc} */
    @Override
    public boolean equals(final Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        final OpenIntToDoubleHashMap that = (OpenIntToDoubleHashMap) o;
        return equalKeys(that) && equalStates(that) && Arrays.equals(values, that.values);
    }

    /** {@inheritDoc} */
    @Override
    public int hashCode() {
        return keysStatesHashCode() + Arrays.hashCode(values);
    }

    /**
     * Remove the value associated with a key.
     * @param key key to which the value is associated
     * @return removed value
     */
    public double remove(final int key) {
        final int index = locate(key);
        if (index < 0) {
            return missingEntries;
        } else {
            final double previous = values[index];
            doRemove(index);
            values[index] = missingEntries;
            return previous;
        }
    }

    /**
     * Put a value associated with a key in the map.
     * @param key key to which value is associated
     * @param value value to put in the map
     * @return previous value associated with the key
     */
    public double put(final int key, final double value) {
        final InsertionHolder ih = put(key);
        final double previous = ih.isExisting() ? values[ih.getIndex()] : missingEntries;
        values[ih.getIndex()] = value;
        return previous;
    }

    /** {@inheritDoc} */
    @Override
    protected int growTable(final int oldIndex) {
        final double[] newValues = new double[RESIZE_MULTIPLIER * values.length];
        final int      newIndex  = doGrowTable(oldIndex, (src, dest) -> newValues[dest] = values[src]);
        values = newValues;
        return newIndex;
    }

    /** Iterator class for the map. */
    public class Iterator extends BaseIterator {

        /** Get the value of current entry.
         * @return value of current entry
         * @exception ConcurrentModificationException if the map is modified during iteration
         * @exception NoSuchElementException if there is no element left in the map
         */
        public double value() throws ConcurrentModificationException, NoSuchElementException {
            return values[getCurrent()];
        }

    }

    /**
     * Read a serialized object.
     * @param stream input stream
     * @throws IOException if object cannot be read
     * @throws ClassNotFoundException if the class corresponding
     * to the serialized object cannot be found
     */
    private void readObject(final ObjectInputStream stream)
        throws IOException, ClassNotFoundException {
        stream.defaultReadObject();
        resetCount();
    }

    /**
     * Replace the instance with a data transfer object for serialization.
     * @return data transfer object that will be serialized
     */
    private Object writeReplace() {
        return new DataTransferObject(missingEntries, getSize(), iterator());
    }

    /** Internal class used only for serialization. */
    private static class DataTransferObject implements Serializable {

        /** Serializable UID. */
        private static final long serialVersionUID = 20240326L;

        /** Return value for missing entries. */
        private final double missingEntries;

        /** Keys table. */
        private final int[] keys;

        /** Values table. */
        private final double[] values;

        /** Simple constructor.
         * @param missingEntries return value for missing entries
         * @param size number of objects in the map
         * @param iterator iterator on serialized map
         */
        DataTransferObject(final double missingEntries, final int size, final Iterator iterator) {
            this.missingEntries = missingEntries;
            this.keys           = new int[size];
            this.values         = new double[size];
            for (int i = 0; i < size; ++i) {
                iterator.advance();
                keys[i]   = iterator.key();
                values[i] = iterator.value();
            }
        }

        /** Replace the deserialized data transfer object with a {@link OpenIntToDoubleHashMap}.
         * @return replacement {@link OpenIntToDoubleHashMap}
         */
        private Object readResolve() {
            final OpenIntToDoubleHashMap map = new OpenIntToDoubleHashMap(missingEntries);
            for (int i = 0; i < keys.length; ++i) {
                map.put(keys[i], values[i]);
            }
            return map;
        }

    }

}
