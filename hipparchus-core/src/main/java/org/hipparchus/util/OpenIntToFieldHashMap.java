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
import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.ConcurrentModificationException;
import java.util.NoSuchElementException;

import org.hipparchus.Field;
import org.hipparchus.FieldElement;

/**
 * Open addressed map from int to FieldElement.
 * <p>This class provides a dedicated map from integers to FieldElements with a
 * much smaller memory overhead than standard <code>java.util.Map</code>.</p>
 * <p>This class is not synchronized. The specialized iterators returned by
 * {@link #iterator()} are fail-fast: they throw a
 * <code>ConcurrentModificationException</code> when they detect the map has been
 * modified during iteration.</p>
 * @param <T> the type of the field elements
 */
public class OpenIntToFieldHashMap<T extends FieldElement<T>> extends AbstractOpenIntHashMap implements Serializable {

    /** Serializable version identifier. */
    private static final long serialVersionUID = 20240326L;

    /** Field to which the elements belong. */
    private final Field<T> field;

    /** Values table. */
    private T[] values;

    /** Return value for missing entries. */
    private final T missingEntries;

    /**
     * Build an empty map with default size and using zero for missing entries.
     * @param field field to which the elements belong
     */
    public OpenIntToFieldHashMap(final Field<T>field) {
        this(field, DEFAULT_EXPECTED_SIZE, field.getZero());
    }

    /**
     * Build an empty map with default size
     * @param field field to which the elements belong
     * @param missingEntries value to return when a missing entry is fetched
     */
    public OpenIntToFieldHashMap(final Field<T>field, final T missingEntries) {
        this(field,DEFAULT_EXPECTED_SIZE, missingEntries);
    }

    /**
     * Build an empty map with specified size and using zero for missing entries.
     * @param field field to which the elements belong
     * @param expectedSize expected number of elements in the map
     */
    public OpenIntToFieldHashMap(final Field<T> field,final int expectedSize) {
        this(field,expectedSize, field.getZero());
    }

    /**
     * Build an empty map with specified size.
     * @param field field to which the elements belong
     * @param expectedSize expected number of elements in the map
     * @param missingEntries value to return when a missing entry is fetched
     */
    public OpenIntToFieldHashMap(final Field<T> field,final int expectedSize,
                                  final T missingEntries) {
        super(expectedSize);
        this.field = field;
        final int capacity = computeCapacity(expectedSize);
        values = buildArray(capacity);
        this.missingEntries = missingEntries;
    }

    /**
     * Copy constructor.
     * @param source map to copy
     */
    public OpenIntToFieldHashMap(final OpenIntToFieldHashMap<T> source) {
        super(source);
        field = source.field;
        values = buildArray(getCapacity());
        System.arraycopy(source.values, 0, values, 0, getCapacity());
        missingEntries = source.missingEntries;
    }

    /**
     * Get the stored value associated with the given key
     * @param key key associated with the data
     * @return data associated with the key
     */
    public T get(final int key) {

        final int hash  = hashOf(key);
        int index = hash & getMask();
        if (containsKey(key, index)) {
            return values[index];
        }

        if (getState(index) == FREE) {
            return missingEntries;
        }

        int j = index;
        for (int perturb = perturb(hash); getState(index) != FREE; perturb >>= PERTURB_SHIFT) {
            j = probe(perturb, j);
            index = j & getMask();
            if (containsKey(key, index)) {
                return values[index];
            }
        }

        return missingEntries;

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
        @SuppressWarnings("unchecked")
        final OpenIntToFieldHashMap<T> that = (OpenIntToFieldHashMap<T>) o;
        return equalKeys(that) &&
               equalStates(that) &&
               field.equals(that.field) &&
               Arrays.equals(values, that.values);
    }

    /** {@inheritDoc} */
    @Override
    public int hashCode() {
        return keysStatesHashCode() + 67 * field.hashCode() + Arrays.hashCode(values);
    }

    /**
     * Remove the value associated with a key.
     * @param key key to which the value is associated
     * @return removed value
     */
    public T remove(final int key) {

        final int hash  = hashOf(key);
        int index = hash & getMask();
        if (containsKey(key, index)) {
            return doRemove(index);
        }

        if (getState(index) == FREE) {
            return missingEntries;
        }

        int j = index;
        for (int perturb = perturb(hash); getState(index) != FREE; perturb >>= PERTURB_SHIFT) {
            j = probe(perturb, j);
            index = j & getMask();
            if (containsKey(key, index)) {
                return doRemove(index);
            }
        }

        return missingEntries;

    }

    /**
     * Remove an element at specified index.
     * @param index index of the element to remove
     * @return removed value
     */
    private T doRemove(int index) {
        topDoRemove(index);
        final T previous = values[index];
        values[index] = missingEntries;
        return previous;
    }

    /**
     * Put a value associated with a key in the map.
     * @param key key to which value is associated
     * @param value value to put in the map
     * @return previous value associated with the key
     */
    public T put(final int key, final T value) {
        final InsertionHolder ih = put(key);
        final T previous = ih.isExisting() ? values[ih.getIndex()] : missingEntries;
        values[ih.getIndex()] = value;
        return previous;
    }

    /**
     * Grow the tables.
     */
    @Override
    protected int growTable(final int oldIndex) {
        final T[] newValues = buildArray(RESIZE_MULTIPLIER * values.length);
        final int newIndex  = doGrowTable(oldIndex, (src, dest) -> newValues[dest] = values[src]);
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
        public T value() throws ConcurrentModificationException, NoSuchElementException {
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

    /** Build an array of elements.
     * @param length size of the array to build
     * @return a new array
     */
    @SuppressWarnings("unchecked") // field is of type T
    private T[] buildArray(final int length) {
        return (T[]) Array.newInstance(field.getRuntimeClass(), length);
    }

}
