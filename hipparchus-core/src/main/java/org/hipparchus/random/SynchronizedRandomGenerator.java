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
package org.hipparchus.random;

/**
 * Any {@link RandomGenerator} implementation can be thread-safe if it
 * is used through an instance of this class.
 * This is achieved by enclosing calls to the methods of the actual
 * generator inside the overridden {@code synchronized} methods of this
 * class.
 */
public class SynchronizedRandomGenerator implements RandomGenerator {

    /** Object to which all calls will be delegated. */
    private final RandomGenerator wrapped;

    /**
     * Creates a synchronized wrapper for the given {@code RandomGenerator}
     * instance.
     *
     * @param rng Generator whose methods will be called through
     * their corresponding overridden synchronized version.
     * To ensure thread-safety, the wrapped generator <em>must</em>
     * not be used directly.
     */
    public SynchronizedRandomGenerator(RandomGenerator rng) {
        wrapped = rng;
    }

    /** {@inheritDoc} */
    @Override
    public void setSeed(int seed) {
        synchronized (wrapped) {
            wrapped.setSeed(seed);
        }
    }

    /** {@inheritDoc} */
    @Override
    public void setSeed(int[] seed) {
        synchronized (wrapped) {
            wrapped.setSeed(seed);
        }
    }

    /** {@inheritDoc} */
    @Override
    public void setSeed(long seed) {
        synchronized (wrapped) {
            wrapped.setSeed(seed);
        }
    }

    /** {@inheritDoc} */
    @Override
    public void nextBytes(byte[] bytes) {
        synchronized (wrapped) {
            wrapped.nextBytes(bytes);
        }
    }

    /** {@inheritDoc} */
    @Override
    public void nextBytes(byte[] bytes, int offset, int len) {
        synchronized (wrapped) {
            wrapped.nextBytes(bytes, offset, len);
        }
    }

    /** {@inheritDoc} */
    @Override
    public int nextInt() {
        synchronized (wrapped) {
            return wrapped.nextInt();
        }
    }

    /** {@inheritDoc} */
    @Override
    public int nextInt(int n) {
        synchronized (wrapped) {
            return wrapped.nextInt(n);
        }
    }

    /** {@inheritDoc} */
    @Override
    public long nextLong() {
        synchronized (wrapped) {
            return wrapped.nextLong();
        }
    }

    /** {@inheritDoc} */
    @Override
    public long nextLong(long n) {
        synchronized (wrapped) {
            return wrapped.nextLong(n);
        }
    }

    /** {@inheritDoc} */
    @Override
    public boolean nextBoolean() {
        synchronized (wrapped) {
            return wrapped.nextBoolean();
        }
    }

    /** {@inheritDoc} */
    @Override
    public float nextFloat() {
        synchronized (wrapped) {
            return wrapped.nextFloat();
        }
    }

    /** {@inheritDoc} */
    @Override
    public double nextDouble() {
        synchronized (wrapped) {
            return wrapped.nextDouble();
        }
    }

    /** {@inheritDoc} */
    @Override
    public double nextGaussian() {
        synchronized (wrapped) {
            return wrapped.nextGaussian();
        }
    }

}
