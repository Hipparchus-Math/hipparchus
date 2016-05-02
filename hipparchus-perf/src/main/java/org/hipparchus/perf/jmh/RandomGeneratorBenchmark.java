/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements. See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache license, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License. You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the license for the specific language governing permissions and
 * limitations under the license.
 */

package org.hipparchus.perf.jmh;

import java.util.Random;

import org.hipparchus.random.ISAACRandom;
import org.hipparchus.random.JDKRandomGenerator;
import org.hipparchus.random.MersenneTwister;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.Well1024a;
import org.hipparchus.random.Well19937a;
import org.hipparchus.random.Well19937c;
import org.hipparchus.random.Well44497a;
import org.hipparchus.random.Well44497b;
import org.hipparchus.random.Well512a;
import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.Level;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;

/**
 * Tests performance of various RandomGenerator implementations.
 */
// ============================== HOW TO RUN THIS TEST: ====================================
//
// single thread:
// java -jar perf/target/benchmarks.jar ".*Random.*" -i 10 -f 1 -wi 5 -bm sample -tu ns
//
// multiple threads (for example, 4 threads):
// java -jar perf/target/benchmarks.jar ".*Random.*" -i 10 -f 1 -wi 5 -bm sample -tu ns -t 4
//
// Usage help:
// java -jar perf/target/benchmarks.jar -help
//
@State(Scope.Thread)
public class RandomGeneratorBenchmark {

    Random jdkRandom;
    RandomGenerator jdkRandomGenerator;
    RandomGenerator mersenneTwister;
    RandomGenerator isaac;
    RandomGenerator well1024a;
    RandomGenerator well19937a;
    RandomGenerator well19937c;
    RandomGenerator well44497a;
    RandomGenerator well44497b;
    RandomGenerator well512a;


    @Setup(Level.Trial)
    public void up() {
        jdkRandom          = new Random();
        jdkRandomGenerator = new JDKRandomGenerator();
        mersenneTwister    = new MersenneTwister();
        isaac              = new ISAACRandom();
        well1024a          = new Well1024a();
        well19937a         = new Well19937a();
        well19937c         = new Well19937c();
        well44497a         = new Well44497a();
        well44497b         = new Well44497b();
        well512a           = new Well512a();
    }

    @Benchmark
    public void baseline() {
    }

    @Benchmark
    public int jdkRandom() {
        return jdkRandom.nextInt();
    }

    @Benchmark
    public int jdkRandomGenerator() {
        return jdkRandomGenerator.nextInt();
    }

    @Benchmark
    public long mersenneTwister() {
        return mersenneTwister.nextInt();
    }

    @Benchmark
    public long isaac() {
        return isaac.nextInt();
    }

    @Benchmark
    public long well1024a() {
        return well1024a.nextInt();
    }

    @Benchmark
    public long well19937a() {
        return well19937a.nextInt();
    }

    @Benchmark
    public long well19937c() {
        return well19937c.nextInt();
    }

    @Benchmark
    public long well44497a() {
        return well44497a.nextInt();
    }

    @Benchmark
    public long well44497b() {
        return well44497b.nextInt();
    }

    @Benchmark
    public long well512a() {
        return well512a.nextInt();
    }

}
