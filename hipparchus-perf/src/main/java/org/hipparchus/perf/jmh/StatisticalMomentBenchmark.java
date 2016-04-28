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

import java.lang.reflect.Constructor;

import org.hipparchus.stat.descriptive.StorelessUnivariateStatistic;
import org.hipparchus.stat.descriptive.moment.Mean;
import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.Level;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;

/**
 * Tests performance of statistical moments.
 * The test checks whether using an immutable Mean instance is still slower
 * than using a FirstMoment directly. Previously with mutable instances
 * it was slower as the jvm could not optimize as aggressive as it can with
 * final instances.
 */
// ============================== HOW TO RUN THIS TEST: ====================================
//
// single thread:
// java -jar perf/target/benchmarks.jar ".*StatisticalMoment.*" -i 10 -f 1 -wi 5 -bm sample -tu ns
//
// multiple threads (for example, 4 threads):
// java -jar perf/target/benchmarks.jar ".*StatisticalMoment.*" -i 10 -f 1 -wi 5 -bm sample -tu ns -t 4
//
// Usage help:
// java -jar perf/target/benchmarks.jar -help
//
@State(Scope.Thread)
public class StatisticalMomentBenchmark {

    StorelessUnivariateStatistic firstMoment;
    StorelessUnivariateStatistic meanNormal;
    StorelessUnivariateStatistic meanOptimized;

    double value = 2.0;

    @Setup(Level.Trial)
    public void up() {
        try {
            // need to get a FirstMoment instance via reflection, as it is package private
            Class<?> firstMomentClass = Class.forName("org.hipparchus.stat.descriptive.moment.FirstMoment");
            Constructor<?> ctor = firstMomentClass.getDeclaredConstructor();
            ctor.setAccessible(true);
            firstMoment = (StorelessUnivariateStatistic) ctor.newInstance();
        } catch (Exception ex) {}
        meanNormal = new Mean();
    }

    @Benchmark
    public double baseline() {
        return value;
    }

    @Benchmark
    public double firstMoment() {
        firstMoment.increment(value);
        return firstMoment.getResult();
    }

    @Benchmark
    public double meanNormal() {
        meanNormal.increment(value);
        return meanNormal.getResult();
    }

}
