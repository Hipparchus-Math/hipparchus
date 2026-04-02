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
package org.hipparchus.optim.nonlinear.vector.constrained;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.Locale;
import java.util.stream.Collectors;

import org.hipparchus.util.FastMath;
import org.hipparchus.optim.nonlinear.scalar.ObjectiveFunction;
import org.junit.jupiter.api.DynamicTest;
import org.junit.jupiter.api.TestFactory;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

class MarosMeszarosQPSolverTest {

    /** Fixed benchmark directory provided by user. */
    private static final Path FIXED_BENCHMARK_DIRECTORY =
            Paths.get("org.hipparchus.optim.nonlinear.vector.constrained.maros_meszaros");

    private static final String VERBOSE_PROPERTY = "maros.meszaros.verbose";

    @TestFactory
    List<DynamicTest> testBenchmarkQpsProblemsWhenConfigured() {
        final Path benchmarkDir = resolveBenchmarkDirectory();
        final List<Path> problemFiles = listProblemFilesUnchecked(benchmarkDir);
        return problemFiles.stream()
                       .map(file -> DynamicTest.dynamicTest("benchmark: " + stripExtension(file.getFileName().toString()),
                               () -> runSingleFile(file)))
                       .collect(Collectors.toList());
    }

    private Path resolveBenchmarkDirectory() {
        final List<Path> candidates = java.util.Arrays.asList(
                FIXED_BENCHMARK_DIRECTORY,
                Paths.get("src", "test", "resources")
                     .resolve(FIXED_BENCHMARK_DIRECTORY),
                Paths.get("hipparchus-optim", "src", "test", "resources")
                     .resolve(FIXED_BENCHMARK_DIRECTORY),
                Paths.get("org", "hipparchus", "optim", "nonlinear", "vector", "constrained", "maros_meszaros"),
                Paths.get("src", "test", "resources", "org", "hipparchus", "optim", "nonlinear", "vector", "constrained", "maros_meszaros"),
                Paths.get("hipparchus-optim", "src", "test", "resources", "org", "hipparchus", "optim", "nonlinear", "vector", "constrained", "maros_meszaros")
        );

        for (Path candidate : candidates) {
            if (Files.isDirectory(candidate)) {
                return candidate;
            }
        }

        fail("QPS benchmark directory not found. Expected folder: " +
                FIXED_BENCHMARK_DIRECTORY +
                " (or same under src/test/resources / hipparchus-optim/src/test/resources)");
        return null;
    }

    private List<Path> listProblemFilesUnchecked(final Path directory) {
        try {
            final List<Path> files = Files.walk(directory)
                                          .filter(Files::isRegularFile)
                                          .filter(this::isQpsFile)
                                          .sorted()
                                          .collect(Collectors.toList());
            assertFalse(files.isEmpty(), "No QPS problem files found in benchmark directory: " + directory);
            return files;
        } catch (IOException e) {
            throw new IllegalStateException(e);
        }
    }

    private boolean isQpsFile(final Path path) {
        final String name = path.getFileName().toString().toLowerCase(Locale.ROOT);
        return name.endsWith(".qps") || name.endsWith(".sif") || name.endsWith(".mps");
    }

    private String stripExtension(final String name) {
        final int idx = name.lastIndexOf('.');
        return idx > 0 ? name.substring(0, idx) : name;
    }

    private void runSingleFile(final Path file) {
        final MarosMeszarosDenseProblemLoader loader = new MarosMeszarosDenseProblemLoader();
        final MarosMeszarosDenseProblemLoader.DenseQPProblem problem;
        try {
            problem = loader.loadDenseProblem(file);
        } catch (IOException e) {
            throw new IllegalStateException("Failed to load QPS file: " + file, e);
        }

        final boolean verbose = Boolean.parseBoolean(System.getProperty(VERBOSE_PROPERTY, "true"));
        if (verbose) {
            System.out.println("[MarosMeszarosQPSolverTest] START " + problem.getName() +
                               " file=" + file.getFileName());
        }

        if (problem.isTooLargeForDense()) {
            if (verbose) {
                System.out.println("[MarosMeszarosQPSolverTest] SKIP " + problem.getName() +
                                   " reason=Problem too large for dense solver materialization: " +
                                   problem.getName() + " n=" + problem.getVariableCount() +
                                   ". Increase -Dmaros.meszaros.max.dense.variables if you want to run it.");
            }
            return;
        }

        try {
            final LagrangeSolution solution = runSingle(problem);
            if (verbose) {
                System.out.println("[MarosMeszarosQPSolverTest] PASS " + problem.getName() +
                                   " expectedValue=" + problem.getExpectedValue() +
                                   " actualValue=" + solution.getValue());
            }
        } catch (AssertionError | RuntimeException ex) {
            if (verbose) {
                System.out.println("[MarosMeszarosQPSolverTest] FAIL " + problem.getName() +
                                   " reason=" + ex.getMessage());
            }
            throw ex;
        }
    }

    private LagrangeSolution runSingle(final MarosMeszarosDenseProblemLoader.DenseQPProblem problem) {
        final QPDualActiveSolverR solver = new QPDualActiveSolverR();
        final QuadraticFunction q = new QuadraticFunction(problem.getH(), problem.getC(), problem.getConstantTerm());
        final ObjectiveFunction objectiveFunction = new ObjectiveFunction(q);

        final LagrangeSolution solution = solve(solver, objectiveFunction, problem);

        if (!problem.hasExpectedValue()) {
            fail("Missing expected objective value for problem " + problem.getName() +
                 ". Check 00readme.qp OPT parsing.");
        }

        if (Boolean.parseBoolean(System.getProperty(VERBOSE_PROPERTY, "true"))) {
            System.out.println("[MarosMeszarosQPSolverTest] " + problem.getName() +
                               " expectedValue=" + problem.getExpectedValue() +
                               " actualValue=" + solution.getValue());
        }

//        if (problem.hasExpectedX()) {
//            assertArrayEquals(problem.getExpectedX(),
//                              solution.getX().toArray(),
//                              problem.getXTolerance(),
//                              "Unexpected x for problem " + problem.getName());
//        }
        final double fExpected = problem.getExpectedValue();
        final double absTol = problem.getValueTolerance();
        final double tolF = absTol * (FastMath.abs(fExpected) + 1.0);
        if (fExpected != 0.0) {
            assertTrue(FastMath.abs(solution.getValue() - fExpected) < absTol  * FastMath.abs(fExpected), problem.getName()+";expected:"+fExpected+"--value:"+solution.getValue());
         } else {
            assertTrue(FastMath.abs(solution.getValue()) < absTol,problem.getName()+";expected:"+fExpected+"--value:"+solution.getValue());
}
       
        return solution;
    }

    private LagrangeSolution solve(final QPDualActiveSolverR solver,
                                   final ObjectiveFunction objectiveFunction,
                                   final MarosMeszarosDenseProblemLoader.DenseQPProblem problem) {
        final boolean hasEq = problem.getAeq().length > 0;
        final boolean hasIq = problem.getAiq().length > 0;

        final LinearEqualityConstraint eq = hasEq ?
                new LinearEqualityConstraint(problem.getAeq(), problem.getBeq()) :
                null;
        final LinearInequalityConstraint iq = hasIq ?
                new LinearInequalityConstraint(problem.getAiq(), problem.getBiq()) :
                null;

        if (hasEq && hasIq) {
            return solver.optimize(objectiveFunction, eq, iq);
        } else if (hasEq) {
            return solver.optimize(objectiveFunction, eq);
        } else if (hasIq) {
            return solver.optimize(objectiveFunction, iq);
        }
        return solver.optimize(objectiveFunction);
    }
}