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

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Loader for Maros-Meszaros QP benchmark files in text QPS format.
 */
class MarosMeszarosDenseProblemLoader {

    /** Cached verbose flag so it is available across all parsing helper methods. */
    private final boolean verbose;

    /** Property for minimum accepted density. */
    static final String MIN_DENSITY_PROPERTY = "maros.meszaros.min.density";

    /** Property for maximum variable count (inclusive). */
    static final String MAX_VARIABLES_PROPERTY = "maros.meszaros.max.variables";

    /** Safety cap to avoid heap exhaustion when building dense matrices. */
    static final String HARD_MAX_VARIABLES_PROPERTY = "maros.meszaros.hard.max.variables";

    /** Verbose trace property for benchmark loading. */
    static final String VERBOSE_PROPERTY = "maros.meszaros.verbose";

    /** Enable strict README-vs-parse validation checks. */
    static final String STRICT_METADATA_VALIDATION_PROPERTY = "maros.meszaros.strict.metadata.validation";

    /** Default variable limit to avoid loading huge benchmarks by default. */
    private static final double DEFAULT_MAX_VARIABLES = 1000.0;

    /** Max variables allowed for dense materialization (safety). */
    static final String MAX_DENSE_VARIABLES_PROPERTY = "maros.meszaros.max.dense.variables";

    /** Default dense materialization limit. */
    private static final int DEFAULT_MAX_DENSE_VARIABLES = 5000;

    /** Tolerance to detect equality after range expansion. */
    private static final double EQUALITY_TOL = 1.0e-12;

    MarosMeszarosDenseProblemLoader() {
        this.verbose = Boolean.parseBoolean(System.getProperty(VERBOSE_PROPERTY, "true"));
    }


    /**
     * Load dense QP problems from a benchmark directory containing {@code .qps} files.
     *
     * @param directory benchmark directory
     * @return parsed dense QP problems
     * @throws IOException if reading fails
     */
    List<DenseQPProblem> loadDenseProblems(final Path directory) throws IOException {
        final Double minDensity = parseOptionalBound(System.getProperty(MIN_DENSITY_PROPERTY));
        final Double configuredMaxVariables = parseOptionalBound(System.getProperty(MAX_VARIABLES_PROPERTY));
        final Double maxVariables = configuredMaxVariables == null ? DEFAULT_MAX_VARIABLES : configuredMaxVariables;
        final Double configuredHardMaxVariables = parseOptionalBound(System.getProperty(HARD_MAX_VARIABLES_PROPERTY));
        final Double hardMaxVariables = configuredHardMaxVariables == null ? maxVariables : configuredHardMaxVariables;
        final Double configuredMaxDenseVariables = parseOptionalBound(System.getProperty(MAX_DENSE_VARIABLES_PROPERTY));
        final int maxDenseVariables = configuredMaxDenseVariables == null ?
                DEFAULT_MAX_DENSE_VARIABLES : configuredMaxDenseVariables.intValue();

        final List<Path> files = Files.walk(directory)
                                      .filter(Files::isRegularFile)
                                      .filter(this::isQpsFile)
                                      .sorted(Comparator.comparing(Path::toString))
                                      .collect(Collectors.toList());

        final Map<String, ProblemMetadata> metadataByName = loadProblemMetadata(directory);
        final List<DenseQPProblem> problems = new ArrayList<>();
        int skippedByHardLimit = 0;

        for (Path file : files) {
            final String problemName = stripExtension(file.getFileName().toString());
            final ProblemMetadata metadata = metadataByName.get(canonicalProblemKey(problemName));

            if (metadata == null && maxVariables != null) {
                final int estimatedVariables = estimateVariableCount(file, maxVariables.intValue() + 1);
                if (estimatedVariables > maxVariables.intValue()) {
                    if (verbose) {
                        log("SKIP (estimated vars from file): " + file + " n~" + estimatedVariables +
                            " > " + maxVariables.intValue());
                    }
                    continue;
                }
            }

            if (hardMaxVariables != null && metadata != null && metadata.n > hardMaxVariables.doubleValue()) {
                skippedByHardLimit++;
                if (verbose) {
                    log("SKIP (hard max vars): " + file + " n=" + metadata.n + " > " + hardMaxVariables.intValue());
                }
                continue;
            }

            if (!matchesMetadataFilters(metadata, maxVariables, minDensity)) {
                if (verbose) {
                    log("SKIP (metadata filter): " + file + describeMetadata(metadata));
                }
                continue;
            }

            if (verbose) {
                log("OPEN: " + file + describeMetadata(metadata));
            }

            // A problem file is always parsed in full; we never partially load a file.
            final DenseQPProblem problem = loadProblem(file, metadata, maxDenseVariables);
            if (matchesVariableCount(problem, maxVariables) && matchesDensity(problem, minDensity)) {
                problems.add(problem);
                if (verbose) {
                    log("LOADED: " + problem.getName() + " n=" + problem.getC().length);
                }
            } else if (verbose) {
                log("SKIP (post-parse filter): " + problem.getName() + " n=" + problem.getC().length);
            }
        }

        if (verbose) {
            log("SUMMARY: loaded=" + problems.size() + ", skippedByHardLimit=" + skippedByHardLimit +
                ", maxVariables=" + maxVariables.intValue() +
                ", hardMax=" + hardMaxVariables.intValue() +
                ", filesFound=" + files.size() + ", maxDenseVariables=" + maxDenseVariables);
        }

        return problems;
    }

    /**
     * Load a single dense QP problem from one benchmark file.
     *
     * @param file benchmark file
     * @return parsed dense QP problem
     * @throws IOException if reading fails
     */
    DenseQPProblem loadDenseProblem(final Path file) throws IOException {
        final Path metadataDirectory = file.getParent() == null ? file.toAbsolutePath().getParent() : file.getParent();
        final Map<String, ProblemMetadata> metadataByName = metadataDirectory == null ?
                java.util.Collections.emptyMap() : loadProblemMetadata(metadataDirectory);

        final String problemName = stripExtension(file.getFileName().toString());
        final ProblemMetadata metadata = metadataByName.get(canonicalProblemKey(problemName));

        final Double configuredMaxVariables = parseOptionalBound(System.getProperty(MAX_VARIABLES_PROPERTY));
        final double maxVariables = configuredMaxVariables == null ? DEFAULT_MAX_VARIABLES : configuredMaxVariables.doubleValue();

        final Double configuredHardMaxVariables = parseOptionalBound(System.getProperty(HARD_MAX_VARIABLES_PROPERTY));
        final double hardMaxVariables = configuredHardMaxVariables == null ? maxVariables : configuredHardMaxVariables.doubleValue();

        final Double configuredMaxDenseVariables = parseOptionalBound(System.getProperty(MAX_DENSE_VARIABLES_PROPERTY));
        final int maxDenseVariables = configuredMaxDenseVariables == null ?
                DEFAULT_MAX_DENSE_VARIABLES : configuredMaxDenseVariables.intValue();

        final int stopAfter = (int) Math.min((double) maxDenseVariables + 1.0,
                                             Math.min(maxVariables, hardMaxVariables) + 1.0);
        final int estimatedVariables = estimateVariableCount(file, stopAfter);
        if (estimatedVariables > maxDenseVariables ||
            estimatedVariables > maxVariables ||
            estimatedVariables > hardMaxVariables) {
            final double expectedValue = metadata == null ? Double.NaN : metadata.opt;
            final int reportedVariables = metadata == null ? estimatedVariables : metadata.n;
            return DenseQPProblem.tooLarge(problemName, reportedVariables, expectedValue);
        }

        if (metadata != null) {
            if (metadata.n > hardMaxVariables || metadata.n > maxVariables || metadata.n > maxDenseVariables) {
                return DenseQPProblem.tooLarge(problemName, metadata.n, metadata.opt);
            }
        }

        return loadProblem(file, metadata, maxDenseVariables);
    }

    private boolean isQpsFile(final Path path) {
        final String name = path.getFileName().toString().toLowerCase(Locale.ROOT);
        return name.endsWith(".qps") || name.endsWith(".sif") || name.endsWith(".mps");
    }

    /**
     * Fast, low-memory estimate of variable count from COLUMNS section.
     *
     * @param file qps file
     * @param stopAfter stop as soon as count reaches this threshold
     * @return number of distinct variables seen (possibly truncated at threshold)
     * @throws IOException if reading fails
     */
    private int estimateVariableCount(final Path file, final int stopAfter) throws IOException {
        boolean inColumns = false;
        final java.util.Set<String> variables = new java.util.HashSet<>();

        try (BufferedReader reader = Files.newBufferedReader(file)) {
            String line;
            while ((line = reader.readLine()) != null) {
                final String rawLine = line;
                final String trimmed = rawLine.trim();
                if (trimmed.isEmpty() || trimmed.startsWith("*") || trimmed.startsWith("&")) {
                    continue;
                }

                final String upper = trimmed.toUpperCase(Locale.ROOT);
                if ("COLUMNS".equals(upper)) {
                    inColumns = true;
                    continue;
                }
                if (inColumns && ("RHS".equals(upper) || "RANGES".equals(upper) ||
                                  "BOUNDS".equals(upper) || "QMATRIX".equals(upper) ||
                                  "QUADOBJ".equals(upper) || "DMATRIX".equals(upper) ||
                                  "ENDATA".equals(upper))) {
                    break;
                }

                if (!inColumns) {
                    continue;
                }

                final String[] tokens = parseDataTokens(rawLine, "COLUMNS");
                if (tokens.length == 0) {
                    continue;
                }
                if (tokens.length >= 2 && "'MARKER'".equalsIgnoreCase(tokens[1])) {
                    continue;
                }

                variables.add(tokens[0]);
                if (variables.size() >= stopAfter) {
                    return variables.size();
                }
            }
        }

        return variables.size();
    }

    private DenseQPProblem loadProblem(final Path file,
                                     final ProblemMetadata metadata,
                                     final int maxDenseVariables) throws IOException {
        final ProblemLayout layout = scanLayout(file);
        final int n = layout.varIndex.size();
        if (n > maxDenseVariables) {
            final double expectedValue = metadata == null ? Double.NaN : metadata.opt;
            return DenseQPProblem.tooLarge(stripExtension(file.getFileName().toString()), n, expectedValue);
        }

        final double[] c = new double[n];
        final double[][] rowA = new double[layout.constraintRowCount][n];
        final double[] rhs = new double[layout.constraintRowCount];
        final boolean[] rhsSeen = new boolean[layout.constraintRowCount];
        final double[] ranges = new double[layout.constraintRowCount];
        final boolean[] rangesSeen = new boolean[layout.constraintRowCount];
        final double[][] lowerQ = new double[n][n];

        final double[] lb = new double[n];
        final double[] ub = new double[n];
        Arrays.fill(lb, 0.0);
        Arrays.fill(ub, Double.POSITIVE_INFINITY);

        final double[] objectiveConstant = new double[] { 0.0 };
        parseAndFill(file, layout, c, rowA, rhs, rhsSeen, ranges, rangesSeen, lowerQ, lb, ub, objectiveConstant);

        final List<double[]> eqRows = new ArrayList<>();
        final List<Double> eqValues = new ArrayList<>();
        final List<double[]> iqRows = new ArrayList<>();
        final List<Double> iqValues = new ArrayList<>();

        for (Map.Entry<String, Integer> entry : layout.constraintRowIndex.entrySet()) {
            final String rowName = entry.getKey();
            final int rowIndex = entry.getValue().intValue();
            final char rowType = layout.rows.get(rowName).charValue();
            final double[] a = rowA[rowIndex];
            final double rhsValue = rhsSeen[rowIndex] ? rhs[rowIndex] : 0.0;

            if (!rangesSeen[rowIndex]) {
                addRowWithoutRange(rowType, a, rhsValue, eqRows, eqValues, iqRows, iqValues);
            } else {
                addRowWithRange(rowType, a, rhsValue, ranges[rowIndex], eqRows, eqValues, iqRows, iqValues);
            }
        }

        final double[][] h = buildSymmetricH(lowerQ, file);

        validateAgainstMetadata(file, metadata, h, layout);

        final double expectedValue = metadata == null ? Double.NaN : metadata.opt;

        return new DenseQPProblem(stripExtension(file.getFileName().toString()),
                                  h,
                                  c,
                                  toMatrix(eqRows, n),
                                  toVector(eqValues),
                                  toMatrix(iqRows, n),
                                  toVector(iqValues),
                                  objectiveConstant[0],
                                  Arrays.copyOf(lb, lb.length),
                                  Arrays.copyOf(ub, ub.length),
                                  null,
                                  expectedValue,
                                  1e-6,
                                  1e-5);
    }

    private void addRowWithoutRange(final char rowType,
                                    final double[] a,
                                    final double rhsValue,
                                    final List<double[]> eqRows,
                                    final List<Double> eqValues,
                                    final List<double[]> iqRows,
                                    final List<Double> iqValues) {
        if (rowType == 'E') {
            eqRows.add(copyRow(a));
            eqValues.add(rhsValue);
        } else if (rowType == 'G') {
            iqRows.add(copyRow(a));
            iqValues.add(rhsValue);
        } else if (rowType == 'L') {
            iqRows.add(negate(a));
            iqValues.add(-rhsValue);
        }
    }

    private void addRowWithRange(final char rowType,
                                 final double[] a,
                                 final double rhsValue,
                                 final double rangeValue,
                                 final List<double[]> eqRows,
                                 final List<Double> eqValues,
                                 final List<double[]> iqRows,
                                 final List<Double> iqValues) {

        final double absRange = Math.abs(rangeValue);
        final double lower;
        final double upper;

        if (rowType == 'G') {
            lower = rhsValue;
            upper = rhsValue + absRange;
        } else if (rowType == 'L') {
            lower = rhsValue - absRange;
            upper = rhsValue;
        } else if (rowType == 'E') {
            if (rangeValue >= 0.0) {
                lower = rhsValue;
                upper = rhsValue + absRange;
            } else {
                lower = rhsValue - absRange;
                upper = rhsValue;
            }
        } else {
            return;
        }

        if (Math.abs(upper - lower) <= EQUALITY_TOL) {
            eqRows.add(copyRow(a));
            eqValues.add(lower);
        } else {
            iqRows.add(copyRow(a));
            iqValues.add(lower);
            iqRows.add(negate(a));
            iqValues.add(-upper);
        }
    }

    private double[][] buildSymmetricH(final double[][] lowerQ, final Path file) {
        final int n = lowerQ.length;
        final double[][] h = new double[n][n];

        for (int i = 0; i < n; i++) {
            h[i][i] = lowerQ[i][i];
            for (int j = 0; j < i; j++) {
                final double lij = lowerQ[i][j];
                final double uji = lowerQ[j][i];
                if (uji != 0.0 && lij != 0.0 && Double.compare(uji, lij) != 0 && verbose) {
                    log("WARN (quadratic duplicate mismatch): using lower-triangular value for (" +
                        i + "," + j + ") lower=" + lij + " upper=" + uji + " file=" + file);
                }
                final double v = (lij != 0.0) ? lij : uji;
                h[i][j] = v;
                h[j][i] = v;
            }
        }

        return h;
    }

    private double[] copyRow(final double[] row) {
        return Arrays.copyOf(row, row.length);
    }

    private ProblemLayout scanLayout(final Path file) throws IOException {
        final ProblemLayout layout = new ProblemLayout();
        String section = "";

        try (BufferedReader reader = Files.newBufferedReader(file)) {
            String line;
            while ((line = reader.readLine()) != null) {
                final String rawLine = line;
                final String trimmed = rawLine.trim();
                if (trimmed.isEmpty() || trimmed.startsWith("*") || trimmed.startsWith("&")) {
                    continue;
                }

                final String nextSection = parseSectionHeader(rawLine, trimmed);
                if (nextSection != null) {
                    section = nextSection;
                    if ("OBJNAME".equals(section)) {
                        final String[] headerTokens = parseDataTokens(rawLine);
                        if (headerTokens.length > 1) {
                            layout.objectiveRow = headerTokens[1];
                        }
                    }
                    if ("ENDATA".equals(section)) {
                        break;
                    }
                    continue;
                }

                if ("ROWS".equals(section) && shouldUseFixedFormat(rawLine, section, parseDataTokens(rawLine))) {
                    layout.fixedFormat = true;
                }

                final String[] tokens = parseDataTokens(rawLine, section, layout.fixedFormat);
                if ("ROWS".equals(section)) {
                    if (tokens.length >= 2) {
                        final char type = Character.toUpperCase(tokens[0].charAt(0));
                        final String name = tokens[1];
                        layout.rows.put(name, Character.valueOf(type));
                        if (type == 'N' && layout.objectiveRow == null) {
                            layout.objectiveRow = name;
                        }
                    }
                } else if ("OBJNAME".equals(section)) {
                    parseObjectiveName(layout, tokens);
                } else if ("COLUMNS".equals(section)) {
                    if (tokens.length > 0 && !(tokens.length >= 2 && "'MARKER'".equalsIgnoreCase(tokens[1]))) {
                        layout.addVariable(tokens[0]);
                    }
                } else if ("BOUNDS".equals(section)) {
                    final int varToken = boundVariableTokenIndex(tokens);
                    if (varToken >= 0) {
                        layout.addVariable(tokens[varToken]);
                    }
                } else if ("QUADOBJ".equals(section) || "QMATRIX".equals(section) || "DMATRIX".equals(section)) {
                    if (tokens.length >= 1) {
                        layout.addVariable(tokens[0]);
                    }
                    for (int i = 1; i + 1 < tokens.length; i += 2) {
                        layout.addVariable(tokens[i]);
                    }
                }
            }
        }

        for (Map.Entry<String, Character> row : layout.rows.entrySet()) {
            final String rowName = row.getKey();
            final char rowType = row.getValue().charValue();
            if (rowName.equals(layout.objectiveRow) || rowType == 'N') {
                continue;
            }
            if (rowType == 'E') {
                layout.eqCount++;
                layout.constraintRowIndex.put(rowName, Integer.valueOf(layout.constraintRowCount++));
            } else if (rowType == 'G' || rowType == 'L') {
                layout.iqCount++;
                layout.constraintRowIndex.put(rowName, Integer.valueOf(layout.constraintRowCount++));
            }
        }

        return layout;
    }

    private void parseAndFill(final Path file,
                              final ProblemLayout layout,
                              final double[] c,
                              final double[][] rowA,
                              final double[] rhs,
                              final boolean[] rhsSeen,
                              final double[] ranges,
                              final boolean[] rangesSeen,
                              final double[][] lowerQ,
                              final double[] lb,
                              final double[] ub,
                              final double[] objectiveConstant) throws IOException {
        String section = "";

        try (BufferedReader reader = Files.newBufferedReader(file)) {
            String line;
            while ((line = reader.readLine()) != null) {
                final String rawLine = line;
                final String trimmed = rawLine.trim();
                if (trimmed.isEmpty() || trimmed.startsWith("*") || trimmed.startsWith("&")) {
                    continue;
                }

                final String nextSection = parseSectionHeader(rawLine, trimmed);
                if (nextSection != null) {
                    section = nextSection;
                    if ("ENDATA".equals(section)) {
                        break;
                    }
                    continue;
                }

                final String[] tokens = parseDataTokens(rawLine, section, layout.fixedFormat);
                if (tokens.length == 0) {
                    continue;
                }

                if ("COLUMNS".equals(section)) {
                    parseColumnsLine(layout, tokens, c, rowA);
                } else if ("RHS".equals(section)) {
                    parseRhsLine(layout, tokens, rhs, rhsSeen, objectiveConstant);
                } else if ("RANGES".equals(section)) {
                    parseRangesLine(layout, tokens, ranges, rangesSeen);
                } else if ("BOUNDS".equals(section)) {
                    parseBoundsLine(layout, tokens, lb, ub);
                } else if ("QUADOBJ".equals(section) || "QMATRIX".equals(section) || "DMATRIX".equals(section)) {
                    parseQuadraticLine(layout, tokens, lowerQ);
                }
            }
        }
    }

    private void parseColumnsLine(final ProblemLayout layout,
                                  final String[] tokens,
                                  final double[] c,
                                  final double[][] rowA) {
        if (tokens.length < 3 || (tokens.length >= 2 && "'MARKER'".equalsIgnoreCase(tokens[1]))) {
            return;
        }

        final Integer varIndex = layout.varIndex.get(tokens[0]);
        if (varIndex == null) {
            return;
        }

        for (int i = 1; i + 1 < tokens.length; i += 2) {
            final String row = tokens[i];
            final Double valueObj = tryParseDouble(tokens[i + 1]);
            if (valueObj == null) {
                continue;
            }
            final double value = valueObj.doubleValue();
            if (row.equals(layout.objectiveRow)) {
                c[varIndex.intValue()] += value;
                continue;
            }
            final Integer rowIndex = layout.constraintRowIndex.get(row);
            if (rowIndex == null) {
                if (verbose) {
                    final Character rowType = layout.rows.get(row);
                    if (rowType != null && rowType.charValue() == 'N') {
                        log("Ignoring rim objective coefficient row='" + row + "' var='" + tokens[0] + "' value=" + value);
                    }
                }
                continue;
            }
            rowA[rowIndex.intValue()][varIndex.intValue()] += value;
        }
    }

    private void parseRhsLine(final ProblemLayout layout,
                              final String[] tokens,
                              final double[] rhs,
                              final boolean[] rhsSeen,
                              final double[] objectiveConstant) {
        for (int i = 1; i + 1 < tokens.length; i += 2) {
            final String row = tokens[i];
            final Double valueObj = tryParseDouble(tokens[i + 1]);
            if (valueObj == null) {
                continue;
            }
            final double value = valueObj.doubleValue();
            if (row.equals(layout.objectiveRow)) {
                objectiveConstant[0] = -value;
                continue;
            }
            final Integer rowIndex = layout.constraintRowIndex.get(row);
            if (rowIndex != null) {
                rhs[rowIndex.intValue()] = value;
                rhsSeen[rowIndex.intValue()] = true;
            }
        }
    }

    private void parseRangesLine(final ProblemLayout layout,
                                 final String[] tokens,
                                 final double[] ranges,
                                 final boolean[] rangesSeen) {
        for (int i = 1; i + 1 < tokens.length; i += 2) {
            final String row = tokens[i];
            final Double valueObj = tryParseDouble(tokens[i + 1]);
            if (valueObj == null) {
                continue;
            }
            final Integer rowIndex = layout.constraintRowIndex.get(row);
            if (rowIndex != null) {
                ranges[rowIndex.intValue()] = valueObj.doubleValue();
                rangesSeen[rowIndex.intValue()] = true;
            }
        }
    }

    private void parseBoundsLine(final ProblemLayout layout,
                                 final String[] tokens,
                                 final double[] lb,
                                 final double[] ub) {
        final int varToken = boundVariableTokenIndex(tokens);
        if (varToken < 0) {
            return;
        }
        final Integer varIndex = layout.varIndex.get(tokens[varToken]);
        if (varIndex == null) {
            return;
        }

        final String type = tokens[0].toUpperCase(Locale.ROOT);
        final int valueToken = boundValueTokenIndex(tokens, type);
        final Double parsed = valueToken >= 0 ? tryParseDouble(tokens[valueToken]) : Double.valueOf(0.0);
        final double value = parsed == null ? 0.0 : parsed.doubleValue();
        final int j = varIndex.intValue();
        if ("LO".equals(type)) {
            lb[j] = value;
        } else if ("UP".equals(type)) {
            ub[j] = value;
        } else if ("FX".equals(type)) {
            lb[j] = value;
            ub[j] = value;
        } else if ("FR".equals(type)) {
            lb[j] = Double.NEGATIVE_INFINITY;
            ub[j] = Double.POSITIVE_INFINITY;
        } else if ("MI".equals(type)) {
            lb[j] = Double.NEGATIVE_INFINITY;
            if (!Double.isFinite(ub[j])) {
                ub[j] = 0.0;
            }
        } else if ("PL".equals(type)) {
            ub[j] = Double.POSITIVE_INFINITY;
        } else if ("BV".equals(type)) {
            lb[j] = 0.0;
            ub[j] = 1.0;
        } else if ("LI".equals(type)) {
            lb[j] = value;
        } else if ("UI".equals(type)) {
            ub[j] = value;
        }
    }


    private int boundVariableTokenIndex(final String[] tokens) {
        if (tokens.length < 2) {
            return -1;
        }
        if (tokens.length == 2) {
            return 1;
        }
        final String type = tokens[0].toUpperCase(Locale.ROOT);
        return tokens.length == 3 && boundTypeRequiresValue(type) ? 1 : 2;
    }

    private int boundValueTokenIndex(final String[] tokens, final String type) {
        if (tokens.length > 3) {
            return 3;
        }
        return tokens.length == 3 && boundTypeRequiresValue(type) ? 2 : -1;
    }

    private boolean boundTypeRequiresValue(final String type) {
        return "LO".equals(type) ||
               "UP".equals(type) ||
               "FX".equals(type) ||
               "LI".equals(type) ||
               "UI".equals(type);
    }

    private void parseObjectiveName(final ProblemLayout layout, final String[] tokens) {
        if (tokens.length > 0) {
            layout.objectiveRow = tokens[0];
        }
    }

    private void parseQuadraticLine(final ProblemLayout layout,
                                    final String[] tokens,
                                    final double[][] lowerQ) {
        if (tokens.length < 3) {
            return;
        }
        final Integer i = layout.varIndex.get(tokens[0]);
        if (i == null) {
            return;
        }

        for (int k = 1; k + 1 < tokens.length; k += 2) {
            final Integer j = layout.varIndex.get(tokens[k]);
            if (j == null) {
                continue;
            }
            final Double valueObj = tryParseDouble(tokens[k + 1]);
            if (valueObj != null) {
                lowerQ[i.intValue()][j.intValue()] += valueObj.doubleValue();
            }
        }
    }

    private void validateAgainstMetadata(final Path file,
                                         final ProblemMetadata metadata,
                                         final double[][] h,
                                         final ProblemLayout layout) {
        if (metadata == null) {
            return;
        }

        final boolean strict = Boolean.parseBoolean(System.getProperty(STRICT_METADATA_VALIDATION_PROPERTY, "false"));
        final int parsedN = h.length;
        final int parsedM = layout.constraintRowCount;

        if (parsedN != metadata.n) {
            final String message = "README mismatch for " + file + ": n(parsed)=" + parsedN +
                    " n(readme)=" + metadata.n;
            if (strict) {
                throw new IllegalArgumentException(message);
            }
            log("WARN " + message);
        }

        if (parsedM != metadata.m) {
            final String message = "README mismatch for " + file + ": m(parsed)=" + parsedM +
                    " m(readme)=" + metadata.m;
            if (strict) {
                throw new IllegalArgumentException(message);
            }
            log("WARN " + message);
        }
    }

    private Map<String, ProblemMetadata> loadProblemMetadata(final Path directory) throws IOException {
        final Path metadataPath = resolveMetadataPath(directory);
        if (metadataPath == null) {
            return java.util.Collections.emptyMap();
        }

        final Map<String, ProblemMetadata> metadataByName = new LinkedHashMap<>();
        boolean inTable = false;

        try (BufferedReader reader = Files.newBufferedReader(metadataPath)) {
            String line;
            while ((line = reader.readLine()) != null) {
                final String trimmed = line.trim();
                if (!inTable) {
                    if (trimmed.startsWith("NAME") && trimmed.contains("M") && trimmed.contains("N") && trimmed.contains("QNZ")) {
                        inTable = true;
                    }
                    continue;
                }

                if (trimmed.isEmpty()) {
                    continue;
                }

                final String[] tokens = trimmed.split("\\s+");
                if (tokens.length < 7) {
                    if (trimmed.startsWith("Acknowledgements") || trimmed.startsWith("-")) {
                        break;
                    }
                    continue;
                }

                final ParsedMetadataRow row = parseMetadataRow(tokens);
                if (row == null) {
                    if (trimmed.startsWith("Acknowledgements") || trimmed.startsWith("-")) {
                        break;
                    }
                    continue;
                }

                metadataByName.put(canonicalProblemKey(row.name),
                        new ProblemMetadata(row.m, row.n, row.nz, row.qn, row.qnz, row.opt));
            }
        }

        return metadataByName;
    }




    private ParsedMetadataRow parseMetadataRow(final String[] tokens) {
        for (int i = 1; i + 5 < tokens.length; i++) {
            final Integer m = tryParseInt(tokens[i]);
            final Integer n = tryParseInt(tokens[i + 1]);
            final Integer nz = tryParseInt(tokens[i + 2]);
            final Integer qn = tryParseInt(tokens[i + 3]);
            final Integer qnz = tryParseInt(tokens[i + 4]);
            final Double opt = tryParseDouble(tokens[i + 5]);

            if (m == null || n == null || nz == null || qn == null || qnz == null || opt == null) {
                continue;
            }

            final StringBuilder nameBuilder = new StringBuilder();
            for (int k = 0; k < i; k++) {
                if (k > 0) {
                    nameBuilder.append(' ');
                }
                nameBuilder.append(tokens[k]);
            }

            final String name = nameBuilder.toString().trim();
            if (name.isEmpty() || !Character.isLetterOrDigit(name.charAt(0))) {
                return null;
            }

            return new ParsedMetadataRow(name, m.intValue(), n.intValue(), nz.intValue(),
                    qn.intValue(), qnz.intValue(), opt.doubleValue());
        }

        return null;
    }

    private Integer tryParseInt(final String value) {
        try {
            return Integer.valueOf(value);
        } catch (NumberFormatException ex) {
            return null;
        }
    }

    private Double tryParseDouble(final String value) {
        try {
            return Double.valueOf(parseDouble(value));
        } catch (NumberFormatException ex) {
            return null;
        }
    }
    private Path resolveMetadataPath(final Path directory) {
        final List<Path> candidates = Arrays.asList(
                directory.resolve("00readme.qp"),
                directory.resolve("00README.QP"),
                directory.getParent() == null ? null : directory.getParent().resolve("00readme.qp"),
                directory.getParent() == null ? null : directory.getParent().resolve("00README.QP")
        );

        for (Path candidate : candidates) {
            if (candidate != null && Files.isRegularFile(candidate)) {
                return candidate;
            }
        }
        return null;
    }

    private String canonicalProblemKey(final String value) {
        final StringBuilder builder = new StringBuilder(value.length());
        for (int i = 0; i < value.length(); i++) {
            final char c = Character.toLowerCase(value.charAt(i));
            if (Character.isLetterOrDigit(c)) {
                builder.append(c);
            }
        }
        return builder.toString();
    }

    private boolean matchesMetadataFilters(final ProblemMetadata metadata,
                                           final Double maxVariables,
                                           final Double minDensity) {
        if (metadata == null) {
            return true;
        }

        if (maxVariables != null && metadata.n > maxVariables.doubleValue()) {
            return false;
        }

        if (minDensity != null) {
            final double estimatedDensity = metadata.estimatedDensity();
            if (estimatedDensity < minDensity.doubleValue()) {
                return false;
            }
        }

        return true;
    }

    private String describeMetadata(final ProblemMetadata metadata) {
        if (metadata == null) {
            return " [no-readme-metadata]";
        }
        return " [n=" + metadata.n + ", m=" + metadata.m + ", nz=" + metadata.nz +
                ", qn=" + metadata.qn + ", qnz=" + metadata.qnz + "]";
    }

    private void log(final String message) {
        System.out.println("[MarosMeszarosLoader] " + message);
    }

    private String parseSectionHeader(final String rawLine,
                                      final String trimmed) {
        if (rawLine.isEmpty() || Character.isWhitespace(rawLine.charAt(0))) {
            return null;
        }
        final String upper = trimmed.toUpperCase(Locale.ROOT);
        final String[] tokens = upper.split("\\s+");
        if (tokens.length == 0) {
            return null;
        }
        final String head = tokens[0];
        if ("OBJECT".equals(head) && tokens.length > 1 && "BOUND".equals(tokens[1])) {
            return "OBJECT BOUND";
        }
        if ("NAME".equals(head) ||
            "OBJSENSE".equals(head) ||
            "OBJNAME".equals(head) ||
            "ROWS".equals(head) ||
            "COLUMNS".equals(head) ||
            "RHS".equals(head) ||
            "RANGES".equals(head) ||
            "BOUNDS".equals(head) ||
            "QUADOBJ".equals(head) ||
            "QMATRIX".equals(head) ||
            "DMATRIX".equals(head) ||
            "ENDATA".equals(head)) {
            return head;
        }
        return null;
    }

    private double parseDouble(final String value) {
        return Double.parseDouble(value.replace('D', 'E').replace('d', 'e'));
    }

    /**
     * Parse one MPS/QPS data line using whitespace tokenization.
     *
     * <p>For Maros-Meszaros files this avoids silent truncation of long names
     * (for example variable names beyond index 999/1000) that can happen with
     * strict fixed-width slicing.</p>
     *
     * @param line raw data line (original line, not trimmed)
     * @return parsed tokens
     */
    private String[] parseDataTokens(final String line) {
        final String trimmed = line.trim();
        if (trimmed.isEmpty()) {
            return new String[0];
        }
        return trimmed.split("\\s+");
    }

    /**
     * Parse one MPS/QPS data line using either whitespace tokenization or fixed
     * MPS columns.
     *
     * <p>Most Maros-Meszaros files use free whitespace-separated names, including
     * names longer than eight characters. Some older fixed-format files, however,
     * use the full eight-character MPS name fields and embed spaces inside names
     * (for example {@code "DEDO3 1R"}). Whitespace tokenization would split such
     * names and silently merge rows/columns. This method keeps the free parser by
     * default, but switches to fixed columns when the free tokens cannot form
     * valid row/value pairs for the current section.</p>
     *
     * @param line raw data line
     * @param section active MPS/QPS section
     * @return parsed tokens
     */
    private String[] parseDataTokens(final String line, final String section) {
        final String[] freeTokens = parseDataTokens(line);
        if (shouldUseFixedFormat(line, section, freeTokens)) {
            return parseFixedDataTokens(line, section);
        }
        return freeTokens;
    }

    private String[] parseDataTokens(final String line,
                                     final String section,
                                     final boolean fixedFormat) {
        if (fixedFormat && isFixedFormatSection(section)) {
            return parseFixedDataTokens(line, section);
        }
        return parseDataTokens(line, section);
    }

    private boolean isFixedFormatSection(final String section) {
        return "ROWS".equals(section) ||
               "COLUMNS".equals(section) ||
               "RHS".equals(section) ||
               "RANGES".equals(section) ||
               "BOUNDS".equals(section) ||
               "QUADOBJ".equals(section) ||
               "QMATRIX".equals(section) ||
               "DMATRIX".equals(section);
    }

    private boolean shouldUseFixedFormat(final String line,
                                         final String section,
                                         final String[] freeTokens) {
        if (line.isEmpty() || !Character.isWhitespace(line.charAt(0)) || section == null) {
            return false;
        }

        if ("ROWS".equals(section)) {
            final String[] fixedTokens = parseFixedDataTokens(line, section);
            return fixedTokens.length == 2 && freeTokens.length != 2;
        }

        if ("BOUNDS".equals(section)) {
            final String[] fixedTokens = parseFixedDataTokens(line, section);
            return fixedTokens.length >= 3 &&
                   fixedTokens[2].indexOf(' ') >= 0 &&
                   (freeTokens.length > 4 || (fixedTokens.length > 3 && freeTokens.length > 3 &&
                    tryParseDouble(freeTokens[3]) == null));
        }

        if ("COLUMNS".equals(section) || "RHS".equals(section) || "RANGES".equals(section) ||
            "QUADOBJ".equals(section) || "QMATRIX".equals(section) || "DMATRIX".equals(section)) {
            final String[] fixedTokens = parseFixedDataTokens(line, section);
            if (fixedTokens.length < 3 || tryParseDouble(fixedTokens[2]) == null) {
                return false;
            }
            if (freeTokens.length < 3 || tryParseDouble(freeTokens[2]) == null) {
                return true;
            }
            return fixedTokens.length > 3 && tryParseDouble(fixedTokens[4]) != null &&
                   (freeTokens.length < 5 || tryParseDouble(freeTokens[4]) == null);
        }

        return false;
    }

    private String[] parseFixedDataTokens(final String line, final String section) {
        final List<String> tokens = new ArrayList<>();

        if ("ROWS".equals(section)) {
            addFixedField(tokens, line, 1, 3);
            addFixedField(tokens, line, 4, 12);
        } else if ("BOUNDS".equals(section)) {
            addFixedField(tokens, line, 1, 3);
            addFixedField(tokens, line, 4, 12);
            addFixedField(tokens, line, 14, 22);
            addFixedField(tokens, line, 24, 36);
        } else {
            addFixedField(tokens, line, 4, 12);
            addFixedField(tokens, line, 14, 22);
            addFixedField(tokens, line, 24, 36);
            addFixedField(tokens, line, 39, 47);
            addFixedField(tokens, line, 49, 61);
        }

        return tokens.toArray(new String[0]);
    }

    private void addFixedField(final List<String> tokens,
                               final String line,
                               final int start,
                               final int end) {
        if (line.length() <= start) {
            return;
        }
        final String token = line.substring(start, Math.min(line.length(), end)).trim();
        if (!token.isEmpty()) {
            tokens.add(token);
        }
    }

    private double[][] toMatrix(final List<double[]> rows, final int cols) {
        final double[][] matrix = new double[rows.size()][cols];
        for (int i = 0; i < rows.size(); i++) {
            matrix[i] = rows.get(i);
        }
        return matrix;
    }

    private double[] toVector(final List<Double> values) {
        final double[] out = new double[values.size()];
        for (int i = 0; i < values.size(); i++) {
            out[i] = values.get(i).doubleValue();
        }
        return out;
    }

    private double[] negate(final double[] row) {
        final double[] out = new double[row.length];
        for (int i = 0; i < row.length; i++) {
            out[i] = -row[i];
        }
        return out;
    }

    boolean matchesVariableCount(final DenseQPProblem problem,
                                 final Double maxVariables) {
        final double variableCount = problem.c.length;
        return maxVariables == null || variableCount <= maxVariables;
    }

    private static Double parseOptionalBound(final String value) {
        if (value == null || value.trim().isEmpty()) {
            return null;
        }
        return Double.valueOf(value.trim());
    }

    boolean matchesDensity(final DenseQPProblem problem, final Double minDensity) {
        final List<Double> densities = new ArrayList<>();
        densities.add(density(problem.h));
        if (problem.aeq.length > 0) {
            densities.add(density(problem.aeq));
        }
        if (problem.aiq.length > 0) {
            densities.add(density(problem.aiq));
        }

        double avg = 0.0;
        for (double d : densities) {
            avg += d;
        }
        avg /= densities.size();

        return minDensity == null || avg >= minDensity;
    }

    private double density(final double[][] matrix) {
        int nonZero = 0;
        int total = 0;
        for (double[] row : matrix) {
            for (double v : row) {
                total++;
                if (Math.abs(v) > 1e-15) {
                    nonZero++;
                }
            }
        }
        return total == 0 ? 0.0 : (double) nonZero / (double) total;
    }

    private String stripExtension(final String name) {
        final int idx = name.lastIndexOf('.');
        return idx > 0 ? name.substring(0, idx) : name;
    }

    private static final class ParsedMetadataRow {
        private final String name;
        private final int m;
        private final int n;
        private final int nz;
        private final int qn;
        private final int qnz;
        private final double opt;

        ParsedMetadataRow(final String name,
                          final int m,
                          final int n,
                          final int nz,
                          final int qn,
                          final int qnz,
                          final double opt) {
            this.name = name;
            this.m = m;
            this.n = n;
            this.nz = nz;
            this.qn = qn;
            this.qnz = qnz;
            this.opt = opt;
        }
    }

    private static final class ProblemLayout {
        private final Map<String, Character> rows = new LinkedHashMap<>();
        private final Map<String, Integer> varIndex = new LinkedHashMap<>();
        private final Map<String, Integer> constraintRowIndex = new LinkedHashMap<>();
        private String objectiveRow;
        private int eqCount;
        private int iqCount;
        private int constraintRowCount;
        private boolean fixedFormat;

        private void addVariable(final String name) {
            if (!varIndex.containsKey(name)) {
                varIndex.put(name, Integer.valueOf(varIndex.size()));
            }
        }
    }

    private static final class ProblemMetadata {
        private final int m;
        private final int n;
        private final int nz;
        private final int qn;
        private final int qnz;
        private final double opt;

        ProblemMetadata(final int m,
                        final int n,
                        final int nz,
                        final int qn,
                        final int qnz,
                        final double opt) {
            this.m = m;
            this.n = n;
            this.nz = nz;
            this.qn = qn;
            this.qnz = qnz;
            this.opt = opt;
        }

        private double estimatedDensity() {
            final double aDensity = m > 0 && n > 0 ? (double) nz / ((double) m * (double) n) : 0.0;
            final double qDiag = qn;
            final double qOffDiag = 2.0 * qnz;
            final double qDensity = n > 0 ? (qDiag + qOffDiag) / ((double) n * (double) n) : 0.0;
            return (aDensity + qDensity) / 2.0;
        }
    }

    /** Container for one parsed dense QP problem. */
    static final class DenseQPProblem {
        private final String name;
        private final double[][] h;
        private final double[] c;
        private final double[][] aeq;
        private final double[] beq;
        private final double[][] aiq;
        private final double[] biq;
        private final double constantTerm;
        private final double[] lowerBound;
        private final double[] upperBound;
        private final double[] expectedX;
        private final double expectedValue;
        private final double xTolerance;
        private final double valueTolerance;
        private final int variableCount;
        private final boolean tooLargeForDense;

        static DenseQPProblem tooLarge(final String name,
                                      final int variableCount,
                                      final double expectedValue) {
            return new DenseQPProblem(name,
                                      new double[0][0],
                                      new double[0],
                                      new double[0][0],
                                      new double[0],
                                      new double[0][0],
                                      new double[0],
                                      0.0,
                                      new double[0],
                                      new double[0],
                                      null,
                                      expectedValue,
                                      1e-6,
                                      1e-5,
                                      variableCount,
                                      true);
        }

        DenseQPProblem(final String name,
                       final double[][] h,
                       final double[] c,
                       final double[][] aeq,
                       final double[] beq,
                       final double[][] aiq,
                       final double[] biq,
                       final double constantTerm,
                       final double[] lowerBound,
                       final double[] upperBound,
                       final double[] expectedX,
                       final double expectedValue,
                       final double xTolerance,
                       final double valueTolerance) {
            this(name, h, c, aeq, beq, aiq, biq, constantTerm, lowerBound, upperBound, expectedX, expectedValue, xTolerance, valueTolerance, c.length, false);
        }

        DenseQPProblem(final String name,
                       final double[][] h,
                       final double[] c,
                       final double[][] aeq,
                       final double[] beq,
                       final double[][] aiq,
                       final double[] biq,
                       final double constantTerm,
                       final double[] lowerBound,
                       final double[] upperBound,
                       final double[] expectedX,
                       final double expectedValue,
                       final double xTolerance,
                       final double valueTolerance,
                       final int variableCount,
                       final boolean tooLargeForDense) {
            this.name = name;
            this.h = h;
            this.c = c;
            this.aeq = aeq;
            this.beq = beq;
            this.aiq = aiq;
            this.biq = biq;
            this.constantTerm = constantTerm;
            this.lowerBound = lowerBound;
            this.upperBound = upperBound;
            this.expectedX = expectedX;
            this.expectedValue = expectedValue;
            this.xTolerance = xTolerance;
            this.valueTolerance = valueTolerance;
            this.variableCount = variableCount;
            this.tooLargeForDense = tooLargeForDense;
        }

        String getName() {
            return name;
        }

        double[][] getH() {
            return h;
        }

        double[] getC() {
            return c;
        }

        double[][] getAeq() {
            return aeq;
        }

        double[] getBeq() {
            return beq;
        }

        double[][] getAiq() {
            return aiq;
        }

        double[] getBiq() {
            return biq;
        }

        double getConstantTerm() {
            return constantTerm;
        }

        double[] getLowerBound() {
            return lowerBound;
        }

        double[] getUpperBound() {
            return upperBound;
        }

        double[] getExpectedX() {
            return expectedX;
        }

        double getExpectedValue() {
            return expectedValue;
        }

        double getXTolerance() {
            return xTolerance;
        }

        double getValueTolerance() {
            return valueTolerance;
        }

        boolean hasExpectedX() {
            return expectedX != null && expectedX.length > 0;
        }

        boolean hasExpectedValue() {
            return !Double.isNaN(expectedValue);
        }

        int getVariableCount() {
            return variableCount;
        }

        boolean isTooLargeForDense() {
            return tooLargeForDense;
        }
    }
}
