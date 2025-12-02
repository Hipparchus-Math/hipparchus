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

import org.hipparchus.util.Precision;

/**
 * Utility class for formatting SQP iteration logs with dynamic precision and aligned columns.
 */
public class SQPLogger {

    /** LS column fixed to 2 digits + 1 space for safety. */
    private static final int LS_WIDTH = 3;

    /** Field start. */
    private static final String FIELD_START = " %";

    /** Field continuation. */
    private static final String FIELD_CONTINUATION = "s |";

    /** Fields width. */
    private int width;

    /** Format of header line. */
    private String headerFormat;

    /** Format of row lines. */
    private String rowFormat;

    /** Debug printer. */
    private DebugPrinter printer;

    /**
     * Constructs a LogFormatter with custom epsilon.
     *
     * @param epsilon convergence threshold
     */
    public SQPLogger(double epsilon) {
        setEps(epsilon);
    }

    /**
     * Updates the formatter precision and formats based on a new epsilon value.
     *
     * @param epsilon convergence threshold
     */
    public void setEps(double epsilon) {
        final int precision = (int) Math.ceil(-Math.log10(epsilon)) + 2;
        this.width = precision + 7; // space for digits, sign, exponent, padding

        final String f = "%%-%ds";
        String col = String.format(f, width);
        String lsCol = String.format(f, LS_WIDTH);
        this.headerFormat = String.format(
            "[SQP] ITER %%2s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
            col, lsCol, col, col, col, col, col, col, col
        );

        final String percent = "%%";
        String fld = String.format(percent + width + "." + precision + "f");
        String intf = String.format(percent + LS_WIDTH + "d");
        this.rowFormat = String.format(
            "[SQP] ITER %%2d | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
            fld, intf, fld, fld, fld, fld, fld, fld, fld
        );
    }

    /** Set debug printer.
     * @param debugPrinter debug printer
     */
    public void setDebugPrinter(DebugPrinter debugPrinter) {
        this.printer = debugPrinter;
    }

    /** Get header line.
     * @return header line
     */
    public String header() {
        return String.format(headerFormat,
            "", "alpha", "LS", "dxNorm", "dx'Hdx", "KKT", "viol", "sigma", "penalty", "f(x)");
    }

    /** Format one row.
     * @param iter     iteration number
     * @param alpha    step length
     * @param lsCount  line search iteration
     * @param dxNorm   || dX ||
     * @param dxHdx    dX H dX
     * @param kkt      Lagrangian norm
     * @param viol     constraints violations
     * @param sigma    solution of the additional variable in QP subproblem
     * @param penalty  penalty
     * @param fx       objective function evaluation
     * @return formatted row
     */
    public String formatRow(final int iter, final double alpha, final int lsCount,
                            final double dxNorm, final double dxHdx, final double kkt,
                            final double viol, final double sigma, final double penalty, final double fx) {
        return String.format(rowFormat,
                             iter, alpha, lsCount, dxNorm, dxHdx, kkt, viol, sigma, penalty, fx);
    }

    /** Log header.
     */
    public void logHeader() {
        if (printer != null) {
            printer.print(header());
        }
    }

    /**
     * Log one row.
     * @param crit2 norm criterion
     * @param crit1 gradient criterion?
     * @param crit0 Lagrangian norm criterion
     * @param crit3 constraints violations criterion
     */
    public void logRow(final boolean crit2, final boolean crit1, final boolean crit0, final boolean crit3) {
        if (printer == null) {
            return;
        }
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("[SQP] ITER %2d |", -1)).
           append(String.format(FIELD_START + width + FIELD_CONTINUATION, "")).
           append(String.format(FIELD_START + LS_WIDTH + FIELD_CONTINUATION, "")).
           append(String.format(FIELD_START + width + FIELD_CONTINUATION, crit2)).
           append(String.format(FIELD_START + width + FIELD_CONTINUATION, crit1)).
           append(String.format(FIELD_START + width + FIELD_CONTINUATION, crit0)).
           append(String.format(FIELD_START + width + FIELD_CONTINUATION, crit3)).
           append(String.format(FIELD_START + width + FIELD_CONTINUATION, "")).
           append(String.format(FIELD_START + width + FIELD_CONTINUATION, "")).
           append(String.format(FIELD_START + width + FIELD_CONTINUATION, ""));
        printer.print(sb.toString());
    }

    /** Log one row.
     * @param iter     iteration number
     * @param alpha    step length
     * @param lsCount  line search iteration
     * @param dxNorm   || dX ||
     * @param dxHdx    dX H dX
     * @param kkt      Lagrangian norm
     * @param viol     constraints violations
     * @param sigma    solution of the additional variable in QP subproblem
     * @param penalty  penalty
     * @param fx       objective function evaluation
     */
    public void logRow(int iter, double alpha, int lsCount,
                       double dxNorm, double dxHdx, double kkt,
                       double viol, double sigma, double penalty, double fx) {
        if (printer != null) {
            printer.print(formatRow(iter, alpha, lsCount, dxNorm, dxHdx, kkt, viol, sigma, penalty, fx));
        }
    }

    /** Get default logger.
     * @return default logger
     */
    public static SQPLogger defaultLogger() {
        return new SQPLogger(Precision.EPSILON);
    }

}
