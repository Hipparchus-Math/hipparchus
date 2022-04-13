/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
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

/** Plotter for complex functions using gnuplot.
 */
package org.hipparchus.samples.complex;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import org.hipparchus.analysis.integration.IterativeLegendreGaussIntegrator;
import org.hipparchus.complex.Complex;
import org.hipparchus.complex.ComplexUnivariateIntegrator;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.special.elliptic.jacobi.FieldJacobiElliptic;
import org.hipparchus.special.elliptic.jacobi.JacobiEllipticBuilder;
import org.hipparchus.special.elliptic.legendre.LegendreEllipticIntegral;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.RyuDouble;

public class GnuplotComplexPlotter {

    /** Elliptic integrals characteristic. */
    private Complex n;

    /** Elliptic integrals parameter. */
    private Complex m;

    /** Jacobi functions computer. */
    private FieldJacobiElliptic<Complex> jacobi;

    /** Functions to plot. */
    private List<Predefined> functions;

    /** Output directory (display to terminal if null). */
    private File output;

    /** Domain coloring. */
    private DomainColoring coloring;

    /** Output width. */
    private int width;

    /** Output height. */
    private int height;

    /** Min x. */
    private double xMin;

    /** Max x. */
    private double xMax;

    /** Min y. */
    private double yMin;

    /** Max y. */
    private double yMax;

    /** Maximum number of integrands evaluations for each integral evaluation. */
    private int maxEval;

    /** Integrator for numerically integrated elliptical integrals. */
    final ComplexUnivariateIntegrator integrator;

    /** Default constructor.
     */
    private GnuplotComplexPlotter() {
        n          = new Complex(3.4, 1.3);
        m          = new Complex(0.64);
        jacobi     = JacobiEllipticBuilder.build(m);
        functions  = new ArrayList<>();
        output     = null;
        width      = 800;
        height     = 800;
        coloring   = new SawToothPhaseModuleValue(1.0, 0.7, 1.0, 15);
        xMin       = -7;
        xMax       = +7;
        yMin       = -7;
        yMax       = +7;
        maxEval    = 100000;
        integrator = new ComplexUnivariateIntegrator(new IterativeLegendreGaussIntegrator(24,
                                                                                          1.0e-6,
                                                                                          1.0e-6));
    }

    /** Main program.
     * @param args program arguments
     */
    public static void main(String[] args) {
        final GnuplotComplexPlotter plotter = new GnuplotComplexPlotter();
        try {
            for (int i = 0; i < args.length; ++i) {
                switch (args[i]) {
                    case "--help" :
                        usage(0);
                        break;
                    case "--output-dir" :
                        plotter.output = new File(args[++i]);
                        if (!(plotter.output.exists() && plotter.output.isDirectory() && plotter.output.canWrite())) {
                            System.err.format(Locale.US, "cannot generate output file in %s%n",
                                              plotter.output.getAbsolutePath());
                            System.exit(1);
                        }
                        break;
                    case "--width" :
                        plotter.width = Integer.parseInt(args[++i]);
                        break;
                    case "--height" :
                        plotter.height = Integer.parseInt(args[++i]);
                        break;
                    case "--color" :
                        switch (args[++i]) {
                            case "classical" :
                                plotter.coloring = new ContinuousModuleValue(1.0);
                                break;
                            case "enhanced-module" :
                                plotter.coloring = new SawToothModuleValue(1.0);
                                break;
                            case "enhanced-phase-module" :
                                plotter.coloring = new SawToothPhaseModuleValue(1.0, 0.7, 1.0, 15);
                                break;
                            default :
                                usage(1);
                        }
                        break;
                    case "--xmin" :
                        plotter.xMin = Double.parseDouble(args[++i]);
                        break;
                    case "--xmax" :
                        plotter.xMax = Double.parseDouble(args[++i]);
                        break;
                    case "--ymin" :
                        plotter.yMin = Double.parseDouble(args[++i]);
                        break;
                    case "--ymax" :
                        plotter.yMax = Double.parseDouble(args[++i]);
                        break;
                    case "--m" : {
                        plotter.m      = new Complex(Double.parseDouble(args[++i]), Double.parseDouble(args[++i]));
                        plotter.jacobi = JacobiEllipticBuilder.build(plotter.m);
                        break;
                    }
                    case "--n" : {
                        plotter.n      = new Complex(Double.parseDouble(args[++i]), Double.parseDouble(args[++i]));
                        break;
                    }
                    case "--maxeval" : {
                        plotter.maxEval = Integer.parseInt(args[++i]);
                        break;
                    }
                    case "--function" :
                        try {
                            plotter.functions.add(Predefined.valueOf(args[++i]));
                        } catch (IllegalArgumentException iae) {
                            System.err.format(Locale.US, "unknown function %s, known functions:%n", args[i]);
                            for (final Predefined predefined : Predefined.values()) {
                                System.err.format(Locale.US, " %s", predefined.name());
                            }
                            System.err.format(Locale.US, "%n");
                            System.exit(1);
                        }
                        break;
                    default :
                        usage(1);
                }
            }
            if (plotter.functions.isEmpty()) {
                usage(1);
            }

            plotter.plot();

        } catch (IndexOutOfBoundsException iobe) {
            usage(1);
        } catch (IOException ioe) {
            System.err.println(ioe.getLocalizedMessage());
            System.exit(1);
        }

        System.exit(0);

    }

    /** Display usage.
     * @param status exit code
     */
    private static void usage(final int status) {
        System.err.println("usage: java org.hipparchus.samples.complex.GnuplotComplexPlotter" +
                           " [--help]" +
                           " [--output-dir directory]" +
                           " [--color {classical|enhanced-module|enhanced-phase-module}]" +
                           " [--xmin xMin] [--xmax xMax] [--ymin yMin] [--ymax yMax]" +
                           " [--m mRe mIm] [--n nRe nIm] [--maxeval maxEval]" +
                           " --function {id|sn|cn|dn|cs|...|sin|cos|...} [--function ...]");
        System.exit(status);
        
    }

    /** Plot the function.
     * @throws IOException if gnuplot process cannot be run
     */
    public void plot() throws IOException {
        for (final Predefined predefined : functions) {
            final ProcessBuilder pb = new ProcessBuilder("gnuplot").
                            redirectOutput(ProcessBuilder.Redirect.INHERIT).
                            redirectError(ProcessBuilder.Redirect.INHERIT);
            pb.environment().remove("XDG_SESSION_TYPE");
            final Process gnuplot = pb.start();
            try (PrintStream out = new PrintStream(gnuplot.getOutputStream(), false, StandardCharsets.UTF_8.name())) {
                if (output == null) {
                    out.format(Locale.US, "set terminal qt size %d, %d title 'complex plotter'%n", width, height);
                } else {
                    out.format(Locale.US, "set terminal pngcairo size %d, %d%n", width, height);
                    out.format(Locale.US, "set output '%s'%n", new File(output, predefined.name() + ".png").getAbsolutePath());
                }
                out.format(Locale.US, "set view map scale 1%n");
                out.format(Locale.US, "set xrange [%f : %f]%n", xMin, xMax);
                out.format(Locale.US, "set yrange [%f : %f]%n", yMin, yMax);
                out.format(Locale.US, "set xlabel 'Re(z)'%n");
                out.format(Locale.US, "set ylabel 'Im(z)'%n");
                out.format(Locale.US, "set key off%n");
                out.format(Locale.US, "unset colorbox%n");
                out.format(Locale.US, "set title '%s'%n", predefined.title(m));
                out.format(Locale.US, "$data <<EOD%n");

                for (int i = 0; i < width; ++i) {
                    final double x = xMin + i * (xMax - xMin) / (width - 1);
                    for (int j = 0; j < height; ++j) {
                        final double y = yMin + j * (yMax - yMin) / (height - 1);
                        Complex z  = Complex.valueOf(x, y);
                        Complex fz;
                        try {
                            fz = predefined.evaluator.value(this, z);
                        } catch (MathIllegalStateException e) {
                            fz = Complex.NaN;
                        }
                        out.format(Locale.US, "%12.9f %12.9f %12.9f %12.9f %12.9f%n",
                                   z.getRealPart(), z.getImaginaryPart(),
                                   coloring.hue(fz), coloring.saturation(fz), coloring.value(fz));
                    }
                    out.format(Locale.US, "%n");
                }
                out.format(Locale.US, "EOD%n");
                out.println("splot $data using 1:2:(hsv2rgb($3,$4,$5)) with pm3d lc rgb variable");
                if (output == null) {
                    out.format(Locale.US, "pause mouse close%n");
                } else {
                    System.out.format(Locale.US, "output written to %s%n",
                                      new File(output, predefined.name() + ".png").getAbsolutePath());
                }
            }
        }
    }

    /** Interface for evaluating complex functions. */
    private static interface Evaluator {

        /** Evaluate complex function.
         * @param plotter associated plotter
         * @param z free variable
         */
        Complex value(GnuplotComplexPlotter plotter, Complex z);

    }

    /** Predefined complex functions for plotting. */
    private static enum Predefined {

        id((plotter, z)             -> z),
        sn((plotter, z)             -> plotter.jacobi.valuesN(z).sn()),
        cn((plotter, z)             -> plotter.jacobi.valuesN(z).cn()),
        dn((plotter, z)             -> plotter.jacobi.valuesN(z).dn()),
        cs((plotter, z)             -> plotter.jacobi.valuesS(z).cs()),
        ds((plotter, z)             -> plotter.jacobi.valuesS(z).ds()),
        ns((plotter, z)             -> plotter.jacobi.valuesS(z).ns()),
        dc((plotter, z)             -> plotter.jacobi.valuesC(z).dc()),
        nc((plotter, z)             -> plotter.jacobi.valuesC(z).nc()),
        sc((plotter, z)             -> plotter.jacobi.valuesC(z).sc()),
        nd((plotter, z)             -> plotter.jacobi.valuesD(z).nd()),
        sd((plotter, z)             -> plotter.jacobi.valuesD(z).sd()),
        cd((plotter, z)             -> plotter.jacobi.valuesD(z).cd()),
        arcsn((plotter, z)          -> plotter.jacobi.arcsn(z)),
        arccn((plotter, z)          -> plotter.jacobi.arccn(z)),
        arcdn((plotter, z)          -> plotter.jacobi.arcdn(z)),
        arccs((plotter, z)          -> plotter.jacobi.arccs(z)),
        arcds((plotter, z)          -> plotter.jacobi.arcds(z)),
        arcns((plotter, z)          -> plotter.jacobi.arcns(z)),
        arcdc((plotter, z)          -> plotter.jacobi.arcdc(z)),
        arcnc((plotter, z)          -> plotter.jacobi.arcnc(z)),
        arcsc((plotter, z)          -> plotter.jacobi.arcsc(z)),
        arcnd((plotter, z)          -> plotter.jacobi.arcnd(z)),
        arcsd((plotter, z)          -> plotter.jacobi.arcsd(z)),
        arccd((plotter, z)          -> plotter.jacobi.arccd(z)),
        K((plotter, z)              -> LegendreEllipticIntegral.bigK(z)),
        KPrime((plotter, z)         -> LegendreEllipticIntegral.bigKPrime(z)),
        Fzm((plotter, z)            -> LegendreEllipticIntegral.bigF(z, plotter.m)),
        integratedFzm((plotter, z)  -> LegendreEllipticIntegral.bigF(z, plotter.m, plotter.integrator, plotter.maxEval)),
        E((plotter, z)              -> LegendreEllipticIntegral.bigE(z)),
        Ezm((plotter, z)            -> LegendreEllipticIntegral.bigE(z, plotter.m)),
        integratedEzm((plotter, z)  -> LegendreEllipticIntegral.bigE(z, plotter.m, plotter.integrator, plotter.maxEval)),
        Pi((plotter, z)             -> LegendreEllipticIntegral.bigPi(plotter.n, z)),
        Pizm((plotter, z)           -> LegendreEllipticIntegral.bigPi(plotter.n, z, plotter.m)),
        integratedPizm((plotter, z) -> LegendreEllipticIntegral.bigPi(plotter.n, z, plotter.m, plotter.integrator, plotter.maxEval)),
        sin((plotter, z)            -> FastMath.sin(z)),
        cos((plotter, z)            -> FastMath.cos(z)),
        tan((plotter, z)            -> FastMath.tan(z)),
        asin((plotter, z)           -> FastMath.asin(z)),
        acos((plotter, z)           -> FastMath.acos(z)),
        atan((plotter, z)           -> FastMath.atan(z)),
        sinh((plotter, z)           -> FastMath.sinh(z)),
        cosh((plotter, z)           -> FastMath.cosh(z)),
        tanh((plotter, z)           -> FastMath.tanh(z)),
        asinh((plotter, z)          -> FastMath.asinh(z)),
        acosh((plotter, z)          -> FastMath.acosh(z)),
        atanh((plotter, z)          -> FastMath.atanh(z));
        
        /** Function evaluator. */
        private final Evaluator evaluator;

        /** Simple constructor.
         * @param evaluator function evaluator
         */
        private Predefined(final Evaluator evaluator) {
            this.evaluator = evaluator;
        }

        /** Get plot title.
         * @param m elliptic parameter
         * @return plot title
         */
        public String title(final Complex m) {
            if (name().endsWith("zm")) {
                return name().substring(0, name().length() - 2) +
                       "(z, m = " +
                       RyuDouble.doubleToString(m.getRealPart()) +
                       (m.getImaginary() >= 0 ? " + " : " - ") +
                       RyuDouble.doubleToString(FastMath.abs(m.getImaginaryPart())) +
                       "i)";
            } else {
                return name() + "(z)";
            }
        }

    }

}
