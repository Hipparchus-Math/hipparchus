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
import java.util.Locale;

import org.hipparchus.complex.Complex;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalStateException;
import org.hipparchus.special.elliptic.jacobi.FieldJacobiElliptic;
import org.hipparchus.special.elliptic.jacobi.JacobiEllipticBuilder;
import org.hipparchus.special.elliptic.legendre.LegendreEllipticIntegral;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.FieldSinCos;
import org.hipparchus.util.MathUtils;

public class GnuplotComplexPlotter {

    /** Jacobi functions computer. */
    private FieldJacobiElliptic<Complex> jacobi;

    /** Elliptic module. */
    private Complex k;

    /** Nome. */
    private Complex q;

    /** Quarter period. */
    private Complex bigK;

    /** Function to plot. */
    private Predefined f;

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

    /** Default constructor.
     */
    private GnuplotComplexPlotter() {
        setK(new Complex(0.8));
        jacobi   = JacobiEllipticBuilder.build(k.multiply(k));
        f        = Predefined.sn;
        output   = null;
        width    = 800;
        height   = 800;
        coloring = new SawToothPhaseModuleValue(1.0, 0.7, 1.0, 15);
        xMin     = -7;
        xMax     = +7;
        yMin     = -7;
        yMax     = +7;
    }

    /** Set the elliptic module.
     * @param k elliptic module
     */
    private void setK(final Complex k) {
        this.k    = k;
        this.q    = LegendreEllipticIntegral.nome(k);
        this.bigK = LegendreEllipticIntegral.bigK(k);
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
                    case "--k" : {
                        plotter.setK(new Complex(Double.parseDouble(args[++i]), Double.parseDouble(args[++i])));
                        break;
                    }
                    case "--function" :
                        try {
                            plotter.f = Predefined.valueOf(args[++i]);
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
                           " [--xmin xmin] [--xmax xmax] [--ymin ymin] [--ymax ymax]" +
                           " [--k kRe kIm]" +
                           " --function {sn|cn|dn|cs|...|sin|cos|...}");
        System.exit(status);
        
    }

    /** Plot the function.
     * @throws IOException if gnuplot process cannot be run
     */
    public void plot() throws IOException {
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
                out.format(Locale.US, "set output '%s'%n", new File(output, f.name() + ".png").getAbsolutePath());
            }
            out.format(Locale.US, "set view map scale 1%n");
            out.format(Locale.US, "set xrange [%f : %f]%n", xMin, xMax);
            out.format(Locale.US, "set yrange [%f : %f]%n", yMin, yMax);
            out.format(Locale.US, "set xlabel 'Re(z)'%n");
            out.format(Locale.US, "set ylabel 'Im(z)'%n");
            out.format(Locale.US, "set key off%n");
            out.format(Locale.US, "unset colorbox%n");
            out.format(Locale.US, "set title '%s(z)'%n", f.name());
            out.format(Locale.US, "$data <<EOD%n");

            for (int i = 0; i < width; ++i) {
                final double x = xMin + i * (xMax - xMin) / (width - 1);
                for (int j = 0; j < height; ++j) {
                    final double y = yMin + j * (yMax - yMin) / (height - 1);
                    final Complex z  = Complex.valueOf(x, y);
                    final Complex fz = f.evaluator.value(this, z);
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
                                  new File(output, f.name() + ".png").getAbsolutePath());
            }
        }
    }

    /** Evaluate the Jacobi theta functions.
     * @param z argument of the functions
     * @return container for the four Jacobi theta functions θ₁(z|τ), θ₂(z|τ), θ₃(z|τ), and θ₄(z|τ)
     */
    public Complex[] values(final Complex q, final Complex z) {

        Complex qSquare = q.multiply(q);
        Complex qFourth = FastMath.sqrt(FastMath.sqrt(q));
        // the computation is based on Fourier series,
        // see Digital Library of Mathematical Functions section 20.2
        // https://dlmf.nist.gov/20.2

        // base angle for Fourier Series
        final FieldSinCos<Complex> sc1 = FastMath.sinCos(z);

        // recursion rules initialization
        double                       sgn   = 1.0;
        Complex                            qNN   = q.getField().getOne();
        Complex                            qTwoN = q.getField().getOne();
        Complex                            qNNp1 = q.getField().getOne();
        FieldSinCos<Complex> sc2n1 = sc1;
        final double                 eps   = FastMath.ulp(q.getField().getOne()).getReal();

        // Fourier series
        Complex sum1 = sc1.sin();
        Complex sum2 = sc1.cos();
        Complex sum3 = Complex.ZERO;
        Complex sum4 = Complex.ZERO;
        for (int n = 1; n < 16; ++n) {

            sgn   = -sgn;                            // (-1)ⁿ⁻¹     ← (-1)ⁿ
            qNN   = qNN.multiply(qTwoN).multiply(q); // q⁽ⁿ⁻¹⁾⁽ⁿ⁻¹⁾ ← qⁿⁿ
            qTwoN = qTwoN.multiply(qSquare);         // q²⁽ⁿ⁻¹⁾     ← q²ⁿ
            qNNp1 = qNNp1.multiply(qTwoN);           // q⁽ⁿ⁻¹⁾ⁿ     ← qⁿ⁽ⁿ⁺¹⁾

            sc2n1 = FieldSinCos.sum(sc2n1, sc1); // {sin|cos}([2n-1] z) ← {sin|cos}(2n z)
            sum3  = sum3.add(sc2n1.cos().multiply(qNN));
            sum4  = sum4.add(sc2n1.cos().multiply(qNN.multiply(sgn)));

            sc2n1 = FieldSinCos.sum(sc2n1, sc1); // {sin|cos}(2n z) ← {sin|cos}([2n+1] z)
            sum1  = sum1.add(sc2n1.sin().multiply(qNNp1.multiply(sgn)));
            sum2  = sum2.add(sc2n1.cos().multiply(qNNp1));

            if (qNNp1.getReal() <= eps) {
                // we have reach convergence
                return new Complex[] {sum1.multiply(qFourth.multiply(2)), sum2.multiply(qFourth.multiply(2)),
                                               sum3.multiply(2).add(1),            sum4.multiply(2).add(1)};
            }

        }

        // we were not able to compute the value
        throw new MathIllegalStateException(LocalizedCoreFormats.CONVERGENCE_FAILED);

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

        sn((plotter, z) -> plotter.jacobi.valuesN(z).sn()),
        cn((plotter, z) -> plotter.jacobi.valuesN(z).cn()),
        dn((plotter, z) -> plotter.jacobi.valuesN(z).dn()),
        cs((plotter, z) -> plotter.jacobi.valuesS(z).cs()),
        ds((plotter, z) -> plotter.jacobi.valuesS(z).ds()),
        ns((plotter, z) -> plotter.jacobi.valuesS(z).ns()),
        dc((plotter, z) -> plotter.jacobi.valuesC(z).dc()),
        nc((plotter, z) -> plotter.jacobi.valuesC(z).nc()),
        sc((plotter, z) -> plotter.jacobi.valuesC(z).sc()),
        sc2((plotter, z) -> {
            final Complex   zeta   = z.multiply(MathUtils.SEMI_PI).divide(plotter.bigK);
            final Complex[] theta0 = plotter.values(plotter.q, Complex.ZERO);
            final Complex[] thetaZ = plotter.values(plotter.q, zeta);

            final Complex t02 = theta0[1];
            final Complex t03 = theta0[2];
            final Complex t04 = theta0[3];
            final Complex tz1 = thetaZ[0];
            final Complex tz2 = thetaZ[1];
            final Complex tz3 = thetaZ[2];
            final Complex tz4 = thetaZ[3];
            Complex sn = t03.multiply(tz1)               .divide(t02.multiply(tz4));
            Complex cn = t04.multiply(tz2)               .divide(t02.multiply(tz4));
            Complex dn = t04.multiply(tz3)               .divide(t03.multiply(tz4));
            Complex sd = t03.multiply(t03).multiply(tz1) .divide(t02.multiply(t04).multiply(tz3));
            Complex cd = t03.multiply(tz2)               .divide(t02.multiply(tz3));
            Complex sc = t03.multiply(tz1)               .divide(t04.multiply(tz2));
            return sc;
        }),
        nd((plotter, z)    -> plotter.jacobi.valuesD(z).nd()),
        sd((plotter, z)    -> plotter.jacobi.valuesD(z).sd()),
        cd((plotter, z)    -> plotter.jacobi.valuesD(z).cd()),
        sin((plotter, z)   -> FastMath.sin(z)),
        cos((plotter, z)   -> FastMath.cos(z)),
        tan((plotter, z)   -> FastMath.tan(z)),
        asin((plotter, z)  -> FastMath.asin(z)),
        acos((plotter, z)  -> FastMath.acos(z)),
        atan((plotter, z)  -> FastMath.atan(z)),
        sinh((plotter, z)  -> FastMath.sinh(z)),
        cosh((plotter, z)  -> FastMath.cosh(z)),
        tanh((plotter, z)  -> FastMath.tanh(z)),
        asinh((plotter, z) -> FastMath.asinh(z)),
        acosh((plotter, z) -> FastMath.acosh(z)),
        atanh((plotter, z) -> FastMath.atanh(z));
        
        /** Function evaluator. */
        private final Evaluator evaluator;

        /** Simple constructor.
         * @param evaluator function evaluator
         */
        private Predefined(final Evaluator evaluator) {
            this.evaluator = evaluator;
        }


    }

}
