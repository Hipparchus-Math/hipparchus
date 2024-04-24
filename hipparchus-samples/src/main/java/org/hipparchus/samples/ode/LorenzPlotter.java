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

package org.hipparchus.samples.ode;

import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.ode.sampling.StepNormalizer;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.Locale;

/** Program plotting the Lorenz attractor.
 */
public class LorenzPlotter {

    /** Duration. */
    private double duration;

    /** Step. */
    private double step;

    /** Sigma. */
    private double sigma;

    /** Rho. */
    private double rho;

    /** Beta. */
    private double beta;

    /** Output directory (display to terminal if null). */
    private File output;

    /** Output width. */
    private int width;

    /** Output height. */
    private int height;

    /** View X rotation. */
    private double viewXRot;

    /** View Z rotation. */
    private double viewZRot;

    /** Default constructor.
     */
    private LorenzPlotter() {
        duration   = 125.0;
        step       = 0.002;
        sigma      = 10;
        rho        = 28;
        beta       = 8.0 / 3.0;
        output     = null;
        width      = 1000;
        height     = 1000;
        viewXRot   = 70;
        viewZRot   = 20;
    }

    /** Main program.
     * @param args program arguments
     */
    public static void main(String[] args) {
        final LorenzPlotter plotter = new LorenzPlotter();
        try {
            for (int i = 0; i < args.length; ++i) {
                switch (args[i]) {
                    case "--help" :
                        usage(0);
                        break;
                    case "--duration" :
                        plotter.duration = Double.parseDouble(args[++i]);
                        break;
                    case "--step" :
                        plotter.step = Double.parseDouble(args[++i]);
                        break;
                    case "--sigma" :
                        plotter.sigma = Double.parseDouble(args[++i]);
                        break;
                    case "--rhon" :
                        plotter.rho = Double.parseDouble(args[++i]);
                        break;
                    case "--beta" :
                        plotter.beta = Double.parseDouble(args[++i]);
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
                    case "--view" :
                        plotter.viewXRot = Double.parseDouble(args[++i]);
                        plotter.viewZRot = Double.parseDouble(args[++i]);
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
        System.err.println("usage: java org.hipparchus.samples.ode.LorenzPlotter" +
                           " [--help]" +
                           " [--duration duration] [--step step]" +
                           " [--sigma sigma] [--rho rho] [--beta beta]" +
                           " [--output-dir directory]" +
                           " [--view xRot zRot]");
        System.exit(status);

    }

    /** Plot the system.
     * @throws IOException if gnuplot process cannot be run
     */
    public void plot() throws IOException {
            final ProcessBuilder pb = new ProcessBuilder("gnuplot").
                            redirectOutput(ProcessBuilder.Redirect.INHERIT).
                            redirectError(ProcessBuilder.Redirect.INHERIT);
            pb.environment().remove("XDG_SESSION_TYPE");
            final Process gnuplot = pb.start();
            try (PrintStream out = new PrintStream(gnuplot.getOutputStream(), false, StandardCharsets.UTF_8.name())) {
                final File outputFile;
                if (output == null) {
                    out.format(Locale.US, "set terminal qt size %d, %d title 'Lorenz plotter'%n", width, height);
                    outputFile = null;
                } else {
                    out.format(Locale.US, "set terminal pngcairo size %d, %d%n", width, height);
                    outputFile = new File(output, "Lorenz-attractor.png");
                    out.format(Locale.US, "set output '%s'%n", outputFile.getAbsolutePath());
                }
                out.format(Locale.US, "set xlabel 'X'%n");
                out.format(Locale.US, "set ylabel 'Y'%n");
                out.format(Locale.US, "set zlabel 'Z'%n");
                out.format(Locale.US, "set key off%n");
                out.format(Locale.US, "unset colorbox%n");
                out.format(Locale.US, "set palette model RGB defined (0 0.45 0.66 0.85,0.25 0.725 0.9 1, 0.33 0.85 0.95 1, 0.33 0.68 0.82 0.65,0.4 0.58 0.75 0.53,0.575 0.95 0.9 0.75, 0.8 0.66 0.52 0.32,0.85 0.66 0.6 0.5, 1 0.95 0.95 0.95)%n");
                out.format(Locale.US, "$data <<EOD%n");

                ODEIntegrator integrator = new DormandPrince853Integrator(1.0e-3, 10.0, 1.0e-12, 1.0e-12);
                integrator.addStepHandler(new StepNormalizer(step, (state, isLast) ->
                    out.format(Locale.US, "%.6f %.3f %.3f %.3f%n",
                               state.getTime() / duration,
                               state.getCompleteState()[0],
                               state.getCompleteState()[1],
                               state.getCompleteState()[2])));
                integrator.integrate(new LorenzOde(),
                                     new ODEState(0, new double[] { -8, 8, rho - 1 }),
                                     duration);
                out.format(Locale.US, "EOD%n");
                out.format(Locale.US, "set view %f, %f%n", viewXRot, viewZRot);
                out.format(Locale.US, "splot $data using 2:3:4:1 with lines lc palette%n");
                if (output == null) {
                    out.format(Locale.US, "pause mouse close%n");
                } else {
                    System.out.format(Locale.US, "output written to %s%n",
                                      outputFile.getAbsolutePath());
                }
            }
    }

    /** Lorenz system. */
    private class LorenzOde implements OrdinaryDifferentialEquation {

        /** {@inheritDoc} */
        @Override
        public int getDimension() {
            return 3;
        }

        /** {@inheritDoc} */
        @Override
        public double[] computeDerivatives(final double t, final double[] y) {
            return new double[] {
                sigma * (y[1] - y[0]),
                y[0] * (rho - y[2]) - y[1],
                y[0] * y[1] - beta * y[2]
            };
        }
    }

}
