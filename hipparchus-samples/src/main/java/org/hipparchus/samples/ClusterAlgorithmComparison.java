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
package org.hipparchus.samples;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.geom.Ellipse2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import javax.swing.JComponent;
import javax.swing.JLabel;

import org.hipparchus.clustering.CentroidCluster;
import org.hipparchus.clustering.Cluster;
import org.hipparchus.clustering.Clusterable;
import org.hipparchus.clustering.Clusterer;
import org.hipparchus.clustering.DBSCANClusterer;
import org.hipparchus.clustering.DoublePoint;
import org.hipparchus.clustering.FuzzyKMeansClusterer;
import org.hipparchus.clustering.KMeansPlusPlusClusterer;
import org.hipparchus.geometry.euclidean.twod.Vector2D;
import org.hipparchus.random.RandomAdaptor;
import org.hipparchus.random.RandomDataGenerator;
import org.hipparchus.random.RandomGenerator;
import org.hipparchus.random.SobolSequenceGenerator;
import org.hipparchus.random.Well19937c;
import org.hipparchus.samples.ExampleUtils.ExampleFrame;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.Pair;
import org.hipparchus.util.SinCos;

/**
 * Plots clustering results for various algorithms and datasets.
 * Based on
 * <a href="http://scikit-learn.org/stable/auto_examples/cluster/plot_cluster_comparison.html">scikit learn</a>.
 */
public class ClusterAlgorithmComparison {

    /** Empty constructor.
     * <p>
     * This constructor is not strictly necessary, but it prevents spurious
     * javadoc warnings with JDK 18 and later.
     * </p>
     * @since 3.0
     */
    public ClusterAlgorithmComparison() { // NOPMD - unnecessary constructor added intentionally to make javadoc happy
        // nothing to do
    }

    /** Make circles patterns.
     * @param samples number of points
     * @param shuffle if true, shuffle points
     * @param noise noise to add to points position
     * @param factor reduction factor from outer to inner circle
     * @param random generator to use
     * @return circle patterns
     */
    public static List<Vector2D> makeCircles(int samples, boolean shuffle, double noise, double factor, final RandomGenerator random) {
        if (factor < 0 || factor > 1) {
            throw new IllegalArgumentException();
        }

        List<Vector2D> points = new ArrayList<Vector2D>();
        double range = 2.0 * FastMath.PI;
        double step = range / (samples / 2.0 + 1);
        for (double angle = 0; angle < range; angle += step) {
            Vector2D outerCircle = buildVector(angle);
            Vector2D innerCircle = outerCircle.scalarMultiply(factor);

            points.add(outerCircle.add(generateNoiseVector(random, noise)));
            points.add(innerCircle.add(generateNoiseVector(random, noise)));
        }

        if (shuffle) {
            Collections.shuffle(points, new RandomAdaptor(random));
        }

        return points;
    }

    /** Make Moons patterns.
     * @param samples number of points
     * @param shuffle if true, shuffle points
     * @param noise noise to add to points position
     * @param random generator to use
     * @return Moons patterns
     */
    public static List<Vector2D> makeMoons(int samples, boolean shuffle, double noise, RandomGenerator random) {

        int nSamplesOut = samples / 2;
        int nSamplesIn = samples - nSamplesOut;

        List<Vector2D> points = new ArrayList<Vector2D>();
        double range = FastMath.PI;
        double step = range / (nSamplesOut / 2.0);
        for (double angle = 0; angle < range; angle += step) {
            Vector2D outerCircle = buildVector(angle);
            points.add(outerCircle.add(generateNoiseVector(random, noise)));
        }

        step = range / (nSamplesIn / 2.0);
        for (double angle = 0; angle < range; angle += step) {
            final SinCos sc = FastMath.sinCos(angle);
            Vector2D innerCircle = new Vector2D(1 - sc.cos(), 1 - sc.sin() - 0.5);
            points.add(innerCircle.add(generateNoiseVector(random, noise)));
        }

        if (shuffle) {
            Collections.shuffle(points, new RandomAdaptor(random));
        }

        return points;
    }

    /** Make blobs patterns.
     * @param samples number of points
     * @param centers number of centers
     * @param clusterStd standard deviation of cluster
     * @param min range min value
     * @param max range max value
     * @param shuffle if true, shuffle points
     * @param random generator to use
     * @return blobs patterns
     */
    public static List<Vector2D> makeBlobs(int samples, int centers, double clusterStd,
                                           double min, double max, boolean shuffle, RandomGenerator random) {

        final RandomDataGenerator randomDataGenerator = RandomDataGenerator.of(random);
        //NormalDistribution dist = new NormalDistribution(random, 0.0, clusterStd);

        double range = max - min;
        Vector2D[] centerPoints = new Vector2D[centers];
        for (int i = 0; i < centers; i++) {
            double x = random.nextDouble() * range + min;
            double y = random.nextDouble() * range + min;
            centerPoints[i] = new Vector2D(x, y);
        }

        int[] nSamplesPerCenter = new int[centers];
        int count = samples / centers;
        Arrays.fill(nSamplesPerCenter, count);

        for (int i = 0; i < samples % centers; i++) {
            nSamplesPerCenter[i]++;
        }

        List<Vector2D> points = new ArrayList<Vector2D>();
        for (int i = 0; i < centers; i++) {
            for (int j = 0; j < nSamplesPerCenter[i]; j++) {
                Vector2D point = new Vector2D(randomDataGenerator.nextNormal(0, clusterStd),
                                              randomDataGenerator.nextNormal(0, clusterStd));
                points.add(point.add(centerPoints[i]));
            }
        }

        if (shuffle) {
            Collections.shuffle(points, new RandomAdaptor(random));
        }

        return points;
    }

    /** Make Sobol patterns.
     * @param samples number of points
     * @return Moons patterns
     */
    public static List<Vector2D> makeSobol(int samples) {
        SobolSequenceGenerator generator = new SobolSequenceGenerator(2);
        generator.skipTo(999999);
        List<Vector2D> points = new ArrayList<Vector2D>();
        for (double i = 0; i < samples; i++) {
            double[] vector = generator.nextVector();
            vector[0] = vector[0] * 2 - 1;
            vector[1] = vector[1] * 2 - 1;
            Vector2D point = new Vector2D(vector);
            points.add(point);
        }

        return points;
    }

    /** Generate a random vector.
     * @param randomGenerator random generator to use
     * @param noise noise level
     * @return random vector
     */
    public static Vector2D generateNoiseVector(RandomGenerator randomGenerator, double noise) {
        final RandomDataGenerator randomDataGenerator = RandomDataGenerator.of(randomGenerator);
        return new Vector2D(randomDataGenerator.nextNormal(0, noise), randomDataGenerator.nextNormal(0, noise));
    }

    /** Normolize points in a rectangular area
     * @param input input points
     * @param minX range min value in X
     * @param maxX range max value in X
     * @param minY range min value in Y
     * @param maxY range max value in Y
     * @return normalized points
     */
    public static List<DoublePoint> normalize(final List<Vector2D> input, double minX, double maxX, double minY, double maxY) {
        double rangeX = maxX - minX;
        double rangeY = maxY - minY;
        List<DoublePoint> points = new ArrayList<DoublePoint>();
        for (Vector2D p : input) {
            double[] arr = p.toArray();
            arr[0] = (arr[0] - minX) / rangeX * 2 - 1;
            arr[1] = (arr[1] - minY) / rangeY * 2 - 1;
            points.add(new DoublePoint(arr));
        }
        return points;
    }

    /**
     * Build the 2D vector corresponding to the given angle.
     * @param alpha angle
     * @return the corresponding 2D vector
     */
    private static Vector2D buildVector(final double alpha) {
        final SinCos sc = FastMath.sinCos(alpha);
        return new Vector2D(sc.cos(), sc.sin());
    }

    /** Display frame. */
    @SuppressWarnings("serial")
    public static class Display extends ExampleFrame {

        /** Simple consructor. */
        public Display() {
            setTitle("Hipparchus: Cluster algorithm comparison");
            setSize(800, 800);

            setLayout(new GridBagLayout());

            int nSamples = 1500;

            RandomGenerator rng = new Well19937c(0);
            List<List<DoublePoint>> datasets = new ArrayList<List<DoublePoint>>();

            datasets.add(normalize(makeCircles(nSamples, true, 0.04, 0.5, rng), -1, 1, -1, 1));
            datasets.add(normalize(makeMoons(nSamples, true, 0.04, rng), -1, 2, -1, 1));
            datasets.add(normalize(makeBlobs(nSamples, 3, 1.0, -10, 10, true, rng), -12, 12, -12, 12));
            datasets.add(normalize(makeSobol(nSamples), -1, 1, -1, 1));

            List<Pair<String, Clusterer<DoublePoint>>> algorithms = new ArrayList<Pair<String, Clusterer<DoublePoint>>>();

            algorithms.add(new Pair<String, Clusterer<DoublePoint>>("KMeans\n(k=2)", new KMeansPlusPlusClusterer<DoublePoint>(2)));
            algorithms.add(new Pair<String, Clusterer<DoublePoint>>("KMeans\n(k=3)", new KMeansPlusPlusClusterer<DoublePoint>(3)));
            algorithms.add(new Pair<String, Clusterer<DoublePoint>>("FuzzyKMeans\n(k=3, fuzzy=2)", new FuzzyKMeansClusterer<DoublePoint>(3, 2)));
            algorithms.add(new Pair<String, Clusterer<DoublePoint>>("FuzzyKMeans\n(k=3, fuzzy=10)", new FuzzyKMeansClusterer<DoublePoint>(3, 10)));
            algorithms.add(new Pair<String, Clusterer<DoublePoint>>("DBSCAN\n(eps=.1, min=3)", new DBSCANClusterer<DoublePoint>(0.1, 3)));

            GridBagConstraints c = new GridBagConstraints();
            c.fill = GridBagConstraints.VERTICAL;
            c.gridx = 0;
            c.gridy = 0;
            c.insets = new Insets(2, 2, 2, 2);

            for (Pair<String, Clusterer<DoublePoint>> pair : algorithms) {
                JLabel text = new JLabel("<html><body>" + pair.getFirst().replace("\n", "<br>"));
                add(text, c);
                c.gridx++;
            }
            c.gridy++;

            for (List<DoublePoint> dataset : datasets) {
                c.gridx = 0;
                for (Pair<String, Clusterer<DoublePoint>> pair : algorithms) {
                    long start = System.currentTimeMillis();
                    List<? extends Cluster<DoublePoint>> clusters = pair.getSecond().cluster(dataset);
                    long end = System.currentTimeMillis();
                    add(new ClusterPlot(clusters, end - start), c);
                    c.gridx++;
                }
                c.gridy++;
            }
        }

    }

    /** Plot component. */
    @SuppressWarnings("serial")
    public static class ClusterPlot extends JComponent {

        /** Padding. */
        private static double PAD = 10;

        /** Clusters. */
        private List<? extends Cluster<DoublePoint>> clusters;

        /** Duration of the computation. */
        private long duration;

        /** Simple constructor.
         * @param clusters clusters to plot
         * @param duration duration of the computation
         */
        public ClusterPlot(final List<? extends Cluster<DoublePoint>> clusters, long duration) {
            this.clusters = clusters;
            this.duration = duration;
        }

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2 = (Graphics2D)g;
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                                RenderingHints.VALUE_ANTIALIAS_ON);

            int w = getWidth();
            int h = getHeight();

            g2.clearRect(0, 0, w, h);

            g2.setPaint(Color.black);
            g2.drawRect(0, 0, w - 1, h - 1);

            int index = 0;
            Color[] colors = new Color[] { Color.red, Color.blue, Color.green.darker() };
            for (Cluster<DoublePoint> cluster : clusters) {
                g2.setPaint(colors[index++]);
                for (DoublePoint point : cluster.getPoints()) {
                    Clusterable p = transform(point, w, h);
                    double[] arr = p.getPoint();
                    g2.fill(new Ellipse2D.Double(arr[0] - 1, arr[1] - 1, 3, 3));
                }

                if (cluster instanceof CentroidCluster) {
                    Clusterable p = transform(((CentroidCluster<?>) cluster).getCenter(), w, h);
                    double[] arr = p.getPoint();
                    Shape s = new Ellipse2D.Double(arr[0] - 4, arr[1] - 4, 8, 8);
                    g2.fill(s);
                    g2.setPaint(Color.black);
                    g2.draw(s);
                }
            }

            g2.setPaint(Color.black);
            g2.drawString(String.format("%.2f s", duration / 1e3), w - 40, h - 5);
        }

        @Override
        public Dimension getPreferredSize() {
            return new Dimension(150, 150);
        }

        private Clusterable transform(Clusterable point, int width, int height) {
            double[] arr = point.getPoint();
            return new DoublePoint(new double[] { PAD + (arr[0] + 1) / 2.0 * (width - 2 * PAD),
                                                  height - PAD - (arr[1] + 1) / 2.0 * (height - 2 * PAD) });
        }
    }

    /** Example entry point.
     * @param args arguments (not used)
     */
    public static void main(String[] args) {
        ExampleUtils.showExampleFrame(new Display());
    }

}
