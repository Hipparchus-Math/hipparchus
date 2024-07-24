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
package org.hipparchus.clustering;

import org.hamcrest.CoreMatchers;
import org.hipparchus.clustering.distance.CanberraDistance;
import org.hipparchus.clustering.distance.DistanceMeasure;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.exception.NullArgumentException;
import org.hipparchus.random.JDKRandomGenerator;
import org.hipparchus.random.RandomGenerator;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.junit.jupiter.api.Assertions.assertThrows;

/**
 * Test cases for FuzzyKMeansClusterer.
 *
 */
public class FuzzyKMeansClustererTest {

    @Test
    public void testCluster() {
        final List<DoublePoint> points = new ArrayList<DoublePoint>();

        // create 10 data points: [1], ... [10]
        for (int i = 1; i <= 10; i++) {
            final DoublePoint p = new DoublePoint(new double[] { i } );
            points.add(p);
        }

        final FuzzyKMeansClusterer<DoublePoint> transformer =
                new FuzzyKMeansClusterer<DoublePoint>(3, 2.0);
        final List<CentroidCluster<DoublePoint>> clusters = transformer.cluster(points);

        // we expect 3 clusters:
        //   [1], [2], [3]
        //   [4], [5], [6], [7]
        //   [8], [9], [10]
        final List<DoublePoint> clusterOne = Arrays.asList(points.get(0), points.get(1), points.get(2));
        final List<DoublePoint> clusterTwo = Arrays.asList(points.get(3), points.get(4), points.get(5), points.get(6));
        final List<DoublePoint> clusterThree = Arrays.asList(points.get(7), points.get(8), points.get(9));

        boolean cluster1Found = false;
        boolean cluster2Found = false;
        boolean cluster3Found = false;
        Assertions.assertEquals(3, clusters.size());
        for (final Cluster<DoublePoint> cluster : clusters) {
            if (cluster.getPoints().containsAll(clusterOne)) {
                cluster1Found = true;
            }
            if (cluster.getPoints().containsAll(clusterTwo)) {
                cluster2Found = true;
            }
            if (cluster.getPoints().containsAll(clusterThree)) {
                cluster3Found = true;
            }
        }
        Assertions.assertTrue(cluster1Found);
        Assertions.assertTrue(cluster2Found);
        Assertions.assertTrue(cluster3Found);
    }

    @Test
    public void testTooSmallFuzzynessFactor() {
        assertThrows(MathIllegalArgumentException.class, () -> {
            new FuzzyKMeansClusterer<DoublePoint>(3, 1.0);
        });
    }

    @Test
    public void testNullDataset() {
        assertThrows(NullArgumentException.class, () -> {
            final FuzzyKMeansClusterer<DoublePoint> clusterer = new FuzzyKMeansClusterer<DoublePoint>(3, 2.0);
            clusterer.cluster(null);
        });
    }

    @Test
    public void testGetters() {
        final DistanceMeasure measure = new CanberraDistance();
        final RandomGenerator random = new JDKRandomGenerator();
        final FuzzyKMeansClusterer<DoublePoint> clusterer =
                new FuzzyKMeansClusterer<DoublePoint>(3, 2.0, 100, measure, 1e-6, random);

        Assertions.assertEquals(3, clusterer.getK());
        Assertions.assertEquals(2.0, clusterer.getFuzziness(), 1e-6);
        Assertions.assertEquals(100, clusterer.getMaxIterations());
        Assertions.assertEquals(1e-6, clusterer.getEpsilon(), 1e-12);
        assertThat(clusterer.getDistanceMeasure(), CoreMatchers.is(measure));
        assertThat(clusterer.getRandomGenerator(), CoreMatchers.is(random));
    }

    @Test
    public void testSingleCluster() {
        final List<DoublePoint> points = new ArrayList<DoublePoint>();
        points.add(new DoublePoint(new double[] { 1, 1 }));

        final FuzzyKMeansClusterer<DoublePoint> transformer =
                new FuzzyKMeansClusterer<DoublePoint>(1, 2.0);
        final List<CentroidCluster<DoublePoint>> clusters = transformer.cluster(points);

        Assertions.assertEquals(1, clusters.size());
    }

    @Test
    public void testClusterCenterEqualsPoints() {
        final List<DoublePoint> points = new ArrayList<DoublePoint>();
        points.add(new DoublePoint(new double[] { 1, 1 }));
        points.add(new DoublePoint(new double[] { 1.00001, 1.00001 }));
        points.add(new DoublePoint(new double[] { 2, 2 }));
        points.add(new DoublePoint(new double[] { 3, 3 }));

        final FuzzyKMeansClusterer<DoublePoint> transformer =
                new FuzzyKMeansClusterer<DoublePoint>(3, 2.0);
        final List<CentroidCluster<DoublePoint>> clusters = transformer.cluster(points);

        Assertions.assertEquals(3, clusters.size());
    }

}
