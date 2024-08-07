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

package org.hipparchus.clustering.evaluation;

import org.hipparchus.clustering.Cluster;
import org.hipparchus.clustering.DoublePoint;
import org.hipparchus.clustering.distance.EuclideanDistance;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

class SumOfClusterVariancesTest {

    private ClusterEvaluator<DoublePoint> evaluator;

    @BeforeEach
    void setUp() {
        evaluator = new SumOfClusterVariances<DoublePoint>(new EuclideanDistance());
    }

    @Test
    void testScore() {
        final DoublePoint[] points1 = new DoublePoint[] {
                new DoublePoint(new double[] { 1 }),
                new DoublePoint(new double[] { 2 }),
                new DoublePoint(new double[] { 3 })
        };

        final DoublePoint[] points2 = new DoublePoint[] {
                new DoublePoint(new double[] { 1 }),
                new DoublePoint(new double[] { 5 }),
                new DoublePoint(new double[] { 10 })
        };

        final List<Cluster<DoublePoint>> clusters = new ArrayList<Cluster<DoublePoint>>();

        final Cluster<DoublePoint> cluster1 = new Cluster<DoublePoint>();
        for (DoublePoint p : points1) {
            cluster1.addPoint(p);
        }
        clusters.add(cluster1);

        assertEquals(1.0/3.0, evaluator.score(clusters), 1e-6);

        final Cluster<DoublePoint> cluster2 = new Cluster<DoublePoint>();
        for (DoublePoint p : points2) {
            cluster2.addPoint(p);
        }
        clusters.add(cluster2);

        assertEquals(6.148148148, evaluator.score(clusters), 1e-6);
    }

    @Test
    void testOrdering() {
        assertTrue(evaluator.isBetterScore(10, 20));
        assertFalse(evaluator.isBetterScore(20, 1));
    }
}
