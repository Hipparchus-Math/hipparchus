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

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.List;

public class MultiKMeansPlusPlusClustererTest {

    @Test
    public void dimension2() {
        MultiKMeansPlusPlusClusterer<DoublePoint> transformer =
            new MultiKMeansPlusPlusClusterer<DoublePoint>(
                    new KMeansPlusPlusClusterer<DoublePoint>(3, 10), 5);

        DoublePoint[] points = new DoublePoint[] {

                // first expected cluster
                new DoublePoint(new int[] { -15,  3 }),
                new DoublePoint(new int[] { -15,  4 }),
                new DoublePoint(new int[] { -15,  5 }),
                new DoublePoint(new int[] { -14,  3 }),
                new DoublePoint(new int[] { -14,  5 }),
                new DoublePoint(new int[] { -13,  3 }),
                new DoublePoint(new int[] { -13,  4 }),
                new DoublePoint(new int[] { -13,  5 }),

                // second expected cluster
                new DoublePoint(new int[] { -1,  0 }),
                new DoublePoint(new int[] { -1, -1 }),
                new DoublePoint(new int[] {  0, -1 }),
                new DoublePoint(new int[] {  1, -1 }),
                new DoublePoint(new int[] {  1, -2 }),

                // third expected cluster
                new DoublePoint(new int[] { 13,  3 }),
                new DoublePoint(new int[] { 13,  4 }),
                new DoublePoint(new int[] { 14,  4 }),
                new DoublePoint(new int[] { 14,  7 }),
                new DoublePoint(new int[] { 16,  5 }),
                new DoublePoint(new int[] { 16,  6 }),
                new DoublePoint(new int[] { 17,  4 }),
                new DoublePoint(new int[] { 17,  7 })

        };
        List<CentroidCluster<DoublePoint>> clusters = transformer.cluster(Arrays.asList(points));

        Assertions.assertEquals(3, clusters.size());
        boolean cluster1Found = false;
        boolean cluster2Found = false;
        boolean cluster3Found = false;
        double epsilon = 1e-6;
        for (CentroidCluster<DoublePoint> cluster : clusters) {
            Clusterable center = cluster.getCenter();
            double[] point = center.getPoint();
            if (point[0] < 0) {
                cluster1Found = true;
                Assertions.assertEquals(8, cluster.getPoints().size());
                Assertions.assertEquals(-14, point[0], epsilon);
                Assertions.assertEquals( 4, point[1], epsilon);
            } else if (point[1] < 0) {
                cluster2Found = true;
                Assertions.assertEquals(5, cluster.getPoints().size());
                Assertions.assertEquals( 0, point[0], epsilon);
                Assertions.assertEquals(-1, point[1], epsilon);
            } else {
                cluster3Found = true;
                Assertions.assertEquals(8, cluster.getPoints().size());
                Assertions.assertEquals(15, point[0], epsilon);
                Assertions.assertEquals(5, point[1], epsilon);
            }
        }
        Assertions.assertTrue(cluster1Found);
        Assertions.assertTrue(cluster2Found);
        Assertions.assertTrue(cluster3Found);

    }

}
