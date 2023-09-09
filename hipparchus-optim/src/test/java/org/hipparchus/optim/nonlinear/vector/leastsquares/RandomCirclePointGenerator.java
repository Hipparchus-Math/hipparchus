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
package org.hipparchus.optim.nonlinear.vector.leastsquares;

import org.hipparchus.geometry.euclidean.twod.Vector2D;
import org.hipparchus.random.RandomDataGenerator;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;

/**
 * Factory for generating a cloud of points that approximate a circle.
 */
public class RandomCirclePointGenerator {
    /** Radius of the circle. */
    private final double radius;
    /** Random data generator */
    private final RandomDataGenerator randomDataGenerator;
    /** Error on the x-coordinate of the circumference points. */
    private final double xSigma;
    /** Error on the y-coordinate of the circumference points. */
    private final double ySigma;
    /** Abscissa of the circle center. */
    private final double x;
    /** Ordinate of the circle center. */
    private final double y;


    /**
     * @param x Abscissa of the circle center.
     * @param y Ordinate of the circle center.
     * @param radius Radius of the circle.
     * @param xSigma Error on the x-coordinate of the circumference points.
     * @param ySigma Error on the y-coordinate of the circumference points.
     * @param seed RNG seed.
     */
    public RandomCirclePointGenerator(double x,
                                      double y,
                                      double radius,
                                      double xSigma,
                                      double ySigma,
                                      long seed) {
        randomDataGenerator = new RandomDataGenerator(seed);
        this.radius = radius;
        this.xSigma = xSigma;
        this.ySigma = ySigma;
        this.x = x;
        this.y = y;
    }

    /**
     * Point generator.
     *
     * @param n Number of points to create.
     * @return the cloud of {@code n} points.
     */
    public Vector2D[] generate(int n) {
        final Vector2D[] cloud = new Vector2D[n];
        for (int i = 0; i < n; i++) {
            cloud[i] = create();
        }
        return cloud;
    }

    /**
     * Create one point.
     *
     * @return a point.
     */
    private Vector2D create() {
        final double t = randomDataGenerator.nextUniform(0,  MathUtils.TWO_PI);
        final double pX = randomDataGenerator.nextNormal(x, xSigma) + radius * FastMath.cos(t);
        final double pY = randomDataGenerator.nextNormal(y, ySigma) + radius * FastMath.sin(t);

        return new Vector2D(pX, pY);
    }
}
