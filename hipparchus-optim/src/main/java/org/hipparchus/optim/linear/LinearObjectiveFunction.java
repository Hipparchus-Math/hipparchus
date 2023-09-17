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
package org.hipparchus.optim.linear;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;

import org.hipparchus.analysis.MultivariateFunction;
import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.RealVector;
import org.hipparchus.optim.OptimizationData;

/**
 * An objective function for a linear optimization problem.
 * <p>
 * A linear objective function has one the form:
 * \[
 * c_1 x_1 + \ldots c_n x_n + d
 * \]
 * The c<sub>i</sub> and d are the coefficients of the equation,
 * the x<sub>i</sub> are the coordinates of the current point.
 * </p>
 *
 */
public class LinearObjectiveFunction
    implements MultivariateFunction,
               OptimizationData,
               Serializable {
    /** Serializable version identifier. */
    private static final long serialVersionUID = -4531815507568396090L;
    /** Coefficients of the linear equation (c<sub>i</sub>). */
    private final transient RealVector coefficients;
    /** Constant term of the linear equation. */
    private final double constantTerm;

    /** Simple constructor.
     * @param coefficients Coefficients for the linear equation being optimized.
     * @param constantTerm Constant term of the linear equation.
     */
    public LinearObjectiveFunction(double[] coefficients, double constantTerm) {
        this(new ArrayRealVector(coefficients), constantTerm);
    }

    /** Simple constructor.
     * @param coefficients Coefficients for the linear equation being optimized.
     * @param constantTerm Constant term of the linear equation.
     */
    public LinearObjectiveFunction(RealVector coefficients, double constantTerm) {
        this.coefficients = coefficients;
        this.constantTerm = constantTerm;
    }

    /**
     * Gets the coefficients of the linear equation being optimized.
     *
     * @return coefficients of the linear equation being optimized.
     */
    public RealVector getCoefficients() {
        return coefficients;
    }

    /**
     * Gets the constant of the linear equation being optimized.
     *
     * @return constant of the linear equation being optimized.
     */
    public double getConstantTerm() {
        return constantTerm;
    }

    /**
     * Computes the value of the linear equation at the current point.
     *
     * @param point Point at which linear equation must be evaluated.
     * @return the value of the linear equation at the current point.
     */
    @Override
    public double value(final double[] point) {
        return value(new ArrayRealVector(point, false));
    }

    /**
     * Computes the value of the linear equation at the current point.
     *
     * @param point Point at which linear equation must be evaluated.
     * @return the value of the linear equation at the current point.
     */
    public double value(final RealVector point) {
        return coefficients.dotProduct(point) + constantTerm;
    }

    /** {@inheritDoc} */
    @Override
    public boolean equals(Object other) {
        if (this == other) {
            return true;
        }
        if (other instanceof LinearObjectiveFunction) {
            LinearObjectiveFunction rhs = (LinearObjectiveFunction) other;
          return (constantTerm == rhs.constantTerm) && coefficients.equals(rhs.coefficients);
        }

        return false;
    }

    /** {@inheritDoc} */
    @Override
    public int hashCode() {
        return Double.valueOf(constantTerm).hashCode() ^ coefficients.hashCode();
    }

    /**
     * Serialize the instance.
     * @param oos stream where object should be written
     * @throws IOException if object cannot be written to stream
     */
    private void writeObject(ObjectOutputStream oos)
        throws IOException {
        oos.defaultWriteObject();
        final int n = coefficients.getDimension();
        oos.writeInt(n);
        for (int i = 0; i < n; ++i) {
            oos.writeDouble(coefficients.getEntry(i));
        }
    }

    /**
     * Deserialize the instance.
     * @param ois stream from which the object should be read
     * @throws ClassNotFoundException if a class in the stream cannot be found
     * @throws IOException if object cannot be read from the stream
     */
    private void readObject(ObjectInputStream ois)
      throws ClassNotFoundException, IOException {
        ois.defaultReadObject();

        // read the vector data
        final int n = ois.readInt();
        final double[] data = new double[n];
        for (int i = 0; i < n; ++i) {
            data[i] = ois.readDouble();
        }

        try {
            // create the instance
            ArrayRealVector vector = new ArrayRealVector(data, false);
            final java.lang.reflect.Field f = getClass().getDeclaredField("coefficients");
            f.setAccessible(true); // NOPMD
            f.set(this, vector);
        } catch (NoSuchFieldException | IllegalArgumentException | IllegalAccessException e) {
            IOException ioe = new IOException();
            ioe.initCause(e);
            throw ioe;
        }
    }
}
