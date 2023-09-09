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

package org.hipparchus.ode.nonstiff;

import org.hipparchus.CalculusFieldElement;
import org.hipparchus.exception.LocalizedCoreFormats;
import org.hipparchus.exception.MathIllegalArgumentException;
import org.hipparchus.ode.LocalizedODEFormats;
import org.hipparchus.util.FastMath;

/** Helper for adaptive stepsize control.
 * @since 2.0
 */

public class StepsizeHelper {

    /** Allowed absolute scalar error. */
    private double scalAbsoluteTolerance;

    /** Allowed relative scalar error. */
    private double scalRelativeTolerance;

    /** Allowed absolute vectorial error. */
    private double[] vecAbsoluteTolerance;

    /** Allowed relative vectorial error. */
    private double[] vecRelativeTolerance;

    /** Main set dimension. */
    private int mainSetDimension;

    /** User supplied initial step. */
    private double initialStep;

    /** Minimal step. */
    private double minStep;

    /** Maximal step. */
    private double maxStep;

    /** Simple constructor.
     * @param minStep minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param maxStep maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param scalAbsoluteTolerance allowed absolute error
     * @param scalRelativeTolerance allowed relative error
     */
    public StepsizeHelper(final double minStep, final double maxStep,
                          final double scalAbsoluteTolerance,
                          final double scalRelativeTolerance) {
        this.minStep     = FastMath.abs(minStep);
        this.maxStep     = FastMath.abs(maxStep);
        this.initialStep = -1;

        this.scalAbsoluteTolerance = scalAbsoluteTolerance;
        this.scalRelativeTolerance = scalRelativeTolerance;
        this.vecAbsoluteTolerance  = null;
        this.vecRelativeTolerance  = null;
    }

    /** Simple constructor..
     * @param minStep minimal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param maxStep maximal step (sign is irrelevant, regardless of
     * integration direction, forward or backward), the last step can
     * be smaller than this
     * @param vecAbsoluteTolerance allowed absolute error
     * @param vecRelativeTolerance allowed relative error
     */
    public StepsizeHelper(final double minStep, final double maxStep,
                          final double[] vecAbsoluteTolerance,
                          final double[] vecRelativeTolerance) {

        this.minStep     = FastMath.abs(minStep);
        this.maxStep     = FastMath.abs(maxStep);
        this.initialStep = -1;

       this.scalAbsoluteTolerance = 0;
       this.scalRelativeTolerance = 0;
       this.vecAbsoluteTolerance  = vecAbsoluteTolerance.clone();
       this.vecRelativeTolerance  = vecRelativeTolerance.clone();

    }

    /** Set main set dimension.
     * @param mainSetDimension dimension of the main set
     * @exception MathIllegalArgumentException if adaptive step size integrators
     * tolerance arrays dimensions are not compatible with equations settings
     */
    protected void setMainSetDimension(final int mainSetDimension) throws MathIllegalArgumentException {
        this.mainSetDimension = mainSetDimension;

        if (vecAbsoluteTolerance != null && vecAbsoluteTolerance.length != mainSetDimension) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   mainSetDimension, vecAbsoluteTolerance.length);
        }

        if (vecRelativeTolerance != null && vecRelativeTolerance.length != mainSetDimension) {
            throw new MathIllegalArgumentException(LocalizedCoreFormats.DIMENSIONS_MISMATCH,
                                                   mainSetDimension, vecRelativeTolerance.length);
        }
    }

    /** Get the main set dimension.
     * @return main set dimension
     */
    public int getMainSetDimension() {
        return mainSetDimension;
    }

    /** Get the relative tolerance for one component.
     * @param i component to select
     * @return relative tolerance for selected component
     */
    public double getRelativeTolerance(final int i) {
        return vecAbsoluteTolerance == null ? scalRelativeTolerance : vecRelativeTolerance[i];
    }

    /** Get the tolerance for one component.
     * @param i component to select
     * @param scale scale factor for relative tolerance (i.e. y[i])
     * @return tolerance for selected component
     */
    public double getTolerance(final int i, final double scale) {
        return vecAbsoluteTolerance == null ?
               scalAbsoluteTolerance   + scalRelativeTolerance   * scale :
               vecAbsoluteTolerance[i] + vecRelativeTolerance[i] * scale;
    }

    /** Get the tolerance for one component.
     * @param i component to select
     * @param scale scale factor for relative tolerance (i.e. y[i])
     * @param <T> type of the field elements
     * @return tolerance for selected component
     */
    public <T extends CalculusFieldElement<T>> T getTolerance(final int i, final T scale) {
        return vecAbsoluteTolerance == null ?
               scale.multiply(scalRelativeTolerance).add(scalAbsoluteTolerance) :
               scale.multiply(vecRelativeTolerance[i]).add(vecAbsoluteTolerance[i]);
    }

    /** Filter the integration step.
     * @param h signed step
     * @param forward forward integration indicator
     * @param acceptSmall if true, steps smaller than the minimal value
     * are silently increased up to this value, if false such small
     * steps generate an exception
     * @return a bounded integration step (h if no bound is reach, or a bounded value)
     * @exception MathIllegalArgumentException if the step is too small and acceptSmall is false
     */
    public double filterStep(final double h, final boolean forward, final boolean acceptSmall)
        throws MathIllegalArgumentException {

        double filteredH = h;
        if (FastMath.abs(h) < minStep) {
            if (acceptSmall) {
                filteredH = forward ? minStep : -minStep;
            } else {
                throw new MathIllegalArgumentException(LocalizedODEFormats.MINIMAL_STEPSIZE_REACHED_DURING_INTEGRATION,
                                                       FastMath.abs(h), minStep, true);
            }
        }

        if (filteredH > maxStep) {
            filteredH = maxStep;
        } else if (filteredH < -maxStep) {
            filteredH = -maxStep;
        }

        return filteredH;

    }

    /** Filter the integration step.
     * @param h signed step
     * @param forward forward integration indicator
     * @param acceptSmall if true, steps smaller than the minimal value
     * are silently increased up to this value, if false such small
     * steps generate an exception
     * @param <T> type of the field elements
     * @return a bounded integration step (h if no bound is reach, or a bounded value)
     * @exception MathIllegalArgumentException if the step is too small and acceptSmall is false
     */
    public <T extends CalculusFieldElement<T>> T filterStep(final T h, final boolean forward, final boolean acceptSmall)
        throws MathIllegalArgumentException {

        T filteredH = h;
        if (h.abs().subtract(minStep).getReal() < 0) {
            if (acceptSmall) {
                filteredH = h.getField().getZero().add(forward ? minStep : -minStep);
            } else {
                throw new MathIllegalArgumentException(LocalizedODEFormats.MINIMAL_STEPSIZE_REACHED_DURING_INTEGRATION,
                                                       FastMath.abs(h.getReal()), minStep, true);
            }
        }

        if (filteredH.subtract(maxStep).getReal() > 0) {
            filteredH = h.getField().getZero().add(maxStep);
        } else if (filteredH.add(maxStep).getReal() < 0) {
            filteredH = h.getField().getZero().add(-maxStep);
        }

        return filteredH;

    }

    /** Set the initial step size.
     * <p>This method allows the user to specify an initial positive
     * step size instead of letting the integrator guess it by
     * itself. If this method is not called before integration is
     * started, the initial step size will be estimated by the
     * integrator.</p>
     * @param initialStepSize initial step size to use (must be positive even
     * for backward integration ; providing a negative value or a value
     * outside of the min/max step interval will lead the integrator to
     * ignore the value and compute the initial step size by itself)
     */
    public void setInitialStepSize(final double initialStepSize) {
        if ((initialStepSize < minStep) || (initialStepSize > maxStep)) {
            initialStep = -1.0;
        } else {
            initialStep = initialStepSize;
        }
    }

    /** Get the initial step.
     * @return initial step
     */
    public double getInitialStep() {
        return initialStep;
    }

    /** Get the minimal step.
     * @return minimal step
     */
    public double getMinStep() {
        return minStep;
    }

    /** Get the maximal step.
     * @return maximal step
     */
    public double getMaxStep() {
        return maxStep;
    }

    /** Get a dummy step size.
     * @return geometric mean of {@link #getMinStep()} and {@link #getMaxStep()}
     */
    public double getDummyStepsize() {
        return FastMath.sqrt(minStep * maxStep);
    }

}
