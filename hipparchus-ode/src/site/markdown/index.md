<!--
 Licensed to the Hipparchus project under one or more
 contributor license agreements.  See the NOTICE file distributed with
 this work for additional information regarding copyright ownership.
 The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 (the "License"); you may not use this file except in compliance with
 the License.  You may obtain a copy of the License at

      https://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
-->
# Ordinary Differential Equations Integration
## Overview
The ode package provides classes to solve Ordinary Differential Equations problems.

This package solves Initial Value Problems of the form \\(y'=f(t,y)\\) with \\(t_0\\)
and \\(y(t_0)=y_0\\) known. The \\(f\\) function is the
[Ordinary Differential Equation](../apidocs/org/hipparchus/ode/OrdinaryDifferentialEquation.html)
to be solved. The provided integrators compute an estimate of \\(y(t)\\) from \\(t=t_0\\)
to \\(t=t_1\\).

All integrators provide [dense output](#Dense_Output). This means that besides computing the state vector
at discrete times, they also provide a cheap mean to get both the state and its derivative
between the time steps. They do so through classes extending the
[ODEStateInterpolator](../apidocs/org/hipparchus/ode/sampling/ODEStateInterpolator.html)
abstract class, which are made available to the user at the end of each step.

All integrators handle multiple [discrete events detection](#Discrete_Events_Handling) based on switching
functions. This means that the integrator can be driven by user specified discrete events
(occurring when the sign of user-supplied _switching function_ changes). The steps are
shortened as needed to ensure the events occur at step boundaries (even if the integrator
is a fixed-step integrator). When the events are triggered, integration can
be stopped (this is called a G-stop facility), the state vector can be changed, or integration
can simply go on. The latter case is useful to handle discontinuities in the differential
equations gracefully and get accurate dense output even close to the discontinuity.

All integrators support setting a maximal number of evaluations of differential
equations function. If this number is exceeded, an exception will be thrown during
integration. This can be used to prevent infinite loops if for example error control or
discrete events create a really large number of extremely small steps. By default, the
maximal number of evaluation is set to `Integer.MAX_VALUE` (i.e. \\(2^{31}-1\\)
or 2147483647). It is recommended to set this maximal number to a value suited to the ODE
problem, integration range, and step size or error control settings.

All integrators support expanding the primary ODE with one or more [secondary ODE](#Secondary_States) to manage
additional states that will be integrated together with the primary state. This can be used
for example to compute [partial derivatives](#Derivatives) of the primary state with respect to either
the initial state or some parameters (see below for an example).

Two parallel APIs are available. The first is devoted to solve ode for which the integration free
variable t and the state y(t) are primitive double and primitive double array respectively. The
second API is devoted to solve ode for which the integration free variable t and the state y(t)
are [CalculusFieldElement](../apidocs/org/hipparchus/CalculusFieldElement.html) and
[CalculusFieldElement](../apidocs/org/hipparchus/CalculusFieldElement.html) array respectively. This
allow for example users to integrate ode where the computation values
are for example [DerivativeStructure](../apidocs/org/hipparchus/analysis/differentiation/DerivativeStructure.html)
elements, hence automatically computing partial derivatives with respect to some equations parameters
without a need to set up the [variational equations](../apidocs/org/hipparchus/ode/VariationalEquation.html)
manually or if derivatives with orders higher than 1 is needed, which cannot be handled by
[variational equations](../apidocs/org/hipparchus/ode/VariationalEquation.html). Another example is to use
[Dfp](../apidocs/org/hipparchus/analysis/dfp/Dfp.html) elements in order to solve ode with extended precision.

## Basic Use
The user should describe his problem in his own classes which should implement the
[OrdinaryDifferentialEquation](../apidocs/org/hipparchus/ode/OrdinaryDifferentialEquation.html)
interface (or  [FieldOrdinaryDifferentialEquation](../apidocs/org/hipparchus/ode/FieldOrdinaryDifferentialEquation.html)
interface). The following example shows how to implement a simple two-dimensional problem using double primitives:
\\[
  \\begin{cases}
   y'_0 (t) &= \\omega \\times (c_1 - y_1 (t))\\\\
   y'_1 (t) &= \\omega \\times (y_0 (t) - c_0)
  \\end{cases}
\\]

with some initial state \\(y(t_0) = (y_0(t_0), y_1(t_0))\\).
In fact, the exact solution of this problem is that \\(y(t)\\) moves along a circle
centered at \\(c = (c_0, c_1)\\) with constant angular rate \\(\omega\\).


    private static class Circle1 implements OrdinaryDifferentialEquation {
    
        private double[] c;
        private double omega;
    
        public Circle1(double[] c, double omega) {
            this.c     = c;
            this.omega = omega;
        }
    
        public int getDimension() {
            return 2;
        }
    
        public double[] computeDerivatives(double t, double[] y) {
            return new double[] {
                omega * (c[1] - y[1]),
                omega * (y[0] - c[0])
            };
        }
    
    }

The [OrdinaryDifferentialEquation](../apidocs/org/hipparchus/ode/OrdinaryDifferentialEquation.html)
interface also defines an `init` method that can be used is some initialization must
be performed at the start of the integration. There is a default implementation that does nothing,
so the method must be implemented only if really needed.

The user must also set up an initial state as an [ODEState](../apidocs/org/hipparchus/ode/ODEState.html) instance.
Then he should pass both the equations and the initial state to the integrator he prefers among all the classes that implement
the [ODEIntegrator](../apidocs/org/hipparchus/ode/ODEIntegrator.html)
interface (or the [FieldODEIntegrator](../apidocs/org/hipparchus/ode/FieldODEIntegrator.html)
interface). Computing the state \\(y(16.0)\\) starting from \\(y(0.0) = (0.0, 1.0)\\) and integrating the ODE
is done as follows (using Dormand-Prince 8(5,3) integrator as an example):

    ODEIntegrator                dp853        = new DormandPrince853Integrator(1.0e-8, 100.0, 1.0e-10, 1.0e-10);
    OrdinaryDifferentialEquation ode          = new Circle1(new double[] { 1.0, 1.0 }, 0.1);
    ODEState                     initialState = new ODEState(0.0, new double[] { 0.0, 1.0 });
    ODEStateAndDerivative        finalState   = dp853.integrate(ode, initialState, 16.0);
    double                       t            = finalState.getTime();
    double[]                     y            = finalState.getPrimaryState();
    System.out.format(Locale.US, "final state at %4.1f: %6.3f %6.3f%n", t, y[0], y[1]);

## Dense Output
The solution of the integration problem is provided by two means. The first one, shown in the previous
example, is aimed towards simple use: the state vector at the end of the integration process is returned
by the `ODEIntegrator.integrate` method, in an
[ODEStateAndDerivative](../apidocs/org/hipparchus/ode/ODEStateAndDerivative.html) instance. The second one
should be used when more in-depth information is needed throughout the integration process and not
only at the end. The user can register an object implementing the
[ODEStepHandler](../apidocs/org/hipparchus/ode/sampling/ODEStepHandler.html) interface or a
[StepNormalizer](../apidocs/org/hipparchus/ode/sampling/StepNormalizer.html) object wrapping
a user-specified object implementing the
[ODEFixedStepHandler](../apidocs/org/hipparchus/ode/sampling/ODEFixedStepHandler.html) interface
into the integrator before calling the `ODEIntegrator.integrate` method. The user object
will be called appropriately during the integration process, allowing the user to monitor intermediate
results. The default step handler does nothing. Considering again the previous example, we want to print the
trajectory of the point to check it really is a circle arc. We simply add the following before the call
to `dp853.integrate`:


    ODEStepHandler stepHandler = new ODEStepHandler() {
        public void handleStep(ODEStateInterpolator interpolator, boolean isLast) {
            double stepStart = interpolator.getPreviousState().getTime();
            double stepEnd   = interpolator.getPreviousState().getTime();
            for (int i = 0; i < 20; ++i) {
                // we want to print 20 points for each step
                double t = ((20 - i) * stepStart + i * stepEnd) / 20;
                double[] y = interpolator.getInterpolatedState(t).getPrimaryState();
                System.out.println(t + " " + y[0] + " " + y[1]);
            }
        }
    };
    dp853.addStepHandler(stepHandler);

The [ODEStepHandler](../apidocs/org/hipparchus/ode/sampling/ODEStepHandler.html)
interface also defines an `init` method that can be used is some initialization must
be performed at the start of the integration. There is a default implementation that does nothing,
so the method must be implemented only if really needed.

[DenseOutputModel](../apidocs/org/hipparchus/ode/DenseOutputModel.html)
is a special-purpose step handler that is able to store all steps and to provide transparent access to
any intermediate result once the integration is over. An important feature of this class is that it
implements the `Serializable` interface. This means that a complete continuous model of the
integrated function throughout the integration range can be serialized and reused later (if stored into
a persistent medium like a file system or a database) or elsewhere (if sent to another application).
Only the result of the integration is stored, there is no reference to the integrated problem by itself.

Another predefined implementations of the [ODEStepHandler](../apidocs/org/hipparchus/ode/sampling/ODEStepHandler.html)
interface, [StepNormalizer](../apidocs/org/hipparchus/ode/sampling/StepNormalizer.html)) is devoted
to normalize steps to a fixed size even if the integrator is a variable step integrator. This allows for
example to print a regular ephemeris for a trajectory.

Specific implementations can be developed by users for specific needs. As an example, if an application
is to be completely driven by the integration process, then most of the application code will be
run inside a step handler specific to this application.

Some integrators (the simple ones) use fixed steps that are set at creation time. The more efficient
integrators use variable steps that are handled internally in order to control the local integration error
of the main state with respect to a specified accuracy (these integrators extend the
[AdaptiveStepsizeIntegrator](../apidocs/org/hipparchus/ode/nonstiff/AdaptiveStepsizeIntegrator.html)
abstract class). The secondary equations are explicitly ignored for step size control, in order to get reproducible
results regardless of the secondary equations being integrated or not. The step handler which is called after each
successful step shows up the variable stepsize. The [StepNormalizer](../apidocs/org/hipparchus/ode/sampling/StepNormalizer.html)
class can be used to convert the variable stepsize into a fixed stepsize that can be handled by classes
implementing the [FixedStepHandler](../apidocs/org/hipparchus/ode/sampling/FixedStepHandler.html)
interface. Adaptive stepsize integrators can automatically compute the initial stepsize by themselves,
however the user can specify it if he prefers to retain full control over the integration or if the
automatic guess is wrong.


## Discrete Events Handling
ODE problems are continuous ones. However, sometimes discrete events may be
taken into account. The most frequent case is the stop condition of the integrator
is not defined by the time t but by a target condition on state \\(y\\) (say \\(y[0] = 1.0\\)
for example).

Discrete events detection is based on switching functions. The user provides
a simple [g(state)](../apidocs/org/hipparchus/ode/events/ODEEventHandler.html)
switching function depending on the current state. The integrator will monitor
the value of the function throughout integration range and will trigger the
event when its sign changes. The magnitude of the value is almost irrelevant.
For the sake of root finding, it should however be continuous (but not necessarily smooth)
at least in the roots vicinity. The steps are shortened as needed to ensure the events occur
at step boundaries (even if the integrator is a fixed-step integrator).

When an event is triggered (i.e. when the sign of the \\(g\\) switching function changes),
the event state and an indicator whether the switching function was increasing or
decreasing at event occurrence are provided to the user. Several different options are
available to deal with the event occurrence.

The first case, G-stop, is the most common one. A typical use case is when an
ODE must be integrated until some target state is reached, with a known value of
the state but an unknown occurrence time. As an example, if we want to monitor
a chemical reaction until some predefined concentration for the first substance,
we can use the following switching function setting:


    public double g(ODEStateAndDerivative state) {
        return state.getPrimaryState()[0] - targetConcentration;
    }
    
    public Action eventOccurred(ODEStateAndDerivative state, boolean increasing) {
        return Action.STOP;
    }

The second case correspond to discontinuous dynamical models, for which state vector
or derivatives must be changed instantaneously when the event occurs. A typical case
would be the motion of a spacecraft when thrusters are fired for orbital maneuvers.
The acceleration is smooth as long as no maneuvers are performed, depending only on
gravity, drag, third body attraction, radiation pressure. Firing a thruster introduces a
discontinuity on the acceleration that must be handled appropriately by the integrator.
In such a case, we would use a switching function setting similar to this:


    public double g(ODEStateAndDerivative state) {
        final double t = state.getTime();
        return (t - tManeuverStart) * (t - tManeuverStop);
    }
    
    public Action eventOccurred(ODEStateAndDerivative state, boolean increasing) {
      return Action.RESET_DERIVATIVES;
    }

The third case is useful mainly for monitoring purposes, when we simply want
to log the event but otherwise do not interfere with the integration process.
A simple example is:


    public double g(ODEStateAndDerivative state) {
        final double[] y = state.getPrimaryState();
        return y[0] - y[1];
    }
    
    public Action eventOccurred(ODEStateAndDerivative state, boolean increasing) {
        logger.log("y0(t) and y1(t) curves cross at t = " + v state.getTime());
        return Action.CONTINUE;
    }

The difference between the cases appears in the return value of the
`eventOccurred` method. The integrator will call the user method when
integration process reaches a point where the `g` switching function
sign changes, then it will call the `eventOccurred` method and its
behaviour will depend on the value returned by this call.

## Available Integrators
The tables below show the various integrators available for non-stiff problems. Note that the
implementations of Adams-Bashforth and Adams-Moulton are adaptive stepsize, not fixed stepsize
as is usual for these multi-step integrators. This is due to the fact the implementation relies
on the Nordsieck vector representation of the state.

| <font size="+1">Fixed Step Integrators</font> |
| --- |
| Name | Order |
| [Euler](../apidocs/org/hipparchus/ode/nonstiff/EulerIntegrator.html) | 1 |
| [Midpoint](../apidocs/org/hipparchus/ode/nonstiff/MidpointIntegrator.html) | 2 |
| [Classical Runge-Kutta](../apidocs/org/hipparchus/ode/nonstiff/ClassicalRungeKuttaIntegrator.html) | 4 |
| [Gill](../apidocs/org/hipparchus/ode/nonstiff/GillIntegrator.html) | 4 |
| [3/8](../apidocs/org/hipparchus/ode/nonstiff/ThreeEighthesIntegrator.html) | 4 |
| [Luther](../apidocs/org/hipparchus/ode/nonstiff/LutherIntegrator.html) | 6 |

| <font size="+1">Adaptive Stepsize Integrators</font> |
| --- |
| Name | Integration Order | Error Estimation Order |
| [Higham and Hall](../apidocs/org/hipparchus/ode/nonstiff/HighamHall54Integrator.html) | 5 | 4 |
| [Dormand-Prince 5(4)](../apidocs/org/hipparchus/ode/nonstiff/DormandPrince54Integrator.html) | 5 | 4 |
| [Dormand-Prince 8(5,3)](../apidocs/org/hipparchus/ode/nonstiff/DormandPrince853Integrator.html) | 8 | 5 and 3 |
| [Gragg-Bulirsch-Stoer](../apidocs/org/hipparchus/ode/nonstiff/GraggBulirschStoerIntegrator.html) | variable (up to 18 by default) | variable |
| [Adams-Bashforth](../apidocs/org/hipparchus/ode/nonstiff/AdamsBashforthIntegrator.html) | variable | variable |
| [Adams-Moulton](../apidocs/org/hipparchus/ode/nonstiff/AdamsMoultonIntegrator.html) | variable | variable |


## Secondary States
In some cases, the ordinary differential equations is split into a primary
state and some additional, secondary states. This can be done by providing
an [ExpandableODE](../apidocs/org/hipparchus/ode/ExpandableODE.html)
instance to the integrator instead of a regular
[Ordinary Differential Equation](../apidocs/org/hipparchus/ode/OrdinaryDifferentialEquation.html)
and add [SecondaryODE](../apidocs/org/hipparchus/ode/SecondaryODE.html) to this expandable
equation.

The [ODEState](../apidocs/org/hipparchus/ode/ODEState.html) initial state and
[ODEStateAndDerivative](../apidocs/org/hipparchus/ode/ODEStateAndDerivative.html) final state
contains data pertaining to both the primary state and all the secondary states.

The adaptive stepsize integrators rely only on the primary state for stepsize control by design.
If some components should be used for step size control, then they must be part of the
primary state and not added as secondary states.

## Derivatives
In some cases, the sensitivity \\(dy/dy_0\\) of the final state with respect to the initial
state \\(y_0\\) (often called state transition matrix) or the sensitivity \\(dy/dp_k\\) of
the final state with respect to some parameters \\(p_k\\) of the ODE are needed, in addition
to the final state \\(y(t)\\) itself.

This can be computed by using [variational equations](../apidocs/org/hipparchus/ode/VariationalEquation.html)
which automatically adds a set of secondary equations to the primary ODE and adds initial secondary
states corresponding to the desired partial derivatives before the integration starts. Then the integration
will propagate the compound state composed of both the primary state and its partial derivatives.
At the end of the integration, the Jacobian matrices are extracted from the integrated secondary state.
As the variational equations are considered to be secondary equations here, variable step integrators ignore
them for step size control: they rely only on the main state. This feature is a design choice. The rationale is
to get exactly the same steps, regardless of the Jacobians being computed or not, hence ensuring reproducible
results in both cases.

If for example the primary state dimension is 6 and there are 3 parameters, the compound state will be a 60
elements array. The first 6 elements will be the primary state, the next 54 elements will be a secondary
state containing first the \\(6 \times 6\\) Jacobian matrix of the final state with respect to the initial state, and
3 sets of 6 elements vectors for the Jacobian matrix of the final state with respect to the 3 parameters.

The [VariationalEquation](../apidocs/org/hipparchus/ode/VariationalEquation.html) handles most of
secondary equations and states management, as long as the local partial derivatives are provided to it.
On both initialization of the matrices and extraction of the intermediate or final results from
[ODEStateAndDerivative](../apidocs/org/hipparchus/ode/ODEStateAndDerivative.html) instances, the state
transition matrix is represented as a square \\(n \times n\\) array (for a dimension \\(n\\) primary state).
The parameters are identified by a name (a simple user defined string) and the associated partial
derivatives are dimension \\(n\\) arrays.

What remains of user responsibility is to provide the local Jacobians \\(df(t, y, p)/dy\\) and
\\(df(t, y, p)/dp_k\\) corresponding the the main ODE \\(y'=f(t, y, p)\\). There are two ways for user to
provide theses local partial derivatives. The first way is to provide a full-fledged
[ODEJacobiansProvider](../apidocs/org/hipparchus/ode/ODEJacobiansProvider.html) which is an
[OrdinaryDifferentialEquation](../apidocs/org/hipparchus/ode/OrdinaryDifferentialEquation.html)
with additional capabilities to compute partial derivatives. The second way is to provide
a simple [OrdinaryDifferentialEquation](../apidocs/org/hipparchus/ode/OrdinaryDifferentialEquation.html)
and to configure a way to compute the required derivatives using finite differences.

Note that if the [VariationalEquation](../apidocs/org/hipparchus/ode/VariationalEquation.html) instance
is known by [step handlers](../apidocs/org/hipparchus/ode/sampling/ODEStepHandler.html) or
[event handlers](../apidocs/org/hipparchus/ode/events/ODEEventHandler.html), it is possible to
retrieve the current values of the state transition matrix and the parameters derivatives
during the integration process.


As the variational equation automatically inserts secondary differential equations in the
[expandable ODE](../apidocs/org/hipparchus/ode/ExpandableODE.html), data for initial state
must also be inserted before integration and matrices result must be extracted after integration.
This implies a precise scheduling of the calls to the various methods of this class. The
proper scheduling is the following one:

    // set up equations
    ODEJacobiansProvider jode       = new MyODE(...);
    ExpandableODE        expandable = new Expandable(jode);
    VariationalEquation  ve         = new VariationalEquation(expandable, jode);

    // set up initial state
    ODEState initWithoutDerivatives = new ODEState(t0, y0);
    ve.setInitialMainStateJacobian(dYdY0); // only needed if the default identity matrix is not suitable
    ve.setInitialParameterJacobian(name, dYdP); // only needed if the default zero matrix is not suitable
    ODEState initWithDerivatives = ve.setUpInitialState(initWithoutDerivatives);

    // perform integration on the expanded equations with the expanded initial state
    ODEStateAndDerivative finalState = integrator.integrate(expandable, initWithDerivatives, finalT);

    // extract Jacobian matrices
    dYdY0 = ve.extractMainSetJacobian(finalState);
    dYdP  = ve.extractParameterJacobian(finalState, name);

The most important part is to not forget to call `setUpInitialState` to add the secondary
state with the initial matrices to the `ODEState` used in the `integrate` method. Forgetting to
do this and passing only a `ODEState` without the secondary state set up will trigger an
error as the state vector will not have the correct dimension.

The following example corresponds to a simple case where all derivatives can be computed analytically. The state is
a 2D point traveling along a circle. There are three parameters: the two coordinates of the center and the
angular velocity.


    private static class Circle2 implements ODEJacobiansProvider {

        public static final String CX    = "cx";
        public static final String CY    = "cy";
        public static final String OMEGA = "omega";

        private double cx;
        private double cy;
        private double omega;

        public Circle2(double cx, double cy, double omega) {
            this.cx    = cx;
            this.cy    = cy;
            this.omega = omega;
        }

        public int getDimension() {
            return 2;
        }

        public double[] computeDerivatives(double t, double[] y) {
            return new double[] {
                omega * (cy - y[1]),
                omega * (y[0] - cx)
            };
        }

        public double[][] computeMainStateJacobian(double t, double[] y, double[] yDot) {
            return new double[][] {
                { 0, -omega },
                { omega, 0 }
            };
        }

        public List<String> getParametersNames() {
            return Arrays.asList(CX, CY, OMEGA);
        }

        public boolean isSupported(String name) {
            return CX.equals(name) || CY.equals(name) || OMEGA.equals(name);
        }

        public double[] computeParameterJacobian(double t, double[] y, double[] yDot, String paramName)
            throws MathIllegalArgumentException {
            if (CX.equals(paramName)) {
                return new double[] { 0, -omega };
            } else if (CY.equals(paramName)) {
                return new double[] { omega, 0 };
            }  else if (OMEGA.equals(paramName)) {
                return new double[] { cy - y[1], y[0] - cx };
            } else {
                throw new MathIllegalArgumentException(LocalizedODEFormats.UNKNOWN_PARAMETER,
                                                       paramName);
            }
        }
    }

This ODE is integrated as follows:


    ODEIntegrator               integ        = new DormandPrince54Integrator(1.0e-8, 100.0, 1.0e-10, 1.0e-10);
    Circle2                     circle       = new Circle2(1.0,  1.0, 0.1);
    ExpandableODE               expandable   = new ExpandableODE(circle);
    VariationalEquation         ve           = new VariationalEquation(expandable, circle);
    final ODEState              initialState = ve.setUpInitialState(new ODEState(0, new double[] { 0.0, 1.0 }));
    final ODEStateAndDerivative finalState   = integ.integrate(expandable, initialState, 18 * FastMath.PI);

    double[][] dYdY0 = ve.extractMainSetJacobian(finalState);
    System.out.format(Locale.US, "state transition matrixt at t = %6.3f%n  %6.3f %6.3f%n  %6.3f %6.3f%n",
                      finalState.getTime(),
                      dYdY0[0][0], dYdY0[0][1],
                      dYdY0[1][0], dYdY0[1][1]);

    double[] dYdCX = ve.extractParameterJacobian(finalState, Circle.CX);
    System.out.format(Locale.US, "dY/dCx =   %6.3f %6.3f%n", dYdCX[0], dYdCX[1]);

    double[] dYdCY = ve.extractParameterJacobian(finalState, Circle.CY);
    System.out.format(Locale.US, "dY/dCy =   %6.3f %6.3f%n", dYdCY[0], dYdCY[1]);

    double[] dYdO  = ve.extractParameterJacobian(finalState, Circle.OMEGA);
    System.out.format(Locale.US, "dY/dω  =   %6.3f %6.3f%n", dYdO[0], dYdO[1]);

We could have done a similar display throughout the integration by putting a similar code
in a [step handler](../apidocs/org/hipparchus/ode/sampling/ODEStepHandler.html).
