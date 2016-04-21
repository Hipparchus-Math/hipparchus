# Building Hipparchus

Hipparchus is a pure java library without any dependencies,
split into several modules. It is designed to be built using
[maven](https://maven/apache/org/) but it can also be built
using an IDE like [Eclipse](http://www.eclipse.org),
[NetBeans](https://netbeans.org), [IntelliJ IDEA](https://www.jetbrains.com/idea/)
or any other. We will describe only a few possibilities here.

### Maven build

As Hipparchus is a multi-modules maven project, its build
is slightly different from a single module project. Once
you have retrieved the complete source tree, either by
cloning the git repository or by unpacking a source distribution,
you will find a parent folder containing a few files and
sub-folders:

  * `LICENSE.txt`
  * `NOTICE.txt`
  * `README.md`
  * `pom.xml`
  * `hipparchus-aggregator`
  * `hipparchus-clustering`
  * `hipparchus-core`
  * `hipparchus-fft`
  * `hipparchus-fitting`
  * `hipparchus-genetics`
  * `hipparchus-geometry`
  * `hipparchus-migration`
  * `hipparchus-ode`
  * `hipparchus-optim`
  * `hipparchus-perf`
  * `hipparchus-release`
  * `hipparchus-samples`
  * `src`

This parent folder contains the maven parent project that
defines elements that will be inherited from the modules
projects. The parent project does not contain any code by itself
(even the `src` folder at this level does not contain code).

The various `hipparchus-xyz` folders are the individual
projects that inherit from the parent.

One specific project is particularly important for building
all Hipparchus modules: the `hipparchus-aggregator` project.
As its name implies, it is a maven `aggregator` project. This
project does not either contain code by itself, but it
defines the organisation of Hipparchus into modules (i.e.
it points to its sibblings `hipparchus-core`, `hipparchus-geometry`,
...) projects and it can build all of them in a consistent
order.

So in order to build Hipparchus, you should go to the
`hipparchus-aggregator` folder and run the following
maven command:

    mvn package

This command will build every single module in the
appropriate order and will move the built artifacts
into the aggregator `target` folder that will be created
automatically. Once the command completes, you will
find the various jar files there.

### Eclipse build

If you want to build Hipparchus using the Eclipse IDE, the
simplest way is to rely on the `m2e` plugin (which is
automatically included in recent Eclipse IDE for Java
developers) as this plugin automatically integrate
maven projects and supports multi-module projects.

Just as in the maven build, you should first retrieve
the Hipparchus source tree (Eclipse is also able to
retrieve this tree directly from git if you want). Once
the source tree is available, you should select
`Import...` entry in the context menu in the Package
Explorer tab in the left hand side view in the Java
perspective, and select `Maven -> Existing Maven Projects`
in the wizard. Then in the `Import Maven Project` wizard,
when selecting the root folder, browse to the
_hipparchus-aggregator_ folder, __not__ to the
parent hipparchus folder. The reason you should select
`hipparchus-aggregator` is because it is the project
that references all the modules, so Eclipse will be
aware of them and will be able to configure everything
properly, including the inter-module dependencies. Once
the `hipparchus-aggregator` project as been selected,
the wizard should show all the projects at once. Pressing
the `Finish` button will configure and compile everything.
