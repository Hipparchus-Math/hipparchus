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
  * `hipparchus-clustering`
  * `hipparchus-core`
  * `hipparchus-fft`
  * `hipparchus-fitting`
  * `hipparchus-genetics`
  * `hipparchus-geometry`
  * `hipparchus-migration`
  * `hipparchus-ode`
  * `hipparchus-optim`
  * `hipparchus-parent`
  * `hipparchus-perf`
  * `hipparchus-release`
  * `hipparchus-samples`
  * `src`

This top level folder contains the maven aggregator project that
defines project organization into modules. The aggregator project does
not contain any code by itself (even the `src` folder at this level
does not contain code).

The `hipparchus-parent` folder contains the other models parent
project that defines elements that will be inherited by all modules
like maven plugins version numbers and configuration. This project
does not either contain code by itself.

The other `hipparchus-xyz` folders are the projects for the various
modules that compose Hipparchus and that inherit from the parent.

In order to build Hipparchus, you should stay in the top level folder
(above the various `hipparchus-xyz` sub-folders) and run the following
maven command:

    mvn package

This command will build every single module in the appropriate order
and will assemble the built artifacts into the `target` folder that will
be created automatically. Once the command completes, you will find
several archive files in `zip` and `tar.bz2` formats that contain
all artifacts. The archive files with

### Eclipse build

If you want to build Hipparchus using the Eclipse IDE, the
simplest way is to rely on the `m2e` plugin (which is
automatically included in recent Eclipse IDE for Java
developers) as this plugin automatically integrate
maven projects and supports multi-module projects.

Just as in the maven build, you should first retrieve the Hipparchus
source tree (Eclipse is also able to retrieve this tree directly from
git if you want). Once the source tree is available, you should select
`Import...` entry in the context menu in the Package Explorer tab in
the left hand side view in the Java perspective, and select `Maven ->
Existing Maven Projects` in the wizard. Then in the `Import Maven
Project` wizard, when selecting the root folder, browse to the top
level hipparchus folder. This folder contains the aggregator project
that references all the modules, so Eclipse will be aware of them and
will be able to configure everything properly, including the
inter-module dependencies. Once the `hipparchus-aggregator` project as
been selected, the wizard should show all the projects at
once. Pressing the `Finish` button will configure and compile
everything.
