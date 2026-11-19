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
## Building Hipparchus

Hipparchus is a pure java library without any dependencies,
split into several modules. It is designed to be built using
[maven](https://maven/apache/org/) but it can also be built
using an IDE like [Eclipse](http://www.eclipse.org),
[NetBeans](https://netbeans.org), [IntelliJ IDEA](https://www.jetbrains.com/idea/)
or any other. We will describe only a few possibilities here.

## Maven build

As Hipparchus is a multi-module maven project, its build
is slightly different from a single module project. Once
you have retrieved the complete source tree, either by
cloning the git repository or by unpacking a source distribution,
you will find a parent folder containing a few files and
sub-folders:

  * `LICENSE.txt`
  * `NOTICE.txt`
  * `README.md`
  * `CONTRIBUTING.md`
  * `release.sh`
  * `pom.xml`
  * `hipparchus-clustering`
  * `hipparchus-core`
  * `hipparchus-coverage`
  * `hipparchus-fft`
  * `hipparchus-filtering`
  * `hipparchus-fitting`
  * `hipparchus-geometry`
  * `hipparchus-ode`
  * `hipparchus-optim`
  * `hipparchus-parent`
  * `hipparchus-samples`
  * `hipparchus-stat`
  * `src`

This top level folder contains the maven aggregator project that
defines project organization into modules. The aggregator project does
not contain any code by itself (even the `src` folder at this level
does not contain code).

The `hipparchus-parent` folder contains the parent project for other
modules. This parent project defines elements that will be inherited
by all modules like maven plugins version numbers and configuration.
This project does not contain code by itself.

The `hipparchus-coverage` folder contains a dummy module for aggregating
coverage reports. This project does not contain code by itself.

The remaining `hipparchus-xyz` folders are the projects for the various
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

If you also want to build the static site containing the project
documentation (general project information, user guide, javadocs,
source reference ...), you will have to run the following maven
command:

    mvn site site:stage

The first goal (`site`) will create separate sub-sites in all modules
(in `target/site` folders), including the aggregator module.  All
inter-module links in these temporary folders fail because they are
not combined. The second goal (`site:stage`) combines all these
sub-sites together in one consistent and fully linked site in a
single staging folder common to all modules. This single staging
folder is created as `target/staging` folder in the aggregator
project. Beware that you should _not_ run `mvn clean site site:stage`,
i.e. the `clean` should not be used in the same command line as
the `site:stage` goal, otherwise maven would clean the staged sites
for the previous modules just before processing the final aggregator
module, so only the top level part of the site would be accessible.

## Eclipse build

If you want to build Hipparchus using the Eclipse IDE, the
simplest way is to rely on the `m2e` plugin (which is
automatically included in recent Eclipse IDE for Java
developers) as this plugin automatically integrates
maven projects and supports multi-module projects.

Just as in the maven build, you should first retrieve the Hipparchus
source tree (Eclipse is also able to retrieve this tree directly from
git if you want). Once the source tree is available, you should select the
`Import...` entry in the context menu in the Package Explorer tab in
the left hand side view in the Java perspective, and select `Maven ->
Existing Maven Projects` in the wizard. Then in the `Import Maven
Project` wizard, when selecting the root folder, browse to the top
level hipparchus folder. This folder contains the aggregator project
that references all the modules, so Eclipse will be aware of them and
will be able to configure everything properly, including the
inter-module dependencies. Once the `hipparchus-aggregator` project has
been selected, the wizard should show all the projects at
once. Pressing the `Finish` button will configure and compile
everything.
