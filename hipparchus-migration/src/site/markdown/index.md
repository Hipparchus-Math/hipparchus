
## Hipparchus Migration

This module is devoted to migrating from [Apache Commons Math](https://commons.apache.org/math/)
to [Hipparchus](../index.html).

The module provides a transition script and a temporary migration jar. The migration
process is composed of several phases:

  * run the migration script on the application source
  * modify the build configuration to use Hipparchus modules instead of Apache Commons Math jar
  * fix the (hopefully few) compilation errors that were not handled by the script
  * check everything runs, ignoring the warnings about deprecated interfaces that are
    provided by the intermediate hipparchus-migration jar
  * progressively get rid of the deprecated interfaces, following the deprecation statements hints
  * remove the dependency to the intermediate hipparchus-migration jar

### Running the migration script

The script can run on the source directory of an application currently depending on
Apache Commons Math in order to update the sources to rely on Hipparchus instead.

A typical way to run the script would be:

    luc@lehrin% script=/path/to/hipparchus/hipparchus-migration/scripts/migrate.py
    luc@lehrin% python $script -h
    usage: migrate [-h] --dir DIR [--ignore IGNORE] --ext EXT [--nosave]
               [--classes-subst CLASSES_SUBST] [--dry-run] [--verbose]

    optional arguments:
      -h, --help            show this help message and exit
      --dir DIR
      --ignore IGNORE
      --ext EXT
      --nosave
      --classes-subst CLASSES_SUBST
      --dry-run
      --verbose, -v
    (lehrin) luc% python $script --dir src --ext .java --ignore .git --verbose
    Processing files in dir src
    processing file src/test/java/org/orekit/Utils.java -> unchanged
    processing file src/test/java/org/orekit/KeyValueFileParser.java -> changed
    processing file src/test/java/org/orekit/OrekitMatchers.java -> changed
    processing file src/test/java/org/orekit/SolarInputs97to05.java -> changed
    processing file src/test/java/org/orekit/models/earth/FixedTroposphericModelTest.java -> changed
    processing file src/test/java/org/orekit/models/earth/SaastamoinenModelTest.java -> changed
    processing file src/test/java/org/orekit/models/earth/KlobucharModelTest.java -> changed
    processing file src/test/java/org/orekit/models/earth/ReferenceEllipsoidTest.java -> changed
    processing file src/test/java/org/orekit/models/earth/EarthITU453AtmosphereRefractionTest.java -> changed
       ... snip ...
    processing file src/main/java/org/orekit/utils/SecularAndHarmonic.java -> changed
    processing file src/main/java/org/orekit/utils/LoveNumbers.java -> unchanged
    (lehrin) luc%

The various options of the script are:

  * `--help`: help
  * `--dir`: specify the top directory containing sources to modify
  * `--ext`: extension of the files to process (typically it will be `--ext .java`)
  * `--ignore`: directory to ignore (typically `--ignore .git` or `--ignore .svn`)
  * `--classes-subst` : substitution names for classes, generally not specified to use the one provided by the module
  * `--dry-run` : allow to see what would be changed, but without changing anything
  * `--nosave` : if specified, original files will be removed instead of being saved with a `.orig` extension,
               it is typically used when the source code is already safe in a source code management system
               like git or subversion
  * `--verbose` : if specified, display more messages when doing the modifications 

The script is not foolproof, but does most of the boring job automatically.

### Modifying the build configuration

As your application depended on Apache Commons Math which as only one jar, your build
system knows about this jar. It has to be replaced by the various jars provided by
Hipparchus. All uses will at least need the hipparchus-core jar. Complex uses may
need all the jars.

The simplest way to identify which jars you need is to add them one by one, looking
at the compilation errors. This is easy and fast as there are only a dozen jars
and you will certainly not depend on some of them (hipparchus-perf and hipparchus-samples
are typically not used). The first jar that will be mandatory is hipparchus-core. If this
one is not sufficient, the second one to add is hipparchus-migration because it adds some
glue code for migration that is often needed. Then add other jars depending on the compilation
errors you get, looking only at the package names of the missing classes.

### Fixing compilation errors

There may be some problems not handled by the script. The known problems
are the following ones:

  * `ExceptionContext` not found

### Checking everything runs

Once everything compiles, you should run your application tests to see if everything
runs. If not, ask to the [developers mailing lists](../mail-lists.html).

Deprecation warnings at this step are normal. They correspond to APIs that were present
in Apache Commons Math and have been modified in Hipparchus. In order to help users
switch, we have moved these APIs in the hipparchus-migration jar and they act as
facades between user code and the new Hipparchus APIs. They are marked as deprecated
so users can easily spot them and change to the newer APIs progressively while still
having a running application.

### Getting rid of the deprecated interfaces

When the application runs in a full Hipparchus environment, you can start updating
your code to use the new APIs provided by Hipparchus. This is just a traditional
maintenance task. The deprecation hints in the hipparchus-migration sources will
help you find the new interfaces. They are often very close to the previous
ones and modification is expected to be straightforward.

### Removing the dependency to the intermediate hipparchus-migration jar

When all deprecation warnings have been fixed, you can safely remove the
hipparchus-migration jar from your build system, everything should work.
