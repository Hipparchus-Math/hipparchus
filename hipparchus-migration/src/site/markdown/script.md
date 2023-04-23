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
## Migration script

There is a python migration script in the `hipparchus-migration/scripts` directory.

The `migrate.py` python migration script can run on the source directory of an
application currently depending on Apache Commons Math in order to update the
sources to rely on Hipparchus instead.

The script has a `--help` option to display its arguments:

    luc@lehrin% script=/path/to/hipparchus/hipparchus-migration/scripts/migrate.py
    luc@lehrin% python $script -h
    usage: migrate [-h] --dir DIR [--ignore IGNORE] --ext EXT [--nosave]
                   [--from-prefix FROM_PREFIX] [--to-prefix TO_PREFIX]
                   [--classes-subst CLASSES_SUBST] [--dry-run] [--verbose]

    optional arguments:
      -h, --help            show this help message and exit
      --dir DIR
      --ignore IGNORE
      --ext EXT
      --nosave
      --from-prefix FROM_PREFIX
      --to-prefix TO_PREFIX
      --classes-subst CLASSES_SUBST
      --dry-run
      --verbose, -v
    luc@lehrin%

The various options of the script are:

  * `-h` or `--help`: help
  * `--dir`: specify the top directory containing sources to modify
  * `--ignore`: directory to ignore (typically `--ignore .git` or `--ignore .svn`)
  * `--ext`: extension of the files to process (typically it will be `--ext .java`)
  * `--nosave`: if specified, original files will be removed instead of being saved with a `.orig` extension,
              it is typically used when the source code is already safe in a source code management system
              like git or subversion
  * `--from-prefix`: top level package prefix of classes to substitute from
                     (default is both `org.apache.commons.math3` and `org.apache.commons.math4`)
  * `--to-prefix`: top level package prefix of classes to substitute to
                   (default is `org.hipparchus`)
  * `--classes-subst`: substitution names for classes, default value is `classes-update-deprecated-exceptions.subst`
  * `--dry-run`: allow to see what would be changed, but without changing anything
  * `-v` or `--verbose`: if specified, display messages about files changed and not changed 

A typical run is shown below:

    luc@lehrin:~/sources/eclipse/orekit$ python ../hipparchus/hipparchus-migration/scripts/migrate.py --dir src --ignore .git --ext .java --ext .md --ext .txt --verbose
    Processing files in dir src
    processing file src/site/markdown/building.md -> unchanged
    processing file src/site/markdown/contact.md -> unchanged
    processing file src/site/markdown/faq.md -> unchanged
    processing file src/site/markdown/configuration.md -> unchanged
    processing file src/site/markdown/downloads.md -> unchanged
    processing file src/site/markdown/contributing.md -> unchanged
    processing file src/site/markdown/training.md -> unchanged
    processing file src/site/markdown/sources.md -> unchanged
       ... snip ...
    processing file src/main/java/org/orekit/orbits/EquinoctialOrbit.java -> changed
    processing file src/main/java/org/orekit/orbits/CartesianOrbit.java -> changed
    processing file src/main/java/org/orekit/orbits/KeplerianOrbit.java -> changed
    processing file src/main/java/org/orekit/orbits/PositionAngle.java -> unchanged
    processing file src/main/java/org/orekit/orbits/CircularOrbit.java -> changed
    processing file src/main/java/org/orekit/orbits/Orbit.java -> changed
    processing file src/main/java/org/orekit/attitudes/TabulatedProvider.java -> unchanged
    processing file src/main/java/org/orekit/attitudes/NadirPointing.java -> changed
    processing file src/main/java/org/orekit/attitudes/TabulatedLofOffset.java -> unchanged
       ... snip ...
    processing file src/tutorials/java/fr/cs/examples/propagation/VisibilityCheck.java -> changed
    processing file src/tutorials/java/fr/cs/examples/propagation/EphemerisMode.java -> changed
    processing file src/tutorials/java/fr/cs/examples/propagation/MasterMode.java -> changed
    processing file src/tutorials/java/fr/cs/examples/propagation/TrackCorridor.java -> changed
    processing file src/tutorials/java/fr/cs/examples/propagation/VisibilityCircle.java -> changed
    luc@lehrin%

The script is not foolproof, but does most of the boring job automatically.

There are two classes substitution patterns available: `classes-update-deprecated-exceptions.subst`
and `classes-keep-deprecated-exceptions.subst`. The default value is the first one (`update`), which
should cover most regular use cases. This default substitutions file replaces all occurrences of
exceptions that have been deprecated in Hipparchus by the new simpler exceptions (`MathIllegalArgumentException`,
`MathIllegalStateException` and `MathRuntimeException`). This is suitable for applications that
only catch the exceptions declared by the mathematical library and do not throw these exceptions
by themselves. It is not suitable for applications that mainly create and throw such exceptions by
themselves (typically when they extend the library classes) and that do rarely catch them. The reason
is that the constructors signatures of the replacement classes do not match the constructors signature
of the deprecated exceptions. For such applications, setting the `--classes-subst` option to
`classes-keep-deprecated-exceptions.subst` may be desirable. With this setting, the deprecated exceptions
will be preserved, as they are indeed available in the `hipparchus-migration` jar. So instead
of compilation errors due to wrong constructor arguments, migrated code will have deprecation
warnings as they use directly these deprecated exceptions. Applications that both create and catch
exceptions from the mathematical library should rather rely on the default (`update`) patterns
and fix the compilation errors in the constructors rather than keeping the deprecated exceptions.
The reason is that the deprecated exceptions are _not_ thrown anymore by Hipparchus and hence the
`catch` statements in application code would never be triggered, resulting in a different code path
being followed in application codes. In other words, the `classes-keep-deprecated-exceptions.subst`
pattern is really intended for expert use.