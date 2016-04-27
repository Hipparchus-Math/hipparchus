## Migration script

There is a python migration script in the `hipparchus-migration/scripts` directory called `migrate.py`.

The migration script can run on the source directory of an application currently depending on
Apache Commons Math in order to update the sources to rely on Hipparchus instead.

The script has a `--help` option to display its arguments:

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
    luc@lehrin%

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

A typical run is shown below:

    luc@lehrin% script=/path/to/hipparchus/hipparchus-migration/scripts/migrate.py
    luc@lehrin% python $script --dir src --ext .java --ignore .git --verbose
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
    luc@lehrin%

The script is not foolproof, but does most of the boring job automatically.