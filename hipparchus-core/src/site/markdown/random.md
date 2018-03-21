<!--
 Licensed to the Hipparchus project under one or more
 contributor license agreements.  See the NOTICE file distributed with
 this work for additional information regarding copyright ownership.
 The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 (the "License"); you may not use this file except in compliance with
 the License.  You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
-->
# Random Data Generation

## Overview
The Hipparchus random package includes utilities for

* generating random numbers
* generating random vectors
* generating random strings
* generating cryptographically secure sequences of random numbers or strings
* generating random samples and permutations
* analyzing distributions of values in an input file and generating values "like" the values in the file
* generating data for grouped frequency distributions or histograms

The source of random data used by the data generation utilities is
pluggable.  In most cases, the default is a Well generator.  Whenever a default is provided, the javadoc indicates what the default is. Other good PRNGs suitable for Monte-Carlo analysis (but <strong>not</strong>
for cryptography) provided by the library in the  [raondom](../apidocs/org/hipparchus/random/package.html).

Sections 2.2-2.6 below show how to use the Hipparchus API to generate
different kinds of random data.  The examples all use the default
JDK-supplied PRNG.  PRNG pluggability is covered in 2.7.  The only
modification required to the examples to use alternative PRNGs is to
replace the argumentless constructor calls with invocations including
a `RandomGenerator` instance as a parameter.


## Random numbers

The [RandomDataGenerator](../apidocs/org/hipparchus/random/RandomDataGenerator.html)
class implements methods for generating random sequences
of numbers. The API contracts of these methods use the following concepts:

### Random sequence of numbers from a probability distribution

There is no such thing as a single "random number."  What can be
generated  are <i>sequences</i> of numbers that appear to be random.  When
using the built-in JDK function `Math.random(),` sequences of
values generated follow the
[Uniform Distribution](http://www.itl.nist.gov/div898/handbook/eda/section3/eda3662.htm),
which means that the values are evenly spread over the interval  between 0 and 1,
with no sub-interval having a greater probability of containing generated values than
any other interval of the same length.  The mathematical concept of a
[probability distribution](http://www.itl.nist.gov/div898/handbook/eda/section3/eda36.htm)
basically amounts to asserting that different ranges in the set of possible values of
a random variable have different probabilities of containing the value. Hipparchus supports
generating random sequences from each of the distributions in the
[distributions](../apidocs/org/hipparchus/distribution/package-summary.html) package.
The javadoc for the `nextXxx` methods in
[RandomDataGenerator](../apidocs/org/hipparchus/random/RandomDataGenerator.html)
describes the algorithms used to generate random deviates.  The `nextXxx` methods allow you to get random deviates directly, without instantiating distributions.  For example, to get a random value following a normal (Gaussian) distribution with mean 3 and standard deviation 1.5, you can use

    RandomDataGenerator randomDataGenerator = new RandomDataGenerator(1000)
    randomDataGenerator.nextNormal(3,1.5)

Here the default Well generator is used the source of randomness and 1000 is passed to it as initial seed.  To generate a sequence of random values to use in a simulation, you should always just create one RandomDataGenerator instance and reuse it.

For user-defined distributions, or Hipparchus distributions not included among the `nextXxx` methods of `RandomDataGenerator`, one can use the `nextDeviate` methods, which take real or integer distribution instances as arguments, implementing a generic inversion-based sampling method for arbitrary distributions.  There are also `nextDeviates` methods that take a distribution and an integer as arguments.  These are handy when you want to generate an array of random values.  Any distribution, including thos covered by the `nextXxx` methods, can be passed to these methods and `RandomDataGenerator` will use the best implementation that it has.

### Cryptographically secure random sequences

It is possible for a sequence of numbers to appear random, but
nonetheless to be predictable based on the algorithm used to generate the
sequence. If in addition to randomness, strong unpredictability is
required, it is best to use a
[secure random number generator](http://www.wikipedia.org/wiki/Cryptographically_secure_pseudo-random_number_generator) to generate values (or strings). The
nextSecureXxx methods implemented by `RandomDataGenerator`
use the JDK `SecureRandom` PRNG to generate cryptographically
secure sequences. The `setSecureAlgorithm` method allows you to
change the underlying PRNG. These methods are __much slower__ than
the corresponding "non-secure" versions, so they should only be used when
cryptographic security is required.  The `ISAACRandom` class implements a
fast cryptographically secure pseudorandom number generator.

### Seeding pseudo-random number generators

By default, the implementation provided in `RandomDataGenerator`
uses a Well19937c generator seeded with the system time and its system identity hashcode. For the non-secure methods, starting with the same seed always produces the
same sequence of values.  Secure sequences started with the same seeds will diverge. To generate sequences of random data values, you should always
instantiate __one__ ` RandomDataGenerator` and use it
repeatedly instead of creating new instances for subsequent values in the
sequence.  For example, the following will generate a random sequence of 50
long integers between 1 and 1,000,000:

    RandomDataGenerator randomData = new RandomDataGenerator(); 
    for (int i = 0; i &lt; 1000; i++) {
        value = randomData.nextLong(1, 1000000);
    }

The following will not in general produce a good random sequence, since the
PRNG is reseeded each time through the loop:

    for (int i = 0; i &lt; 1000; i++) {
        RandomDataGenerator randomData = new RandomDataGenerator(); 
        value = randomData.nextLong(1, 1000000);
    }

The following will produce the same random sequence each time it is executed:

    RandomDataGenerator randomData = new RandomDataGenerator(1000); 
    for (int i = 0; i = 1000; i++) {
        value = randomData.nextLong(1, 1000000);
    }

The following will produce a different random sequence each time it is executed.

    RandomDataGenerator randomData = new RandomDataGenerator(); 
    randomData.reSeedSecure(1000);
    for (int i = 0; i &lt; 1000; i++) {
        value = randomData.nextSecureLong(1, 1000000);
    }


## Random Vectors

Some algorithms require random vectors instead of random scalars. When the
components of these vectors are uncorrelated, they may be generated simply
one at a time and packed together in the vector. The
[UncorrelatedRandomVectorGenerator](../apidocs/org.hipparchus/random/UncorrelatedRandomVectorGenerator.html)
class simplifies this process by setting the mean and deviation of each
component once and generating complete vectors. When the components are
correlated however, generating them is much more difficult. The
[CorrelatedRandomVectorGenerator](../apidocs/org/hipparchus/random/CorrelatedRandomVectorGenerator.html)
class provides this service. In this case, the user must set up a complete
covariance matrix instead of a simple standard deviations vector. This matrix
gathers both the variance and the correlation information of the probability law.

The main use for correlated random vector generation is for Monte-Carlo
simulation of physical problems with several variables, for example to
generate error vectors to be added to a nominal vector. A particularly
common case is when the generated vector should be drawn from a
[Multivariate Normal Distribution](http://en.wikipedia.org/wiki/Multivariate_normal_distribution).

### Generating random vectors from a bivariate normal distribution

    // Create and seed a RandomGenerator (could use any of the generators in the random package here)
    RandomGenerator rg = new JDKRandomGenerator();
    rg.setSeed(17399225432l);  // Fixed seed means same results every time
    
    // Create a GassianRandomGenerator using rg as its source of randomness
    GaussianRandomGenerator rawGenerator = new GaussianRandomGenerator(rg);
    
    // Create a CorrelatedRandomVectorGenerator using rawGenerator for the components
    CorrelatedRandomVectorGenerator generator = 
        new CorrelatedRandomVectorGenerator(mean, covariance,
                                            1.0e-12 * covariance.getNorm(),
                                            rawGenerator);
    
    // Use the generator to generate correlated vectors
    double[] randomVector = generator.nextVector();
    ... 

The `mean` argument is a double[] array holding the means of the random vector
components.  In the bivariate case, it must have length 2.  The `covariance` argument
is a RealMatrix, which needs to be 2 x 2.  The main diagonal elements are the
variances of the vector components and the off-diagonal elements are the covariances.
For example, if the means are 1 and 2 respectively, and the desired standard deviations
are 3 and 4, respectively, then we need to use

    double[] mean = {1, 2};
    double[][] cov = {{9, c}, {c, 16}};
    RealMatrix covariance = MatrixUtils.createRealMatrix(cov);

where c is the desired covariance. If you are starting with a desired correlation,
you need to translate this to a covariance by multiplying it by the product of the
standard deviations.  For example, if you want to generate data that will give Pearson's
R of 0.5, you would use c = 3 * 4 * .5 = 6.

In addition to multivariate normal distributions, correlated vectors from multivariate uniform
distributions can be generated by creating a
[UniformRandomGenerator](../apidocs/org/hipparchus/random/UniformRandomGenerator.html)
in place of the `GaussianRandomGenerator` above.  More generally, any
[NormalizedRandomGenerator](../apidocs/org/hipparchus/random/NormalizedRandomGenerator.html)
may be used.


### Low discrepancy sequences

Several quasi-random sequences exist with the property that for all values of N, the subsequence
x<sub>1</sub>, ..., x<sub>N</sub> has low discrepancy, which results in equi-distributed samples.
While their quasi-randomness makes them unsuitable for most applications (i.e. the sequence of values
is completely deterministic), their unique properties give them an important advantage for quasi-Monte Carlo simulations.

Currently, the following low-discrepancy sequences are supported:

* generating random numbers
* generating random vectors
* generating random strings
* generating cryptographically secure sequences of random numbers or strings
* generating random samples and permutations
* analyzing distributions of values in an input file and generating values "like" the values in the file
* generating data for grouped frequency distributions or histograms


    // Create a Sobol sequence generator for 2-dimensional vectors
    RandomVectorGenerator generator = new SobolSequence(2);
    
    // Use the generator to generate vectors
    double[] randomVector = generator.nextVector();
    ... 

The figure below illustrates the unique properties of low-discrepancy sequences when generating N samples
in the interval [0, 1]. Roughly speaking, such sequences "fill" the respective space more evenly which leads to faster convergence in quasi-Monte Carlo simulations.

![Comparison of low-discrepancy sequences](images/userguide/low_discrepancy_sequences.png)

## Random Strings

The methods `nextHexString` and `nextSecureHexString`
can be used to generate random strings of hexadecimal characters.  Both
of these methods produce sequences of strings with good dispersion
properties.  The difference between the two methods is that the second can be
cryptographically secure (depending on the quality of the secure algorithm provider).  The implementation of
`nextHexString(n)` in `RandomDataGenerator` uses the
following simple algorithm to generate a string of `n` hex digits:

* `n/2+1` binary bytes are generated using the underlying RandomGenerator
* Each binary byte is translated into 2 hex digits
    
The `RandomDataGenerator` implementation of the "secure" version,
`nextSecureHexString` generates hex characters in 40-byte
"chunks" using a 3-step process:

* generating random numbers
* generating random vectors
* generating random strings
* generating cryptographically secure sequences of random numbers or strings
* generating random samples and permutations
* analyzing distributions of values in an input file and generating values "like" the values in the file
* generating data for grouped frequency distributions or histograms

Similarly to the secure random number generation methods,
`nextSecureHexString` is __much slower__ than
the non-secure version.  It should be used only for applications such as
generating unique session or transaction ids where predictability of
subsequent ids based on observation of previous values is a security
concern.  If all that is needed is an even distribution of hex characters
in the generated strings, the non-secure method should be used.

## Random permutations, combinations, sampling

To select a random sample of objects in a collection, you can use the
`nextSample` method implemented by `RandomDataGenerator`.
Specifically,  if `c` is a collection containing at least
`k` objects, and `randomData` is a `RandomDataGenerator` instance
`randomData.nextSample(c, k)` will return an `object[]` array of
length `k` consisting of elements randomly selected from the collection.
If `c` contains duplicate references, there may be duplicate
references in the returned array; otherwise returned elements will be
unique -- i.e., the sampling is without replacement among the object
references in the collection.

If `randomData` is a `RandomDataGenerator` instance, and
`n` and `k` are integers with `k <= n`, then
`randomData.nextPermutation(n, k)` returns an `int[]`
array of length `k` whose whose entries are selected randomly,
without repetition, from the integers `0` through
`n-1` (inclusive), i.e., `randomData.nextPermutation(n, k)` returns
a random permutation of  `n` taken `k` at a time.


## Generating data 'like' an input file

Using the `ValueServer` class, you can generate data based on
the values in an input file in one of two ways:


### Replay Mode

The following code will read data from `url` (a `java.net.URL` instance),
cycling through the values in the file in sequence, reopening and starting at
the beginning again when all values have been read.

          ValueServer vs = new ValueServer();
          vs.setValuesFileURL(url); 
          vs.setMode(ValueServer.REPLAY_MODE);
          vs.resetReplayFile();
          double value = vs.getNext();
          // ...Generate and use more values...
          vs.closeReplayFile();

The values in the file are not stored in memory, so it does not matter
how large the file is, but you do need to explicitly close the file
as above.  The expected file format is \n -delimited (i.e. one per line)
strings representing valid floating point numbers.


### Digest Mode

When used in Digest Mode, the ValueServer reads the entire input file
and estimates a probability density function based on data from the file.
The estimation method is essentially the
[Variable Kernel Method](http://nedwww.ipac.caltech.edu/level5/March02/Silverman/Silver2_6.html)
with Gaussian smoothing.  Once the density
has been estimated, `getNext()` returns random values whose
probability distribution matches the empirical distribution -- i.e., if
you generate a large number of such values, their distribution should
"look like" the distribution of the values in the input file.  The values
are not stored in memory in this case either, so there is no limit to the
size of the input file.  Here is an example:

          ValueServer vs = new ValueServer();
          vs.setValuesFileURL(url); 
          vs.setMode(ValueServer.DIGEST_MODE);
          vs.computeDistribution(500); //Read file and estimate distribution using 500 bins
          double value = vs.getNext();
          // ...Generate and use more values...

See the javadoc for `ValueServer` and `EmpiricalDistribution` for more details.
Note that `computeDistribution()` opens and closes the input file by itself.


## PRNG Pluggability

To enable alternative PRNGs to be "plugged in" to the Hipparchus data
generation utilities and to provide a generic means to replace
`java.util.Random` in applications, a random generator
adaptor framework has been added to Hipparchus.  The
[RandomGenerator](../apidocs/org/hipparchus/random/RandomGenerator.html)
interface abstracts the public interface of `java.util.Random` and any
implementation of this interface can be used as the source of random data
for the Hipparchus data generation classes.  An abstract base class,
[AbstractRandomGenerator](../apidocs/org/hipparchus/random/AbstractRandomGenerator.html)
is provided to make implementation easier. This class provides default
implementations of "derived" data generation methods based on the primitive,
`nextDouble()`. To support generic replacement of `java.util.Random`, the
[RandomAdaptor](../apidocs/org/hipparchus/random/RandomAdaptor.html)
class is provided, which extends `java.util.Random` and wraps and delegates calls to
a `RandomGenerator` instance.

Hipparchus provides by itself several implementations of the
[RandomGenerator](../apidocs/org.hipparchus/random/RandomGenerator.html) interface:

* generating random numbers
* generating random vectors
* generating random strings
* generating cryptographically secure sequences of random numbers or strings
* generating random samples and permutations
* analyzing distributions of values in an input file and generating values "like" the values in the file
* generating data for grouped frequency distributions or histograms

The JDK provided generator is a simple one that can be used only for very simple needs.
The Mersenne Twister is a fast generator with very good properties well suited for
Monte-Carlo simulation. It is equidistributed for generating vectors up to dimension 623
and has a huge period: 2<sup>19937</sup> - 1 (which is a Mersenne prime). This generator
is described in a paper by Makoto Matsumoto and Takuji Nishimura in 1998: <a
href="http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/ARTICLES/mt.pdf">Mersenne Twister:
A 623-Dimensionally Equidistributed Uniform Pseudo-Random Number Generator</a>, ACM
Transactions on Modeling and Computer Simulation, Vol. 8, No. 1, January 1998, pp 3--30.
The WELL generators are a family of generators with period ranging from 2<sup>512</sup> - 1
to 2<sup>44497</sup> - 1 (this last one is also a Mersenne prime) with even better properties
than Mersenne Twister. These generators are described in a paper by Fran&#231;ois Panneton,
Pierre L'Ecuyer and Makoto Matsumoto <a
href="http://www.iro.umontreal.ca/~lecuyer/myftp/papers/wellrng.pdf">Improved Long-Period
Generators Based on Linear Recurrences Modulo 2</a> ACM Transactions on Mathematical Software,
32, 1 (2006). The errata for the paper are in <a
href="http://www.iro.umontreal.ca/~lecuyer/myftp/papers/wellrng-errata.txt">wellrng-errata.txt</a>.


For simple sampling, any of these generators is sufficient. For Monte-Carlo simulations the
JDK generator does not have any of the good mathematical properties of the other generators,
so it should be avoided. The Mersenne twister and WELL generators have equidistribution properties
proven according to their bits pool size which is directly linked to their period (all of them
have maximal period, i.e. a generator with size n pool has a period 2<sup>n</sup>-1). They also
have equidistribution properties for 32 bits blocks up to s/32 dimension where s is their pool size.
So WELL19937c for exemple is equidistributed up to dimension 623 (19937/32). This means a Monte-Carlo
simulation generating a vector of n variables at each iteration has some guarantees on the properties
of the vector as long as its dimension does not exceed the limit. However, since we use bits from two
successive 32 bits generated integers to create one double, this limit is smaller when the variables are
of type double. so for Monte-Carlo simulation where less the 16 doubles are generated at each round,
WELL1024 may be sufficient. If a larger number of doubles are needed a generator with a larger pool
would be useful.


The WELL generators are more modern then MersenneTwister (the paper describing than has been published
in 2006 instead of 1998) and fix some of its (few) drawbacks. If initialization array contains many
zero bits, MersenneTwister may take a very long time (several hundreds of thousands of iterations to
reach a steady state with a balanced number of zero and one in its bits pool). So the WELL generators
are better to <i>escape zeroland</i> as explained by the WELL generators creators. The Well19937a and
Well44497a generator are not maximally equidistributed (i.e. there are some dimensions or bits blocks
size for which they are not equidistributed). The Well512a, Well1024a, Well19937c and Well44497b are
maximally equidistributed for blocks size up to 32 bits (they should behave correctly also for double
based on more than 32 bits blocks, but equidistribution is not proven at these blocks sizes).


The MersenneTwister generator uses a 624 elements integer array, so it consumes less than 2.5 kilobytes.
The WELL generators use 6 integer arrays with a size equal to the pool size, so for example the
WELL44497b generator uses about 33 kilobytes. This may be important if a very large number of
generator instances were used at the same time.


All generators are quite fast. As an example, here are some comparisons, obtained on a 64 bits JVM on a
linux computer with a 2008 processor (AMD phenom Quad 9550 at 2.2 GHz). The generation rate for
MersenneTwister was between 25 and 27 millions doubles per second (remember we generate two 32 bits integers for
each double). Generation rates for other PRNG, relative to MersenneTwister:


| <font size="+1">Example of performances</font> |
| --- |
| Name | generation rate (relative to MersenneTwister) |
| [MersenneTwister](../apidocs/org/hipparchus/random/MersenneTwister.html) | 1 |
| [JDKRandomGenerator](../apidocs/org/hipparchus/random/JDKRandomGenerator.html) | between 0.96 and 1.16 |
| [Well512a](../apidocs/org/hipparchus/random/Well512a.html) | between 0.85 and 0.88 |
| [Well1024a](../apidocs/org/hipparchus/random/Well1024a.html) | between 0.63 and 0.73 |
| [Well19937a](../apidocs/org/hipparchus/random/Well19937a.html) | between 0.70 and 0.71 |
| [Well19937c](../apidocs/org/hipparchus/random/Well19937c.html) | between 0.57 and 0.71 |
| [Well44497a](../apidocs/org/hipparchus/random/Well44497a.html) | between 0.69 and 0.71 |
| [Well44497b](../apidocs/org/hipparchus/random/Well44497b.html) | between 0.65 and 0.71 |


So for most simulation problems, the better generators like <a
href="../apidocs/org.hipparchus/random/Well19937c.html">Well19937c</a> and <a
href="../apidocs/org.hipparchus/random/Well44497b.html">Well44497b</a> are probably very good choices.

Note that *none* of these generators are suitable for cryptography. They are devoted
to simulation, and to generate very long series with strong properties on the series as a whole
(equidistribution, no correlation ...). They do not attempt to create small series but with
very strong properties of unpredictability as needed in cryptography.

Examples:

### Create a RandomGenerator based on RngPack's Mersenne Twister

To create a RandomGenerator using the RngPack Mersenne Twister PRNG
as the source of randomness, extend `AbstractRandomGenerator`
overriding the derived methods that the RngPack implementation provides:

    import edu.cornell.lassp.houle.RngPack.RanMT;

    /**
     * AbstractRandomGenerator based on RngPack RanMT generator.
     */
    public class RngPackGenerator extends AbstractRandomGenerator {
        
        private RanMT random = new RanMT();
        
        public void setSeed(long seed) {
           random = new RanMT(seed);
        }
        
        public double nextDouble() {
            return random.raw();
        }
        
        public double nextGaussian() {
            return random.gaussian();
        }
        
        public int nextInt(int n) {
            return random.choose(n);
        }
        
        public boolean nextBoolean() {
            return random.coin();
        }
    }


### Use the Mersenne Twister RandomGenerator in place of `java.util.Random` in `RandomData`

    RandomData randomData = new RandomDataImpl(new RngPackGenerator());

### Create an adaptor instance based on the Mersenne Twister generator that can be used in place of a `Random`

     RandomGenerator generator = new RngPackGenerator();
     Random random = RandomAdaptor.createAdaptor(generator);
     // random can now be used in place of a Random instance, data generation
     // calls will be delegated to the wrapped Mersenne Twister
