# 1 Statistics
## 1.1 Overview
The statistics package provides frameworks and implementations for
basic Descriptive statistics, frequency distributions, bivariate regression,
and t-, chi-square and ANOVA test statistics.

[Descriptive statistics](#a1.2_Descriptive_statistics)<br/>
[Frequency distributions](#a1.3_Frequency_distributions)<br/>
[Simple Regression](#a1.4_Simple_regression)<br/>
[Multiple Regression](#a1.5_Multiple_linear_regression)<br/>
[Rank transformations](#a1.6_Rank_transformations)<br/>
[Covariance and correlation](#a1.7_Covariance_and_correlation)<br/>
[Statistical Tests](#a1.8_Statistical_tests)<br/>


## 1.2 Descriptive statistics
The stat package includes a framework and default implementations for
the following Descriptive statistics:
* arithmetic and geometric means
* variance and standard deviation
* sum, product, log sum, sum of squared values
* minimum, maximum, median, and percentiles
* skewness and kurtosis
* first, second, third and fourth moments


With the exception of percentiles and the median, all of these
statistics can be computed without maintaining the full list of input
data values in memory.  The stat package provides interfaces and
implementations that do not require value storage as well as
implementations that operate on arrays of stored values.

The top level interface is
[          UnivariateStatistic](../apidocs/org/hipparchus/stat/descriptive/UnivariateStatistic.html)
.
This interface, implemented by all statistics, consists of
`evaluate()` methods that take double[] arrays as arguments
and return the value of the statistic.   This interface is extended by
[          StorelessUnivariateStatistic](../apidocs/org/hipparchus/stat/descriptive/StorelessUnivariateStatistic.html)
, which adds `increment(),`
`getResult()` and associated methods to support
"storageless" implementations that maintain counters, sums or other
state information as values are added using the `increment()`
method.

Abstract implementations of the top level interfaces are provided in
[          AbstractUnivariateStatistic](../apidocs/org/hipparchus/stat/descriptive/AbstractUnivariateStatistic.html)
and
[          AbstractStorelessUnivariateStatistic](../apidocs/org/hipparchus/stat/descriptive/AbstractStorelessUnivariateStatistic.html)
respectively.

Each statistic is implemented as a separate class, in one of the
subpackages (moment, rank, summary) and each extends one of the abstract
classes above (depending on whether or not value storage is required to
compute the statistic). There are several ways to instantiate and use statistics.
Statistics can be instantiated and used directly,  but it is generally more convenient
(and efficient) to access them using the provided aggregates,
[           DescriptiveStatistics](../apidocs/org/hipparchus/stat/descriptive/DescriptiveStatistics.html)
and
[           SummaryStatistics.](../apidocs/org/hipparchus/stat/descriptive/SummaryStatistics.html)


`DescriptiveStatistics` maintains the input data in memory
and has the capability of producing "rolling" statistics computed from a
"window" consisting of the most recently added values.

`SummaryStatistics` does not store the input data values
in memory, so the statistics included in this aggregate are limited to those
that can be computed in one pass through the data without access to
the full array of values.

| Aggregate | Statistics Included | Values stored? | "Rolling" capability? |
| --- | --- | --- | --- |
| [ DescriptiveStatistics](../apidocs/org/hipparchus/stat/descriptive/DescriptiveStatistics.html) | min, max, mean, geometric mean, n, sum, sum of squares, standard deviation, variance, percentiles, skewness, kurtosis, median | Yes | Yes |
| [ SummaryStatistics](../apidocs/org/hipparchus/stat/descriptive/SummaryStatistics.html) | min, max, mean, geometric mean, n, sum, sum of squares, standard deviation, variance | No | No |

`SummaryStatistics` can be aggregated using
[          AggregateSummaryStatistics.](../apidocs/org/hipparchus/stat/descriptive/AggregateSummaryStatistics.html)
This class can be used to concurrently
gather statistics for multiple datasets as well as for a combined sample
including all of the data.

`MultivariateSummaryStatistics` is similar to
`SummaryStatistics` but handles n-tuple values instead of
scalar values. It can also compute the full covariance matrix for the
input data.

Neither `DescriptiveStatistics` nor `SummaryStatistics`
is thread-safe.
[           SynchronizedDescriptiveStatistics](../apidocs/org/hipparchus/stat/descriptive/SynchronizedDescriptiveStatistics.html)
and
[            SynchronizedSummaryStatistics](../apidocs/org/hipparchus/stat/descriptive/SynchronizedSummaryStatistics.html)
, respectively, provide thread-safe
versions for applications that require concurrent access to statistical
aggregates by multiple threads.
[            SynchronizedMultivariateSummaryStatistics](../apidocs/org/hipparchus/stat/descriptive/SynchronizedMultiVariateSummaryStatistics.html)
provides thread-safe
`MultivariateSummaryStatistics.`

There is also a utility class,
[          StatUtils](../apidocs/org/hipparchus/stat/StatUtils.html)
, that provides static methods for computing statistics
directly from double[] arrays.

Here are some examples showing how to compute Descriptive statistics.

### Compute summary statistics for a list of double values

Using the `DescriptiveStatistics` aggregate
(values are stored in memory):

    // Get a DescriptiveStatistics instance
    DescriptiveStatistics stats = new DescriptiveStatistics();
    
    // Add the data from the array
    for( int i = 0; i &lt; inputArray.length; i++) {
            stats.addValue(inputArray[i]);
    }
    
    // Compute some statistics
    double mean = stats.getMean();
    double std = stats.getStandardDeviation();
    double median = stats.getPercentile(50);

Using the `SummaryStatistics` aggregate (values are
<strong>not</strong> stored in memory):

    // Get a SummaryStatistics instance
    SummaryStatistics stats = new SummaryStatistics();
    
    // Read data from an input stream,
    // adding values and updating sums, counters, etc.
    while (line != null) {
            line = in.readLine();
            stats.addValue(Double.parseDouble(line.trim()));
    }
    in.close();
    
    // Compute the statistics
    double mean = stats.getMean();
    double std = stats.getStandardDeviation();
    //double median = stats.getMedian(); &lt;-- NOT AVAILABLE

Using the `StatUtils` utility class:

    // Compute statistics directly from the array
    // assume values is a double[] array
    double mean = StatUtils.mean(values);
    double std = FastMath.sqrt(StatUtils.variance(values));
    double median = StatUtils.percentile(values, 50);
    
    // Compute the mean of the first three values in the array
    mean = StatUtils.mean(values, 0, 3);

### Maintain a "rolling mean" of the most recent 100 values from an input stream

Use a `DescriptiveStatistics` instance with
window size set to 100

    // Create a DescriptiveStats instance and set the window size to 100
    DescriptiveStatistics stats = new DescriptiveStatistics();
    stats.setWindowSize(100);
    
    // Read data from an input stream,
    // displaying the mean of the most recent 100 observations
    // after every 100 observations
    long nLines = 0;
    while (line != null) {
            line = in.readLine();
            stats.addValue(Double.parseDouble(line.trim()));
            if (nLines == 100) {
                    nLines = 0;
                    System.out.println(stats.getMean());
           }
    }
    in.close();

### Compute statistics in a thread-safe manner

Use a `SynchronizedDescriptiveStatistics` instance

    // Create a SynchronizedDescriptiveStatistics instance and
    // use as any other DescriptiveStatistics instance
    DescriptiveStatistics stats = new SynchronizedDescriptiveStatistics();

### Compute statistics for multiple samples and overall statistics concurrently

There are two ways to do this using `AggregateSummaryStatistics.`
The first is to use an `AggregateSummaryStatistics` instance
to accumulate overall statistics contributed by `SummaryStatistics`
instances created using
[        AggregateSummaryStatistics.createContributingStatistics()](../apidocs/org/hipparchus/stat/descriptive/AggregateSummaryStatistics.html#createContributingStatistics())
:

    // Create a AggregateSummaryStatistics instance to accumulate the overall statistics 
    // and AggregatingSummaryStatistics for the subsamples
    AggregateSummaryStatistics aggregate = new AggregateSummaryStatistics();
    SummaryStatistics setOneStats = aggregate.createContributingStatistics();
    SummaryStatistics setTwoStats = aggregate.createContributingStatistics();
    // Add values to the subsample aggregates
    setOneStats.addValue(2);
    setOneStats.addValue(3);
    setTwoStats.addValue(2);
    setTwoStats.addValue(4);
    ...
    // Full sample data is reported by the aggregate
    double totalSampleSum = aggregate.getSum();
The above approach has the disadvantages that the `addValue` calls must be synchronized on the
`SummaryStatistics` instance maintained by the aggregate and each value addition updates the
aggregate as well as the subsample. For applications that can wait to do the aggregation until all values
have been added, a static
[          aggregate](../apidocs/org/hipparchus/stat/descriptive/AggregateSummaryStatistics.html#aggregate(java.util.Collection))
method is available, as shown in the following example.
This method should be used when aggregation needs to be done across threads.

    // Create SummaryStatistics instances for the subsample data
    SummaryStatistics setOneStats = new SummaryStatistics();
    SummaryStatistics setTwoStats = new SummaryStatistics();
    // Add values to the subsample SummaryStatistics instances
    setOneStats.addValue(2);
    setOneStats.addValue(3);
    setTwoStats.addValue(2);
    setTwoStats.addValue(4);
    ...
    // Aggregate the subsample statistics
    Collection&lt;SummaryStatistics&gt; aggregate = new ArrayList&lt;SummaryStatistics&gt;();
    aggregate.add(setOneStats);
    aggregate.add(setTwoStats);
    StatisticalSummary aggregatedStats = AggregateSummaryStatistics.aggregate(aggregate);
    
    // Full sample data is reported by aggregatedStats
    double totalSampleSum = aggregatedStats.getSum();


## 1.3 Frequency distributions
[          Frequency](../apidocs/org/hipparchus/stat/Frequency.html)

provides a simple interface for maintaining counts and percentages of discrete
values.

Strings, integers, longs and chars are all supported as value types,
as well as instances of any class that implements `Comparable.`
The ordering of values used in computing cumulative frequencies is by
default the <i>natural ordering,</i> but this can be overridden by supplying a
`Comparator` to the constructor. Adding values that are not
comparable to those that have already been added results in an
`IllegalArgumentException.`

Here are some examples.

###Compute a frequency distribution based on integer values

Mixing integers, longs, Integers and Longs:

     Frequency f = new Frequency();
     f.addValue(1);
     f.addValue(new Integer(1));
     f.addValue(new Long(1));
     f.addValue(2);
     f.addValue(new Integer(-1));
     System.out.prinltn(f.getCount(1));   // displays 3
     System.out.println(f.getCumPct(0));  // displays 0.2
     System.out.println(f.getPct(new Integer(1)));  // displays 0.6
     System.out.println(f.getCumPct(-2));   // displays 0
     System.out.println(f.getCumPct(10));  // displays 1

### Count string frequencies

Using case-sensitive comparison, alpha sort order (natural comparator):

    Frequency f = new Frequency();
    f.addValue("one");
    f.addValue("One");
    f.addValue("oNe");
    f.addValue("Z");
    System.out.println(f.getCount("one")); // displays 1
    System.out.println(f.getCumPct("Z"));  // displays 0.5
    System.out.println(f.getCumPct("Ot")); // displays 0.25

Using case-insensitive comparator:

    Frequency f = new Frequency(String.CASE_INSENSITIVE_ORDER);
    f.addValue("one");
    f.addValue("One");
    f.addValue("oNe");
    f.addValue("Z");
    System.out.println(f.getCount("one"));  // displays 3
    System.out.println(f.getCumPct("z"));  // displays 1



## 1.4 Simple regression
[          SimpleRegression](../apidocs/org/hipparchus/stat/regression/SimpleRegression.html)
provides ordinary least squares regression with
one independent variable estimating the linear model:

` y = intercept + slope * x  `

or

` y = slope * x `

Standard errors for `intercept` and `slope` are
available as well as ANOVA, r-square and Pearson's r statistics.

Observations (x,y pairs) can be added to the model one at a time or they
can be provided in a 2-dimensional array.  The observations are not stored
in memory, so there is no limit to the number of observations that can be
added to the model.

* arithmetic and geometric means
* variance and standard deviation
* sum, product, log sum, sum of squared values
* minimum, maximum, median, and percentiles
* skewness and kurtosis
* first, second, third and fourth moments


* arithmetic and geometric means
* variance and standard deviation
* sum, product, log sum, sum of squared values
* minimum, maximum, median, and percentiles
* skewness and kurtosis
* first, second, third and fourth moments


Here are some examples.

### Estimate a model based on observations added one at a time

Instantiate a regression instance and add data points

    regression = new SimpleRegression();
    regression.addData(1d, 2d);
    // At this point, with only one observation,
    // all regression statistics will return NaN
    
    regression.addData(3d, 3d);
    // With only two observations,
    // slope and intercept can be computed
    // but inference statistics will return NaN
    
    regression.addData(3d, 3d);
    // Now all statistics are defined.

Compute some statistics based on observations added so far

    System.out.println(regression.getIntercept());
    // displays intercept of regression line
    
    System.out.println(regression.getSlope());
    // displays slope of regression line
    
    System.out.println(regression.getSlopeStdErr());
    // displays slope standard error

Use the regression model to predict the y value for a new x value

    System.out.println(regression.predict(1.5d)
    // displays predicted y value for x = 1.5
More data points can be added and subsequent getXxx calls will incorporate
additional data in statistics.

### Estimate a model from a double[][] array of data points

Instantiate a regression object and load dataset

    double[][] data = { { 1, 3 }, {2, 5 }, {3, 7 }, {4, 14 }, {5, 11 }};
    SimpleRegression regression = new SimpleRegression();
    regression.addData(data);

Estimate regression model based on data

    System.out.println(regression.getIntercept());
    // displays intercept of regression line
    
    System.out.println(regression.getSlope());
    // displays slope of regression line
    
    System.out.println(regression.getSlopeStdErr());
    // displays slope standard error
More data points -- even another double[][] array -- can be added and subsequent
getXxx calls will incorporate additional data in statistics.

### Estimate a model from a double[][] array of data points, <em>excluding</em> the intercept

Instantiate a regression object and load dataset

    double[][] data = { { 1, 3 }, {2, 5 }, {3, 7 }, {4, 14 }, {5, 11 }};
    SimpleRegression regression = new SimpleRegression(false);
    //the argument, false, tells the class not to include a constant
    regression.addData(data);

Estimate regression model based on data

    System.out.println(regression.getIntercept());
    // displays intercept of regression line, since we have constrained the constant, 0.0 is returned
    
    System.out.println(regression.getSlope());
    // displays slope of regression line
    
    System.out.println(regression.getSlopeStdErr());
    // displays slope standard error
    
    System.out.println(regression.getInterceptStdErr() );
    // will return Double.NaN, since we constrained the parameter to zero
Caution must be exercised when interpreting the slope when no constant is being estimated. The slope
may be biased.



## 1.5 Multiple linear regression
[         OLSMultipleLinearRegression](../apidocs/org/hipparchus/stat/regression/OLSMultipleLinearRegression.html)
and
[         GLSMultipleLinearRegression](../apidocs/org/hipparchus/stat/regression/GLSMultipleLinearRegression.html)
provide least squares regression to fit the linear model:

` Y=X*b+u `

where Y is an n-vector <b>regressand</b>, X is a [n,k] matrix whose k columns are called
<b>regressors</b>, b is k-vector of <b>regression parameters</b> and u is an n-vector
of <b>error terms</b> or <b>residuals</b>.

[          OLSMultipleLinearRegression](../apidocs/org/hipparchus/stat/regression/OLSMultipleLinearRegression.html)
provides Ordinary Least Squares Regression, and
[          GLSMultipleLinearRegression](../apidocs/org/hipparchus/stat/regression/GLSMultipleLinearRegression.html)
implements Generalized Least Squares.  See the javadoc for these
classes for details on the algorithms and formulas used.

Data for OLS models can be loaded in a single double[] array, consisting of concatenated rows of data, each containing
the regressand (Y) value, followed by regressor values; or using a double[][] array with rows corresponding to
observations. GLS models also require a double[][] array representing the covariance matrix of the error terms.  See
[           AbstractMultipleLinearRegression#newSampleData(double[],int,int)](../apidocs/org/hipparchus/stat/regression/AbstractMultipleLinearRegression.html#newSampleData(double[], int, int))
,
[           OLSMultipleLinearRegression#newSampleData(double[], double[][])](../apidocs/org/hipparchus/stat/regression/OLSMultipleLinearRegression.html#newSampleData(double[], double[][]))
and
[           GLSMultipleLinearRegression#newSampleData(double[],double[][],double[][])](../apidocs/org/hipparchus/stat/regression/GLSMultipleLinearRegression.html#newSampleData(double[], double[][], double[][]))
for details.

* arithmetic and geometric means
* variance and standard deviation
* sum, product, log sum, sum of squared values
* minimum, maximum, median, and percentiles
* skewness and kurtosis
* first, second, third and fourth moments


Here are some examples.

###OLS regression

Instantiate an OLS regression object and load a dataset:

    OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();
    double[] y = new double[]{11.0, 12.0, 13.0, 14.0, 15.0, 16.0};
    double[][] x = new double[6][];
    x[0] = new double[]{0, 0, 0, 0, 0};
    x[1] = new double[]{2.0, 0, 0, 0, 0};
    x[2] = new double[]{0, 3.0, 0, 0, 0};
    x[3] = new double[]{0, 0, 4.0, 0, 0};
    x[4] = new double[]{0, 0, 0, 5.0, 0};
    x[5] = new double[]{0, 0, 0, 0, 6.0};          
    regression.newSampleData(y, x);

Get regression parameters and diagnostics:

    double[] beta = regression.estimateRegressionParameters();       
    
    double[] residuals = regression.estimateResiduals();
    
    double[][] parametersVariance = regression.estimateRegressionParametersVariance();
    
    double regressandVariance = regression.estimateRegressandVariance();
    
    double rSquared = regression.calculateRSquared();
    
    double sigma = regression.estimateRegressionStandardError();

### GLS regression

Instantiate a GLS regression object and load a dataset:

    GLSMultipleLinearRegression regression = new GLSMultipleLinearRegression();
    double[] y = new double[]{11.0, 12.0, 13.0, 14.0, 15.0, 16.0};
    double[][] x = new double[6][];
    x[0] = new double[]{0, 0, 0, 0, 0};
    x[1] = new double[]{2.0, 0, 0, 0, 0};
    x[2] = new double[]{0, 3.0, 0, 0, 0};
    x[3] = new double[]{0, 0, 4.0, 0, 0};
    x[4] = new double[]{0, 0, 0, 5.0, 0};
    x[5] = new double[]{0, 0, 0, 0, 6.0};          
    double[][] omega = new double[6][];
    omega[0] = new double[]{1.1, 0, 0, 0, 0, 0};
    omega[1] = new double[]{0, 2.2, 0, 0, 0, 0};
    omega[2] = new double[]{0, 0, 3.3, 0, 0, 0};
    omega[3] = new double[]{0, 0, 0, 4.4, 0, 0};
    omega[4] = new double[]{0, 0, 0, 0, 5.5, 0};
    omega[5] = new double[]{0, 0, 0, 0, 0, 6.6};
    regression.newSampleData(y, x, omega); 


## 1.6 Rank transformations
Some statistical algorithms require that input data be replaced by ranks.
The [         org.hipparchus.stat.ranking](../apidocs/org/hipparchus/stat/ranking/package-summary.html)
package provides rank transformation.
[         RankingAlgorithm](../apidocs/org/hipparchus/stat/ranking/RankingAlgorithm.html)
defines the interface for ranking.
[         NaturalRanking](../apidocs/org/hipparchus/stat/ranking/NaturalRanking.html)
provides an implementation that has two configuration options.
* arithmetic and geometric means
* variance and standard deviation
* sum, product, log sum, sum of squared values
* minimum, maximum, median, and percentiles
* skewness and kurtosis
* first, second, third and fourth moments


Examples:

    NaturalRanking ranking = new NaturalRanking(NaNStrategy.MINIMAL,
    TiesStrategy.MAXIMUM);
    double[] data = { 20, 17, 30, 42.3, 17, 50,
                      Double.NaN, Double.NEGATIVE_INFINITY, 17 };
    double[] ranks = ranking.rank(exampleData);
results in `ranks` containing `{6, 5, 7, 8, 5, 9, 2, 2, 5}.`

    new NaturalRanking(NaNStrategy.REMOVED,TiesStrategy.SEQUENTIAL).rank(exampleData);   
returns `{5, 2, 6, 7, 3, 8, 1, 4}.`

The default `NaNStrategy` is NaNStrategy.MAXIMAL.  This makes `NaN`
values larger than any other value (including `Double.POSITIVE_INFINITY`). The
default `TiesStrategy` is `TiesStrategy.AVERAGE,` which assigns tied
values the average of the ranks applicable to the sequence of ties.  See the
[        NaturalRanking](../apidocs/org/hipparchus/stat/ranking/NaturalRanking.html)
for more examples and [        TiesStrategy](../apidocs/org/hipparchus/stat/ranking/TiesStrategy.html)
and [NaNStrategy](../apidocs/org/hipparchus/stat/ranking/NaNStrategy.html)
for details on these configuration options.


## 1.7 Covariance and correlation
The [          org.hipparchus.stat.correlation](../apidocs/org/hipparchus/stat/correlation/package-summary.html)
package computes covariances
and correlations for pairs of arrays or columns of a matrix.
[          Covariance](../apidocs/org/hipparchus/stat/correlation/Covariance.html)
computes covariances,
[          PearsonsCorrelation](../apidocs/org/hipparchus/stat/correlation/PearsonsCorrelation.html)
provides Pearson's Product-Moment correlation coefficients,
[          SpearmansCorrelation](../apidocs/org/hipparchus/stat/correlation/SpearmansCorrelation.html)
computes Spearman's rank correlation and
[          KendallsCorrelation](../apidocs/org/hipparchus/stat/correlation/KendallsCorrelation.html)
computes Kendall's tau rank correlation.

<strong>Implementation Notes</strong>
* arithmetic and geometric means
* variance and standard deviation
* sum, product, log sum, sum of squared values
* minimum, maximum, median, and percentiles
* skewness and kurtosis
* first, second, third and fourth moments


<strong>Examples:</strong>
### <strong>Covariance of 2 arrays</strong>

To compute the unbiased covariance between 2 double arrays,
`x` and `y`, use:

    new Covariance().covariance(x, y)
For non-bias-corrected covariances, use

    covariance(x, y, false)

### <strong>Covariance matrix</strong>

A covariance matrix over the columns of a source matrix `data`
can be computed using

    new Covariance().computeCovarianceMatrix(data)
The i-jth entry of the returned matrix is the unbiased covariance of the ith and jth
columns of `data.` As above, to get non-bias-corrected covariances,
use

    computeCovarianceMatrix(data, false)

### <strong>Pearson's correlation of 2 arrays</strong><

To compute the Pearson's product-moment correlation between two double arrays
`x` and `y`, use:

    new PearsonsCorrelation().correlation(x, y)

### <strong>Pearson's correlation matrix</strong>

A (Pearson's) correlation matrix over the columns of a source matrix `data`
can be computed using

    new PearsonsCorrelation().computeCorrelationMatrix(data)
The i-jth entry of the returned matrix is the Pearson's product-moment correlation between the
ith and jth columns of `data.`

### <strong>Pearson's correlation significance and standard errors</strong>

To compute standard errors and/or significances of correlation coefficients
associated with Pearson's correlation coefficients, start by creating a
`PearsonsCorrelation` instance

    PearsonsCorrelation correlation = new PearsonsCorrelation(data);
where `data` is either a rectangular array or a `RealMatrix.`
Then the matrix of standard errors is

    correlation.getCorrelationStandardErrors();
The formula used to compute the standard error is <br/>
`SE<sub>r</sub> = ((1 - r<sup>2</sup>) / (n - 2))<sup>1/2</sup>`<br/>
where `r` is the estimated correlation coefficient and
`n` is the number of observations in the source dataset.<br/><br/>
<strong>p-values</strong> for the (2-sided) null hypotheses that elements of
a correlation matrix are zero populate the RealMatrix returned by

    correlation.getCorrelationPValues()
`getCorrelationPValues().getEntry(i,j)` is the
probability that a random variable distributed as `t<sub>n-2</sub>` takes
a value with absolute value greater than or equal to <br/>
`|r<sub>ij</sub>|((n - 2) / (1 - r<sub>ij</sub><sup>2</sup>))<sup>1/2</sup>`,
where `r<sub>ij</sub>` is the estimated correlation between the ith and jth
columns of the source array or RealMatrix. This is sometimes referred to as the
<i>significance</i> of the coefficient.<br/><br/>
For example, if `data` is a RealMatrix with 2 columns and 10 rows, then

    new PearsonsCorrelation(data).getCorrelationPValues().getEntry(0,1)
is the significance of the Pearson's correlation coefficient between the two columns
of `data`.  If this value is less than .01, we can say that the correlation
between the two columns of data is significant at the 99% level.

### <strong>Spearman's rank correlation coefficient</strong>

To compute the Spearman's rank-moment correlation between two double arrays
`x` and `y`:

    new SpearmansCorrelation().correlation(x, y)
This is equivalent to

    RankingAlgorithm ranking = new NaturalRanking();
    new PearsonsCorrelation().correlation(ranking.rank(x), ranking.rank(y))

### <strong>Kendalls's tau rank correlation coefficient</strong>

To compute the Kendall's tau rank correlation between two double arrays
`x` and `y`:

    new KendallsCorrelation().correlation(x, y)



## 1.8 Statistical tests
The [          org.hipparchus.stat.inference](../apidocs/org/hipparchus/stat/inference/)
package provides implementations for
[          Student's t](http://www.itl.nist.gov/div898/handbook/prc/section2/prc22.htm)
,
[          Chi-Square](http://www.itl.nist.gov/div898/handbook/eda/section3/eda35f.htm)
,
[G Test](http://en.wikipedia.org/wiki/G-test),
[          One-Way ANOVA](http://www.itl.nist.gov/div898/handbook/prc/section4/prc43.htm)
,
[          Mann-Whitney U](http://www.itl.nist.gov/div898/handbook/prc/section3/prc35.htm)
,
[          Wilcoxon signed rank](http://en.wikipedia.org/wiki/Wilcoxon_signed-rank_test)
and
[          Binomial](http://en.wikipedia.org/wiki/Binomial_test)
test statistics as well as
[          p-values](http://www.cas.lancs.ac.uk/glossary_v1.1/hyptest.html#pvalue)
associated with `t-`,
`Chi-Square`, `G`, `One-Way ANOVA`, `Mann-Whitney U`
`Wilcoxon signed rank`, and `Kolmogorov-Smirnov` tests.
The respective test classes are
[          TTest](../apidocs/org/hipparchus/stat/inference/TTest.html)
,
[          ChiSquareTest](../apidocs/org/hipparchus/stat/inference/ChiSquareTest.html)
,
[          GTest](../apidocs/org/hipparchus/stat/inference/GTest.html)
,
[          OneWayAnova](../apidocs/org/hipparchus/stat/inference/OneWayAnova.html)
,
[          MannWhitneyUTest](../apidocs/org/hipparchus/stat/inference/MannWhitneyUTest.html)
,
[          WilcoxonSignedRankTest](../apidocs/org/hipparchus/stat/inference/WilcoxonSignedRankTest.html)
,
[          BinomialTest](../apidocs/org/hipparchus/stat/inference/BinomialTest.html)
and
[          KolmogorovSmirnovTest](../apidocs/org/hipparchus/stat/inference/KolmogorovSmirnovTest.html)
.
The [          TestUtils](../apidocs/org/hipparchus/stat/inference/TestUtils.html)
class provides static methods to get test instances or
to compute test statistics directly.  The examples below all use the
static methods in `TestUtils` to execute tests.  To get
test object instances, either use e.g., `TestUtils.getTTest()`
or use the implementation constructors directly, e.g. `new TTest()`.

<strong>Implementation Notes</strong>
* arithmetic and geometric means
* variance and standard deviation
* sum, product, log sum, sum of squared values
* minimum, maximum, median, and percentiles
* skewness and kurtosis
* first, second, third and fourth moments


<strong>Examples:</strong>

### <strong>One-sample `t` tests</strong>

To compare the mean of a double[] array to a fixed value:

    double[] observed = {1d, 2d, 3d};
    double mu = 2.5d;
    System.out.println(TestUtils.t(mu, observed));
The code above will display the t-statistic associated with a one-sample
t-test comparing the mean of the `observed` values against
`mu.`

To compare the mean of a dataset described by a
[          StatisticalSummary](../apidocs/org/hipparchus/stat/descriptive/StatisticalSummary.html)
to a fixed value:

    double[] observed ={1d, 2d, 3d};
    double mu = 2.5d;
    SummaryStatistics sampleStats = new SummaryStatistics();
    for (int i = 0; i &lt; observed.length; i++) {
        sampleStats.addValue(observed[i]);
    }
    System.out.println(TestUtils.t(mu, observed));

To compute the p-value associated with the null hypothesis that the mean
of a set of values equals a point estimate, against the two-sided alternative that
the mean is different from the target value:

    double[] observed = {1d, 2d, 3d};
    double mu = 2.5d;
    System.out.println(TestUtils.tTest(mu, observed));
The snippet above will display the p-value associated with the null
hypothesis that the mean of the population from which the
`observed` values are drawn equals `mu.`

To perform the test using a fixed significance level, use:

    TestUtils.tTest(mu, observed, alpha);
where `0 &lt; alpha &lt; 0.5` is the significance level of
the test.  The boolean value returned will be `true` iff the
null hypothesis can be rejected with confidence `1 - alpha`.
To test, for example at the 95% level of confidence, use
`alpha = 0.05`

### <strong>Two-Sample t-tests</strong>

<strong>Example 1:</strong> Paired test evaluating
the null hypothesis that the mean difference between corresponding
(paired) elements of the `double[]` arrays
`sample1` and `sample2` is zero.
To compute the t-statistic:

    TestUtils.pairedT(sample1, sample2);

To compute the p-value:

    TestUtils.pairedTTest(sample1, sample2);

To perform a fixed significance level test with alpha = .05:

    TestUtils.pairedTTest(sample1, sample2, .05);

The last example will return `true` iff the p-value
returned by `TestUtils.pairedTTest(sample1, sample2)`
is less than `.05`

### <strong>Example 2: </strong> unpaired, two-sided, two-sample t-test using
`StatisticalSummary` instances, without assuming that
subpopulation variances are equal.

First create the `StatisticalSummary` instances.  Both
`DescriptiveStatistics` and `SummaryStatistics`
implement this interface.  Assume that `summary1` and
`summary2` are `SummaryStatistics` instances,
each of which has had at least 2 values added to the (virtual) dataset that
it describes.  The sample sizes do not have to be the same -- all that is required
is that both samples have at least 2 elements.

<strong>Note:</strong> The `SummaryStatistics` class does
not store the dataset that it describes in memory, but it does compute all
statistics necessary to perform t-tests, so this method can be used to
conduct t-tests with very large samples.  One-sample tests can also be
performed this way.
(See [Descriptive statistics](#1.2 Descriptive statistics) for details
on the `SummaryStatistics` class.)

To compute the t-statistic:

    TestUtils.t(summary1, summary2);

To compute the p-value:

    TestUtils.tTest(sample1, sample2);

To perform a fixed significance level test with alpha = .05:

    TestUtils.tTest(sample1, sample2, .05);

In each case above, the test does not assume that the subpopulation
variances are equal.  To perform the tests under this assumption,
replace "t" at the beginning of the method name with "homoscedasticT"

### <strong>Chi-square tests</strong>

To compute a chi-square statistic measuring the agreement between a
`long[]` array of observed counts and a `double[]`
array of expected counts, use:

    long[] observed = {10, 9, 11};
    double[] expected = {10.1, 9.8, 10.3};
    System.out.println(TestUtils.chiSquare(expected, observed));
the value displayed will be
`sum((expected[i] - observed[i])^2 / expected[i])`

To get the p-value associated with the null hypothesis that
`observed` conforms to `expected` use:

    TestUtils.chiSquareTest(expected, observed);

To test the null hypothesis that `observed` conforms to
`expected` with `alpha` significance level
(equiv. `100 * (1-alpha)%` confidence) where `
0 &lt; alpha &lt; 1 ` use:

    TestUtils.chiSquareTest(expected, observed, alpha);
The boolean value returned will be `true` iff the null hypothesis
can be rejected with confidence `1 - alpha`.

To compute a chi-square statistic statistic associated with a
[          chi-square test of independence](http://www.itl.nist.gov/div898/handbook/prc/section4/prc45.htm)
based on a two-dimensional (long[][])
`counts` array viewed as a two-way table, use:

    TestUtils.chiSquareTest(counts);
The rows of the 2-way table are
`count[0], ... , count[count.length - 1]. `<br/>
The chi-square statistic returned is
`sum((counts[i][j] - expected[i][j])^2/expected[i][j])`
where the sum is taken over all table entries and
`expected[i][j]` is the product of the row and column sums at
row `i`, column `j` divided by the total count.

To compute the p-value associated with the null hypothesis that
the classifications represented by the counts in the columns of the input 2-way
table are independent of the rows, use:

     TestUtils.chiSquareTest(counts);

To perform a chi-square test of independence with `alpha`
significance level (equiv. `100 * (1-alpha)%` confidence)
where `0 &lt; alpha &lt; 1 ` use:

    TestUtils.chiSquareTest(counts, alpha);
The boolean value returned will be `true` iff the null
hypothesis can be rejected with confidence `1 - alpha`.

### <strong>G tests</strong>

G tests are an alternative to chi-square tests that are recommended
when observed counts are small and / or incidence probabilities for
some cells are small. See Ted Dunning's paper,
[          Accurate Methods for the Statistics of Surprise and Coincidence](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.14.5962)
for
background and an empirical analysis showing now chi-square
statistics can be misleading in the presence of low incidence probabilities.
This paper also derives the formulas used in computing G statistics and the
root log likelihood ratio provided by the `GTest` class.

To compute a G-test statistic measuring the agreement between a
`long[]` array of observed counts and a `double[]`
array of expected counts, use:

    double[] expected = new double[]{0.54d, 0.40d, 0.05d, 0.01d};
    long[] observed = new long[]{70, 79, 3, 4};
    System.out.println(TestUtils.g(expected, observed));
the value displayed will be
`2 * sum(observed[i]) * log(observed[i]/expected[i])`

To get the p-value associated with the null hypothesis that
`observed` conforms to `expected` use:

    TestUtils.gTest(expected, observed);

To test the null hypothesis that `observed` conforms to
`expected` with `alpha` siginficance level
(equiv. `100 * (1-alpha)%` confidence) where `
0 &lt; alpha &lt; 1 ` use:

    TestUtils.gTest(expected, observed, alpha);
The boolean value returned will be `true` iff the null hypothesis
can be rejected with confidence `1 - alpha`.

To evaluate the hypothesis that two sets of counts come from the
same underlying distribution, use long[] arrays for the counts and
`gDataSetsComparison` for the test statistic

    long[] obs1 = new long[]{268, 199, 42};
    long[] obs2 = new long[]{807, 759, 184};
    System.out.println(TestUtils.gDataSetsComparison(obs1, obs2)); // G statistic
    System.out.println(TestUtils.gTestDataSetsComparison(obs1, obs2)); // p-value

For 2 x 2 designs, the `rootLogLikelihoodRatio` method
computes the
[          signed root log likelihood ratio.](http://tdunning.blogspot.com/2008/03/surprise-and-coincidence.html)
For example, suppose that for two events
A and B, the observed count of AB (both occurring) is 5, not A and B (B without A)
is 1995, A not B is 0; and neither A nor B is 10000.  Then

    new GTest().rootLogLikelihoodRatio(5, 1995, 0, 100000);
returns the root log likelihood associated with the null hypothesis that A
and B are independent.

### <strong>One-Way ANOVA tests</strong>

    double[] classA =
       {93.0, 103.0, 95.0, 101.0, 91.0, 105.0, 96.0, 94.0, 101.0 };
    double[] classB =
       {99.0, 92.0, 102.0, 100.0, 102.0, 89.0 };
    double[] classC =
       {110.0, 115.0, 111.0, 117.0, 128.0, 117.0 };
    List classes = new ArrayList();
    classes.add(classA);
    classes.add(classB);
    classes.add(classC);
Then you can compute ANOVA F- or p-values associated with the
null hypothesis that the class means are all the same
using a `OneWayAnova` instance or `TestUtils`
methods:

    double fStatistic = TestUtils.oneWayAnovaFValue(classes); // F-value
    double pValue = TestUtils.oneWayAnovaPValue(classes);     // P-value
To test perform a One-Way ANOVA test with significance level set at 0.01
(so the test will, assuming assumptions are met, reject the null
hypothesis incorrectly only about one in 100 times), use

    TestUtils.oneWayAnovaTest(classes, 0.01); // returns a boolean
                                              // true means reject null hypothesis

### <strong>Kolmogorov-Smirnov tests</strong>

Given a double[] array `data` of values, to evaluate the
null hypothesis that the values are drawn from a unit normal distribution

    final NormalDistribution unitNormal = new NormalDistribution(0d, 1d);
    TestUtils.kolmogorovSmirnovTest(unitNormal, sample, false)
returns the p-value and

    TestUtils.kolmogorovSmirnovStatistic(unitNormal, sample)
returns the D-statistic.

If `y` is a double array, to evaluate the null hypothesis that
`x` and `y` are drawn from the same underlying distribution,
use

    TestUtils.kolmogorovSmirnovStatistic(x, y)
to compute the D-statistic and

    TestUtils.kolmogorovSmirnovTest(x, y)
for the p-value associated with the null hypothesis that `x` and
`y` come from the same distribution. By default, here and above strict
inequality is used in the null hypothesis - i.e., we evaluate \(H_0 : D_{n,m} > d \).
To make the inequality above non-strict, add `false` as an actual parameter
above. For large samples, this parameter makes no difference.
<br/>
To force exact computation of the p-value (overriding the selection of estimation
method), first compute the d-statistic and then use the `exactP` method

    final double d = TestUtils.kolmogorovSmirnovStatistic(x, y);
    TestUtils.exactP(d, x.length, y.length, false)
assuming that the non-strict form of the null hypothesis is desired. Note, however,
that exact computation for large samples takes a long time.
