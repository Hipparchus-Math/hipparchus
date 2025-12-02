#include <boost/math/distributions/inverse_gamma.hpp>
  using boost::math::inverse_gamma_distribution;

#include <iostream>
  using std::cout; using std::endl; using std::left; using std::showpoint; using std::noshowpoint;
#include <iomanip>
  using std::setw; using std::setprecision;
#include <limits>
  using std::numeric_limits;

int main()
{
  cout << "Inverse Gamma distribution with Boost" << endl;
  cout << "Using Boost "
       << BOOST_VERSION / 100000     << "."  // major version
       << BOOST_VERSION / 100 % 1000 << "."  // minor version
       << BOOST_VERSION % 100                // patch level
       << endl;

  try
  {
    int precision = 17;
    cout.precision(precision);
    // Construct an inverse Gamma distribution s
    double shape = 4.0;
    double scale = 2.0;
    inverse_gamma_distribution s = inverse_gamma_distribution(shape, scale);
    cout << "Inverse Gamma distribution with shape parameter " << shape << endl;
    cout << "Inverse Gamma distribution with scale parameter " << scale << endl;
    cout << "Mean = "<< mean(s) << ", variance = " << variance(s) << endl;
    cout << endl;
    double pset[] = {0.001, 0.01, 0.025, 0.05, 0.1, 0.999,0.990, 0.975, 0.950, 0.900};
    cout << "Percentile | Probability | Computed cumulative distribution " << endl;
    for (unsigned i = 0; i < sizeof(pset)/sizeof(pset[0]); ++i)
    {
      double x = quantile(s, pset[i]);
      cout << x << " "  << pset[i] << " "  << cdf(s, x) << endl;
    }
    cout << endl;


    double xset[] = {0.5, 1.0, 15.50};
    cout << "Evaluation points | Computed cumulative distribution " << endl;
    for (unsigned i = 0; i < sizeof(xset)/sizeof(xset[0]); ++i)
    {
      double x = xset[i];
      cout << x << " " << cdf(s, x) << endl;
    }
    cout << endl;

    double xpset[] = {0.1, 0.5, 1, 2, 5};
    cout << "Evaluation points | Computed probability density function" << endl;
    for (unsigned i = 0; i < sizeof(xpset)/sizeof(xpset[0]); ++i)
    {
      double x = xpset[i];
      cout << x << " "  << pdf(s, x) << endl;
    }
    cout << endl;

  }
  catch(const std::exception& e)
  {
    std::cout << "\n""Message from thrown exception was:\n   " << e.what() << std::endl;
  }
  return 0;
}

