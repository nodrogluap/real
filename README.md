# real
Reproduction number Estimation Accounting for Lead time

A minimally parameterized algorithm for estimating on average how many cases stem from an existing case in a viral pandemic, for each date during the pandemic. Rather than fitting an exponential curve to the entire dataset, a weighted average of later-dated cases is used to derive an instantaneous [reproduction number](https://en.wikipedia.org/wiki/Basic_reproduction_number#Effective_reproduction_number). This weighting is [empiricially estimated](https://www.medrxiv.org/content/10.1101/2020.10.21.20217042v1) from the virus biology assuming that most cases are unobserved and asymptomatic (probably true for COVID-19 in many jurisdictions).

This script generates very similar Rt numbers to the R [EpiEstim](https://cran.r-project.org/web/packages/EpiEstim/index.html) package when the Serial Interval function is parameterized to be similar to the ```viral_shedding_proportions``` in this script. The main difference is that this script produces confidence intervals for Rt in a different manner, checking for Gaussian or exponential distribution of the Rt estimate in a rolling window. It also normalizes the data for frequent changes in the number of tests performed per date. This aligns more closely with observed SARS-CoV-2 testing results in a jurisdiction with high testing rates.

This script also greatly simplifies the process of plotting multiple Rt estimates (e.g. different "zones") together for comparison.

Test datasets for graphing are provided based on data scrapped from the [Alberta Government COVID-19 data explorer](https://www.alberta.ca/stats/covid-19-alberta-statistics.htm).
