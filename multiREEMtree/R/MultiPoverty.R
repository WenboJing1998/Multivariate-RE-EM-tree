#' Multi-dimensional Poverty Data Set (2012-2017)
#'
#' A data set containing indices of development of 109 countries during 2012 - 2017.
#'
#' @format A data frame with 541 rows and 10 variables:
#' \describe{
#'   \item{Access.to.drinking.water}{People using at least basic drinking water services (% of population).}
#'   \item{Access.to.electricity}{People with access to electricity (% of population).}
#'   \item{School.enrollment}{School enrollment, secondary (% gross) : the ratio of total enrollment,
#'   regardless of age, to the population of the age group that officially corresponds to the level of
#'   education shown.}
#'   \item{Survival.rate}{The estimated probability (%) of a newborn baby surviving to age five.}
#'   \item{Country}{Country.}
#'   \item{Year}{The year in which the sample was collected.}
#'   \item{Health.expenditure}{Current health expenditure per capita (current US$).}
#'   \item{Health.expenditure.percentage}{The percentage of health expenditure in GDP.}
#'   \item{Population.density}{Population density (people per sq. km of land area).}
#'   \item{Gini.index}{Gini index (World Bank estimate).}
#'   \item{Education.expenditure}{Current education expenditure per capita (current US$).}
#'   \item{Education.expenditure.percentage}{The percentage of Education Expenditure in GDP.}
#'   \item{GDP.per.capita}{Gross Domestic Product per capita.}
#'   \item{Rural.population}{Rural population (% of total population).}
#'   \item{Unemployment.rate}{Unemployment (% of total labor force) (modeled ILO estimate).}
#'   \item{Logistics.performance.index}{Quality of trade and transport-related infrastructure (1=low to 5=high).}
#'   \item{Corruption.perceptions.index}{An index indicating the perceived levels
#'   of public sector corruption.(0-100, lower index represents more corruption.)}
#'   \item{Democracy.index}{An index measuring the state of democracy.}
#'   \item{Production.index}{An index indicating the production quantity of agriculture.}
#'   \item{Temperature.change}{Mean surface temperature change relative to a baseline climatology in degrees Celsius.}
#'   \item{Deaths.in.violence.per.mille}{Deaths in violence conflict per mille.}
#'
#' }
#' @source \url{https://databank.worldbank.org/source/world-development-indicators},
#'         \url{https://tcdata360.worldbank.org/indicators/h345264a2?country=BRA&indicator=32534&viz=line_chart&years=2012,2020},
#'         \url{https://www.gapminder.org/data/documentation/democracy-index/},
#'         \url{http://www.fao.org/faostat/en/#data/QI},
#'         \url{http://www.fao.org/faostat/en/#data/ET},
#'         \url{https://datacatalog.worldbank.org/search/dataset/0041070}.
#'
#' @usage data(MultiPoverty)
#'
"MultiPoverty"
