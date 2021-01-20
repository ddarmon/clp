#' Risk Factors Associated with Low Infant Birth Weight
#'
#' The `lbw` data frame has 189 rows and 6 columns. The data were collected at
#' Baystate Medical Center, Springfield, Mass during 1986. The original data
#' set is from the `MASS` package.
#'
#' @format A data frame with 189 rows and 6 columns:
#' \describe{
#'   \item{low}{factor for whether the birth weight is <= 2.5 kg or > 2.5 kg.}
#'   \item{age}{age of the mother in years.}
#'   \item{weight}{weight of mother in kg at last menstrual period.}
#'   \item{black}{indicator for mother's race being black.}
#'   \item{other}{indicator for mother's race being other.}
#'   \item{smoker}{indicator for mother's smoking status during pregnancy.}
#' }
#' @source Hosmer, D.W. and Lemeshow, S. (1989) *Applied Logistic Regression*. New York: Wiley
#' @references Venables, W. N. and Ripley, B. D. (2002) *Modern Applied Statistics with S*. Fourth edition. Springer.
"lbw"

#' Weight Loss Due to Healthy Low Carb versus Healthy Low Fat Diet
#'
#' The `dietstudy` data frame has 609 rows and 2 columns. The data are a
#' **synthetic** data set made to match the summary statistics from the study
#'
#' "Effect of Low-Fat vs Low-Carbohydrate Diet on 12-Month Weight Loss in
#' Overweight Adults and the Association With Genotype Pattern
#' or Insulin Secretion"
#'
#' from JAMA.
#'
#' The data consist of weight changes (negative indicates weight loss)
#' 12 months after direction to follow a Healthy Low Carb (HLC) or Healthy
#' Low Fat (HLF) diet.
#'
#' @format A data frame with 609 rows and 2 columns:
#' \describe{
#'   \item{weightchange}{weight change after 12 months in pounds.}
#'   \item{diet}{factor indicating prescribed diet of participant}
#' }
#' @source Christopher D. Gardner, John F. Trepanowski, Liana C. Del Gobbo, Michelle E. Hauser, Joseph Rigdon, John P. A. Ioannidis, Manisha Desai, and Abby C. King. "Effect of low-fat vs low-carbohydrate diet on 12-month weight loss in overweight adults and the association with genotype pattern or insulin secretion: the DIETFITS randomized clinical trial." *JAMA* 319, no. 7 (2018): 667-679.
"dietstudy"

#' Body Measurements for Predicting Body Fat Percentage in Males
#'
#' The `fat` data frame has 247 rows and 15 columns. These consist of physical
#' measurements of 247 men with ages ranging from 22 to 81.
#'
#' This data set can be used to demonstrate prediction of body fat
#' percentage (as estimated by water displacement) via easier body measurements
#' like height, weight, various circumferences, etc.
#'
#' @format A data frame with 247 rows and 15 columns:
#' \describe{
#'   \item{body.fat}{estimated body fat percentage using water displacement.}
#'   \item{age}{age in years.}
#'   \item{weight}{weight in pounds.}
#'   \item{height}{height in inches.}
#'   \item{BMI}{Body Mass Index.}
#'   \item{neck}{neck circumference in cm.}
#'   \item{chest}{chest circumference in cm.}
#'   \item{abdomen}{abdomen circumference in cm.}
#'   \item{hip}{hip circumference in cm.}
#'   \item{thigh}{thigh circumference in cm.}
#'   \item{knee}{knee circumference in cm.}
#'   \item{ankle}{ankle circumference in cm.}
#'   \item{bicep}{bicep circumference in cm.}
#'   \item{forearm}{forearm circumference in cm.}
#'   \item{wrist}{wrist circumference in cm.}
#' }
#' @source This data set comes from the collection of the Journal of Statistics
#' Education at \url{http://www.amstat.org/publications/jse/datasets/fat.txt}.
#' The data set was contributed by Roger W. Johnson.
#' @references The source of the data is attributed to Dr. A. Garth Fisher,
#' Human Performance Research Center, Brigham Young University,
#' Provo, Utah 84602,
"fat"
