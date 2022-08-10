#' @title Perform stepwise regression with verifying assumptions and identifying possible Box-Cox transformation
#' @name mlrpro-package
#' @aliases mlrpro
#' @description A tool for multiple regression, select independent variables,
#' check multiple linear regression assumptions and identify possible.
#' @usage mlrpro(Data,Y,Column_Y,Alpha)
#' @param Data a data frame containing the variables in the model.
#' @param Y the response variable.
#' @param Column_Y  the column response variable.
#' @param Alpha  significance level.
#' @return An object of class \code{mlrpro} is a list containing at least the following components:
#' \item{coefficients}{a named vector of coefficients.}
#' \item{residuals}{the residuals, that is response minus fitted values.}
#' \item{fitted.values}{the fitted mean values.}
#' \item{rank}{the numeric rank of the fitted linear model.}
#' \item{df.residual}{the residual degrees of freedom.}
#' \item{call}{the matched call.}
#' \item{terms}{the terms object used.}
#' \item{model}{if requested (the default), the model frame used.}
#' \item{lambda}{lambda value utilized in the data conversion.}
#' @importFrom dplyr mutate %>%
#' @importFrom MASS boxcox
#' @importFrom car leveneTest
#' @importFrom stats filter lm median shapiro.test step
#' @examples
#' data(trees)
#' Model1 <- mlrpro(Data = trees,Y = trees$Volume, Column_Y = 3, Alpha = 0.05)
#' ## or ##
#' data(mtcars)
#' Model2 <- mlrpro(Data = mtcars,Y = mtcars$mpg, Column_Y = 1 , Alpha = 0.01)
#' @export

#install.packages("dplyr");library(dplyr)
#install.packages("car");library(car)
#install.packages("MASS");library(MASS)


mlrpro <- function(Data,Y,Column_Y,Alpha) {
  Newdata <- Data
  Newdata[Column_Y] <- NULL
  y <- Y
  Newdata <- data.frame(y,Newdata)
  fit <- suppressWarnings(step(lm(y~.,data=Newdata ,direct="both"),trace = 0))
  sumfit <- suppressWarnings(summary(fit))
  Number_Beta <- as.numeric(nrow(sumfit$coefficients))

  Decision <- 0
  i <- 1
  while (i <= Number_Beta) {
    x <- ifelse(sumfit$coefficients[i,4]<=Alpha,"Sig","NoSig")
    Decision[i] <- x
    i = i+1
  }

  view <- data.frame(sumfit$coefficients)
  view <- mutate(view, Decision )
  NoSig <- subset(view,view$Decision == "NoSig")
  RowNoSig <- as.numeric(nrow(NoSig))

  Sig <-  subset(view,view$Decision == "Sig")
  RowSig <- as.numeric(nrow(Sig))

  if (RowSig <= 0) {
    Newdata1 <- Newdata
    Find.Lambda <- (boxcox(y~.,data= Newdata1))
    Find.Lambda2 <- data.frame(Find.Lambda$x,Find.Lambda$y)
    Find.Lambda3 <- subset(Find.Lambda2,Find.Lambda$y==max(Find.Lambda$y))
    Optimal.lambda <- Find.Lambda3$Find.Lambda.x
    Optimal.lambda.true <- Find.Lambda3$Find.Lambda.x
    ifelse (Optimal.lambda >= 0.75 && Optimal.lambda < 1.5,Optimal.lambda <- 1,
            Optimal.lambda <- Optimal.lambda)

    if (Optimal.lambda == 1 ) {

      writeLines(c("-------------------------",
                   "Stepwise regression model ",
                   "-------------------------"))

      print(suppressWarnings(summary(fit)))

      writeLines(c("----------------------",
                   "Checking  Assumptions",
                   "---------------------"))

      writeLines("[1] The errors do not follow a normal distribution.")
      writeLines("[1] The variance of the errors is not constant (Heteroscedastic).")


      writeLines("------------------------\n Box cox transformation \n------------------------")
      writeLines("[1] Optimal lambda approximate to 1.")

      message("The multiple linear regression may not be appropriate for this data.")

      plot(fit,1)
      plot(fit,2)
      list(
        coefficients = fit$coefficients,
        residuals = fit$residuals,
        fitted.values = fit$fitted.values,
        rank = fit$rank,
        df.residual = fit$df.residual,
        call = fit$call,
        terms = fit$terms,
        model = fit$model
      )

    }
    else if (Optimal.lambda != 1) {

      writeLines(c("--------------------------------------------------",
                   "Regression model derived from data transformations",
                   "--------------------------------------------------"))


      if (Optimal.lambda > -0.25 && Optimal.lambda < 0.25){
        lambda <- 0

        writeLines(c("The lambda value utilized in the data conversion is 0 ",
                     "and transformation the dependent variable in the following : y=log(y)"))

        y.prime <- log(y)
        Newdata1$y <- NULL
        Newdata1 <- data.frame(y.prime,Newdata1)
        fit_end <- suppressWarnings(lm(y.prime~.,data=Newdata1))
        sum_fit_end <- suppressWarnings(summary(fit_end))

      }
      else if (Optimal.lambda >= 0.25 && Optimal.lambda < 0.75 ) {
        lambda <- 0.5

        writeLines(c("The lambda value utilized in the data conversion is 0.5",
                     "and transformation the dependent variable in the following : y=sqrt(y)"))

        y.prime <- sqrt(y)
        Newdata1$y <- NULL
        Newdata1 <- data.frame(y.prime,Newdata1)
        fit_end <- suppressWarnings(lm(y.prime~., data=Newdata1))
        sum_fit_end <- suppressWarnings(summary(fit_end))

      }
      else if (Optimal.lambda >= 1.5 && Optimal.lambda <= 2 ) {
        lambda <- 2

        writeLines(c("The lambda value utilized in the data conversion is 2 ",
                     "and transformation the dependent variable in the following : y=y^2"))

        y.prime <- (y^2)
        Newdata1$y <- NULL
        Newdata1 <- data.frame(y.prime,Newdata1)
        fit_end <- suppressWarnings(lm(y.prime~.,data=Newdata1))
        sum_fit_end <- suppressWarnings(summary(fit_end))

      }
      else if (Optimal.lambda <= -0.25 && Optimal.lambda > -0.75 ) {
        lambda <- -0.5

        writeLines(c("The lambda value utilized in the data conversion is -0.5 ",
                     " and transformation the dependent variable in the following : y=1/sqrt(y)"))

        y.prime <- 1/sqrt(y)
        Newdata1$y <- NULL
        Newdata1 <- data.frame(y.prime,Newdata1)
        fit_end <- suppressWarnings(lm(y.prime~.,data=Newdata1))
        sum_fit_end <- suppressWarnings(summary(fit_end))

      }
      else if (Optimal.lambda <= -0.75 && Optimal.lambda > -1.5 ) {
        lambda <- -1

        writeLines(c("The lambda value utilized in the data conversion is -1 ",
                     "and transformation the dependent variable in the following : y=1/y"))

        y.prime <- 1/y
        Newdata1$y <- NULL
        Newdata1 <- data.frame(y.prime,Newdata1)
        fit_end <- suppressWarnings(lm(y.prime~.,data=Newdata1))
        sum_fit_end <- suppressWarnings(summary(fit_end))

      }
      else if (Optimal.lambda <= -1.5 && Optimal.lambda >= -2 ) {
        lambda <- -2

        writeLines(c("The lambda value utilized in the data conversion is ",
                     " and transformation the dependent variable in the following : y=1/y^2"))

        y.prime <- 1/(y^2)
        Newdata1$y <- NULL
        Newdata1 <- data.frame(y.prime,Newdata1)
        fit_end <- suppressWarnings(lm(y.prime~.,data=Newdata1))
        sum_fit_end <- suppressWarnings(summary(fit_end))

      }

      error <- fit_end$residuals
      error.Group <- factor(error<=median(error))
      Normal <- shapiro.test(error)
      variance <- leveneTest(error,group = error.Group)
      variance_p <- variance$`Pr(>F)`[1]

      plot(fit_end,1)
      plot(fit_end,2)

      print(suppressWarnings(summary(fit_end)))

      writeLines(c ("----------------------",
                    "Checking  Assumptions",
                    "---------------------"))


      ifelse(Normal$p.value <= Alpha,
             print.noquote (" The errors do not follow a normal distribution."),
             print.noquote (" The errors follow a normal distribution."))

      ifelse( variance_p <= Alpha,
              print.noquote(" The variance of the errors is not constant (Heteroscedastic)."),
              print.noquote(" The variance of the errors is constant (Homoscedastic)."))

      list(
        coefficients = fit_end$coefficients,
        residuals = fit_end$residuals,
        fitted.values = fit_end$fitted.values,
        rank = fit_end$rank,
        df.residual = fit_end$df.residual,
        call = fit_end$call,
        terms = fit_end$terms,
        model = fit_end$model,
        lambda = lambda
      )

    }

  }
  else {
    beta0 <- sumfit$coefficients[1,4]
    RowNoSig <- as.numeric(nrow(NoSig))

    if (beta0 <= Alpha && RowNoSig > 0)  {
      delete.Intercept <- row.names(Sig)
      delete.Intercept <- delete.Intercept[-1]
      Newdata1 <- Newdata[delete.Intercept]

      while(RowNoSig > 0) {
        Newdata1 <- Newdata[delete.Intercept]
        Newdata1 <- data.frame(y,Newdata1)
        fit1 <- suppressWarnings(lm(y~., data=Newdata1))
        sumfit1 <- suppressWarnings(summary(fit1))
        Number_Beta1 <- as.numeric(nrow(sumfit1$coefficients))

        Decision1 <- 0
        i <- 1
        while( i <= Number_Beta1 ) {
          x <- ifelse(sumfit1$coefficients[i,4]<= Alpha,"Sig","NoSig")
          Decision1[i] <- x
          i = i+1
        }

        view1 <- data.frame(sumfit1$coefficients)
        view1 <- mutate(view1, Decision1)
        Sig <- subset(view1,view1$Decision == "Sig")

        delete.Intercept <- row.names(Sig)
        delete.Intercept <- delete.Intercept[-1]
        NoSig <- subset(view1,view1$Decision == "NoSig")
        RowNoSig <- as.numeric(nrow(NoSig))
        RowNoSig = RowNoSig+0
      }

      error <- fit1$residuals
      error.Group <- factor(error<=median(error))
      Normal <- shapiro.test(error)
      variance <- leveneTest(error,group = error.Group)
      variance_p <- variance$`Pr(>F)`[1]

      if (Normal$p.value <= Alpha || variance_p <=Alpha) {

        Find.Lambda <- (boxcox(y~.,data= Newdata1))
        Find.Lambda2 <- data.frame(Find.Lambda$x,Find.Lambda$y)
        Find.Lambda3 <- subset(Find.Lambda2,Find.Lambda$y==max(Find.Lambda$y))
        Optimal.lambda <- Find.Lambda3$Find.Lambda.x
        Optimal.lambda.true <- Find.Lambda3$Find.Lambda.x

        ifelse (Optimal.lambda >= 0.75 && Optimal.lambda < 1.5,Optimal.lambda <- 1,
                Optimal.lambda <- Optimal.lambda )

        if (Optimal.lambda == 1 ) {

          writeLines("---------------------------\n Stepwise regression model \n---------------------------")


          print(suppressWarnings(summary(fit1)))

          writeLines ("----------------------\n Checking Assumptions \n----------------------")

          writeLines("[1] The errors do not follow a normal distribution.")
          writeLines("[1] The variance of the errors is not constant (Heteroscedastic).")

          writeLines("------------------------\n Box cox transformation \n------------------------")
          writeLines("Optimal lambda approximate to 1.")

          message("The multiple linear regression may not be appropriate for this data.")

          plot(fit1,1)
          plot(fit1,2)
          list(
            coefficients = fit1$coefficients,
            residuals = fit1$residuals,
            fitted.values = fit1$fitted.values,
            rank = fit1$rank,
            df.residual = fit1$df.residual,
            call = fit1$call,
            terms = fit1$terms,
            model = fit1$model
          )

        } else if (Optimal.lambda != 1) {


          writeLines(c("--------------------------------------------------",
                       "Regression model derived from data transformations",
                       "--------------------------------------------------"))

          if (Optimal.lambda > -0.25 && Optimal.lambda < 0.25){
            lambda <- 0

            writeLines(c("The lambda value utilized in the data conversion is 0",
                         "and transformation the dependent variable in the following : y = log(y)"))

            y.prime <- log(y)
            Newdata1$y <- NULL
            Newdata1 <- data.frame(y.prime,Newdata1)
            fit_end <- suppressWarnings(lm(y.prime~.,data=Newdata1))
            sum_fit_end <- suppressWarnings(summary(fit_end))

          }
          else if (Optimal.lambda >= 0.25 && Optimal.lambda < 0.75 ) {
            lambda <- 0.5

            writeLines(c("The lambda value utilized in the data conversion is 0.5 ",
                         "and transformation the dependent variable in the following : y=sqrt(y)"))

            y.prime <- sqrt(y)
            Newdata1$y <- NULL
            Newdata1 <- data.frame(y.prime,Newdata1)
            fit_end <- suppressWarnings(lm(y.prime~.,data=Newdata1))
            sum_fit_end <- suppressWarnings(summary(fit_end))

          }
          else if (Optimal.lambda >= 1.5 && Optimal.lambda <= 2 ) {
            lambda <- 2

            writeLines(c("The lambda value utilized in the data conversion is 2 ",
                         "and transformation the dependent variable in the following : y=y^2"))

            y.prime <- (y^2)
            Newdata1$y <- NULL
            Newdata1 <- data.frame(y.prime,Newdata1)
            fit_end <- suppressWarnings(lm(y.prime~.,data=Newdata1))
            sum_fit_end <- suppressWarnings(summary(fit_end))

          }
          else if (Optimal.lambda <= -0.25 && Optimal.lambda > -0.75 ) {
            lambda <- -0.5

            writeLines(c("The lambda value utilized in the data conversion is -0.5 ",
                         "and transformation the dependent variable in the following : y=1/sqrt(y)"))

            y.prime <- 1/sqrt(y)
            Newdata1$y <- NULL
            Newdata1 <- data.frame(y.prime,Newdata1)
            fit_end <- suppressWarnings(lm(y.prime~.,data=Newdata1))
            sum_fit_end <- suppressWarnings(summary(fit_end))

          }
          else if (Optimal.lambda <= -0.75 && Optimal.lambda > -1.5 ) {
            lambda <- -1

            writeLines(c("The lambda value utilized in the data conversion is -1 ",
                         "and transformation the dependent variable in the following : y=1/y"))

            y.prime <- 1/y
            Newdata1$y <- NULL
            Newdata1 <- data.frame(y.prime,Newdata1)
            fit_end <- suppressWarnings(lm(y.prime~.,data=Newdata1))
            sum_fit_end <- suppressWarnings(summary(fit_end))

          }
          else if (Optimal.lambda <= -1.5 && Optimal.lambda>= -2 ) {
            lambda <- -2

            writeLines(c("The lambda value utilized in the data conversion is -2",
                         "and transformation the dependent variable in the following : y=1/y^2"))

            y.prime <- 1/(y^2)
            Newdata1$y <- NULL
            Newdata1 <- data.frame(y.prime,Newdata1)
            fit_end <- suppressWarnings(lm(y.prime~.,data=Newdata1))
            sum_fit_end <- suppressWarnings(summary(fit_end))

          }

          Number_Beta <- as.numeric(nrow(sum_fit_end$coefficients))
          Decision <- 0
          i <- 1
          while (i <= Number_Beta) {
            x <- ifelse(sum_fit_end$coefficients[i,4]<=Alpha,"Sig","NoSig")
            Decision[i] <- x
            i = i+1
          }

          view <- data.frame(sum_fit_end$coefficients)
          view <- mutate(view, Decision )
          NoSig <- subset(view,view$Decision == "NoSig")
          RowNoSig <- as.numeric(nrow(NoSig))

          Sig <-  subset(view,view$Decision == "Sig")

          beta0 <- sum_fit_end$coefficients[1,4]

          if (beta0 <= Alpha && RowNoSig > 0) {
            delete.Intercept <- row.names(Sig)
            delete.Intercept <- delete.Intercept[-1]
            Newdata1 <- Newdata[delete.Intercept]

            while(RowNoSig > 0) {
              Newdata1 <- Newdata[delete.Intercept]
              Newdata1 <- data.frame(y,Newdata1)
              fit1 <- suppressWarnings(lm(y~., data=Newdata1))
              sumfit1 <- suppressWarnings(summary(fit1))
              Number_Beta1 <- as.numeric(nrow(sumfit1$coefficients))

              Decision1 <- 0
              i <- 1
              while( i <= Number_Beta1 ) {
                x <- ifelse(sumfit1$coefficients[i,4]<= Alpha,"Sig","NoSig")
                Decision1[i] <- x
                i = i+1
              }

              view1 <- data.frame(sumfit1$coefficients)
              view1 <- mutate(view1, Decision1)
              Sig <- subset(view1,view1$Decision == "Sig")

              delete.Intercept <- row.names(Sig)
              delete.Intercept <- delete.Intercept[-1]
              NoSig <- subset(view1,view1$Decision == "NoSig")
              RowNoSig <- as.numeric(nrow(NoSig))
              RowNoSig = RowNoSig+0
            }

            error <- fit1$residuals
            error.Group <- factor(error<=median(error))
            Normal <- shapiro.test(error)
            variance <- leveneTest(error,group = error.Group)
            variance_p <- variance$`Pr(>F)`[1]

            plot(fit1,1)
            plot(fit1,2)

            print(suppressWarnings(summary(fit1)))

            writeLines(c ("----------------------",
                          "Checking  Assumptions",
                          "---------------------"))


            ifelse(Normal$p.value <= Alpha,
                   print.noquote (" The errors do not follow a normal distribution."),
                   print.noquote (" The errors follow a normal distribution."))

            ifelse( variance_p <= Alpha,
                    print.noquote(" The variance of the errors is not constant (Heteroscedastic)."),
                    print.noquote(" The variance of the errors is constant (Homoscedastic)."))

            list(
              coefficients = fit1$coefficients,
              residuals = fit1$residuals,
              fitted.values = fit1$fitted.values,
              rank = fit1$rank,
              df.residual = fit1$df.residual,
              call = fit1$call,
              terms = fit1$terms,
              model = fit1$model,
              lambda = lambda
            )
          }
          else if (beta0 <= Alpha && RowNoSig <= 0) {
            error_end <- fit_end$residuals
            error.Group <- factor(error_end<=median(error_end))
            Normal <- shapiro.test(error_end)
            variance <- leveneTest(error_end,group = error.Group)
            variance_p <- variance$`Pr(>F)`[1]
            plot(fit_end,1)
            plot(fit_end,2)

            print(suppressWarnings(summary(fit_end)))

            writeLines(c ("----------------------",
                          "Checking  Assumptions",
                          "---------------------"))



            ifelse(Normal$p.value <= Alpha,
                   print.noquote (" The errors do not follow a normal distribution."),
                   print.noquote (" The errors follow a normal distribution."))

            ifelse( variance_p <= Alpha,
                    print.noquote(" The variance of the errors is not constant (Heteroscedastic)."),
                    print.noquote(" The variance of the errors is constant (Homoscedastic)."))

            list(
              coefficients = fit_end$coefficients,
              residuals = fit_end$residuals,
              fitted.values = fit_end$fitted.values,
              rank = fit_end$rank,
              df.residual = fit_end$df.residual,
              call = fit_end$call,
              terms = fit_end$terms,
              model = fit_end$model,
              lambda = lambda
            )
          }
          else if (beta0 > Alpha && RowNoSig > 0) {
            delete.Intercept <- row.names(Sig)
            Newdata1 <- Newdata[delete.Intercept]

            while(RowNoSig > 0) {
              Newdata1 <- Newdata[delete.Intercept]
              Newdata1 <- data.frame(y,Newdata1)
              fit1 <- suppressWarnings(lm(y~0+.,data=Newdata1))
              sumfit1 <- suppressWarnings(summary(fit1))
              Number_Beta1 <- as.numeric(nrow(sumfit1$coefficients))

              Decision1 <- 0
              i <- 1
              while( i <= Number_Beta1 ) {
                x <- ifelse(sumfit1$coefficients[i,4]<= Alpha,"Sig","NoSig")
                Decision1[i] <- x
                i = i+1
              }

              view1 <- data.frame(sumfit1$coefficients)
              view1 <- mutate(view1, Decision1)
              Sig <- subset(view1,view1$Decision == "Sig")

              delete.Intercept <- row.names(Sig)
              NoSig <- subset(view1,view1$Decision == "NoSig")
              RowNoSig <- as.numeric(nrow(NoSig))
              RowNoSig = RowNoSig+0
            }

            error <- fit1$residuals
            error.Group <- factor(error<=median(error))
            Normal <- shapiro.test(error)
            variance <- leveneTest(error,group = error.Group)
            variance_p <- variance$`Pr(>F)`[1]

            plot(fit1,1)
            plot(fit1,2)

            print(suppressWarnings(summary(fit1)))

            writeLines(c ("----------------------",
                          "Checking  Assumptions",
                          "---------------------"))


            ifelse(Normal$p.value <= Alpha,
                   print.noquote (" The errors do not follow a normal distribution."),
                   print.noquote (" The errors follow a normal distribution."))

            ifelse( variance_p <= Alpha,
                    print.noquote(" The variance of the errors is not constant (Heteroscedastic)."),
                    print.noquote(" The variance of the errors is constant (Homoscedastic)."))

            list(
              coefficients = fit1$coefficients,
              residuals = fit1$residuals,
              fitted.values = fit1$fitted.values,
              rank = fit1$rank,
              df.residual = fit1$df.residual,
              call = fit1$call,
              terms = fit1$terms,
              model = fit1$model,
              lambda = lambda
            )


          }
          else {
            delete.Intercept <- row.names(Sig)
            Newdata1 <- Newdata[delete.Intercept]
            Newdata1 <- data.frame(y,Newdata1)
            fit1 <- suppressWarnings(lm(y~0+.,data=Newdata1))
            sumfit1 <- suppressWarnings(summary(fit1))

            error <- sumfit1$residuals
            error.Group <- factor(error<=median(error))
            Normal <- shapiro.test(error)
            variance <- leveneTest(error,group = error.Group)
            variance_p <- variance$`Pr(>F)`[1]

            plot(fit1,1)
            plot(fit1,2)

            print(suppressWarnings(summary(fit1)))

            writeLines(c ("----------------------",
                          "Checking  Assumptions",
                          "---------------------"))



            ifelse(Normal$p.value <= Alpha,
                   print.noquote (" The errors do not follow a normal distribution."),
                   print.noquote (" The errors follow a normal distribution."))

            ifelse( variance_p <= Alpha,
                    print.noquote(" The variance of the errors is not constant (Heteroscedastic)."),
                    print.noquote(" The variance of the errors is constant (Homoscedastic)."))

            list(
              coefficients = fit1$coefficients,
              residuals = fit1$residuals,
              fitted.values = fit1$fitted.values,
              rank = fit1$rank,
              df.residual = fit1$df.residual,
              call = fit1$call,
              terms = fit1$terms,
              model = fit1$model,
              lambda = lambda
            )


          }

        }




      }
      else if (Normal$p.value >= Alpha && variance_p >= Alpha) {

        writeLines(c("-------------------------",
                     "Stepwise regression model ",
                     "-------------------------"))

        print(suppressWarnings(summary(fit1)))
        plot(fit1,1)
        plot(fit1,2)

        writeLines(c("----------------------",
                     "Checking  Assumptions",
                     "---------------------"))

        writeLines("[1] The errors follow a normal distribution.")
        writeLines("[1] The variance of the errors is constant (Homoscedastic).")

        list(
          coefficients = fit1$coefficients,
          residuals = fit1$residuals,
          fitted.values = fit1$fitted.values,
          rank = fit1$rank,
          df.residual = fit1$df.residual,
          call = fit1$call,
          terms = fit1$terms,
          model = fit1$model

        )
      }

    }
    else if (beta0 <= Alpha && RowNoSig <= 0) {
      delete.Intercept <- row.names(Sig)
      delete.Intercept <- delete.Intercept[-1]
      Newdata1 <- Newdata[delete.Intercept]

      error <- sumfit$residuals
      error.Group <- factor(error<=median(error))
      Normal <- shapiro.test(error)
      variance <- leveneTest(error,group = error.Group)
      variance_p <- variance$`Pr(>F)`[1]

      if (Normal$p.value <= Alpha || variance_p <= Alpha) {
        Newdata1 <- data.frame(y,Newdata1)
        Find.Lambda <- (boxcox(y~., data=Newdata1))
        Find.Lambda2 <- data.frame(Find.Lambda$x,Find.Lambda$y)
        Find.Lambda3 <- subset(Find.Lambda2,Find.Lambda$y==max(Find.Lambda$y))
        Optimal.lambda <- Find.Lambda3$Find.Lambda.x
        Optimal.lambda.true <- Find.Lambda3$Find.Lambda.x

        ifelse (Optimal.lambda >= 0.75 && Optimal.lambda < 1.5,Optimal.lambda <- 1,
                Optimal.lambda <- Optimal.lambda)

        if (Optimal.lambda == 1 ) {

          writeLines(c("-------------------------",
                       "Stepwise regression model ",
                       "-------------------------"))

          print(suppressWarnings(summary(fit)))

          writeLines(c("----------------------",
                       "Checking  Assumptions",
                       "---------------------"))

          writeLines("[1] The errors do not follow a normal distribution.")
          writeLines("[1] The variance of the errors is not constant (Heteroscedastic).")


          writeLines("------------------------\n Box cox transformation \n------------------------")
          writeLines("[1] Optimal lambda approximate to 1.")

          message("The multiple linear regression may not be appropriate for this data.")

          plot(fit,1)
          plot(fit,2)
          list(
            coefficients = fit$coefficients,
            residuals = fit$residuals,
            fitted.values = fit$fitted.values,
            rank = fit$rank,
            df.residual = fit$df.residual,
            call = fit$call,
            terms = fit$terms,
            model = fit$model
          )

        }
        else if (Optimal.lambda != 1) {

          writeLines(c("--------------------------------------------------",
                       "Regression model derived from data transformations",
                       "--------------------------------------------------"))


          if (Optimal.lambda > -0.25 && Optimal.lambda < 0.25){
            lambda <- 0

            writeLines(c("The lambda value utilized in the data conversion is 0 ",
                         "and transformation the dependent variable in the following : y=log(y)"))

            y.prime <- log(y)
            Newdata1$y <- NULL
            Newdata1 <- data.frame(y.prime,Newdata1)
            fit_end <- suppressWarnings(lm(y.prime~.,data=Newdata1))
            sum_fit_end <- suppressWarnings(summary(fit_end))

          }
          else if (Optimal.lambda >= 0.25 && Optimal.lambda < 0.75 ) {
            lambda <- 0.5

            writeLines(c("The lambda value utilized in the data conversion is 0.5",
                         "and transformation the dependent variable in the following : y=sqrt(y)"))

            y.prime <- sqrt(y)
            Newdata1$y <- NULL
            Newdata1 <- data.frame(y.prime,Newdata1)
            fit_end <- suppressWarnings(lm(y.prime~., data=Newdata1))
            sum_fit_end <- suppressWarnings(summary(fit_end))

          }
          else if (Optimal.lambda >= 1.5 && Optimal.lambda <= 2 ) {
            lambda <- 2

            writeLines(c("The lambda value utilized in the data conversion is 2 ",
                         "and transformation the dependent variable in the following : y=y^2"))

            y.prime <- (y^2)
            Newdata1$y <- NULL
            Newdata1 <- data.frame(y.prime,Newdata1)
            fit_end <- suppressWarnings(lm(y.prime~.,data=Newdata1))
            sum_fit_end <- suppressWarnings(summary(fit_end))

          }
          else if (Optimal.lambda <= -0.25 && Optimal.lambda > -0.75 ) {
            lambda <- -0.5

            writeLines(c("The lambda value utilized in the data conversion is -0.5 ",
                         " and transformation the dependent variable in the following : y=1/sqrt(y)"))

            y.prime <- 1/sqrt(y)
            Newdata1$y <- NULL
            Newdata1 <- data.frame(y.prime,Newdata1)
            fit_end <- suppressWarnings(lm(y.prime~.,data=Newdata1))
            sum_fit_end <- suppressWarnings(summary(fit_end))

          }
          else if (Optimal.lambda <= -0.75 && Optimal.lambda > -1.5 ) {
            lambda <- -1

            writeLines(c("The lambda value utilized in the data conversion is -1 ",
                         "and transformation the dependent variable in the following : y=1/y"))

            y.prime <- 1/y
            Newdata1$y <- NULL
            Newdata1 <- data.frame(y.prime,Newdata1)
            fit_end <- suppressWarnings(lm(y.prime~.,data=Newdata1))
            sum_fit_end <- suppressWarnings(summary(fit_end))

          }
          else if (Optimal.lambda <= -1.5 && Optimal.lambda >= -2 ) {
            lambda <- -2

            writeLines(c("The lambda value utilized in the data conversion is ",
                         " and transformation the dependent variable in the following : y=1/y^2"))

            y.prime <- 1/(y^2)
            Newdata1$y <- NULL
            Newdata1 <- data.frame(y.prime,Newdata1)
            fit_end <- suppressWarnings(lm(y.prime~.,data=Newdata1))
            sum_fit_end <- suppressWarnings(summary(fit_end))

          }

          Number_Beta <- as.numeric(nrow(sum_fit_end$coefficients))
          Decision <- 0
          i <- 1
          while (i <= Number_Beta) {
            x <- ifelse(sum_fit_end$coefficients[i,4]<=Alpha,"Sig","NoSig")
            Decision[i] <- x
            i = i+1
          }

          view <- data.frame(sum_fit_end$coefficients)
          view <- mutate(view, Decision )
          NoSig <- subset(view,view$Decision == "NoSig")
          RowNoSig <- as.numeric(nrow(NoSig))

          Sig <-  subset(view,view$Decision == "Sig")

          beta0 <- sum_fit_end$coefficients[1,4]

          if (beta0 <= Alpha && RowNoSig > 0) {
            delete.Intercept <- row.names(Sig)
            delete.Intercept <- delete.Intercept[-1]
            Newdata1 <- Newdata[delete.Intercept]

            while(RowNoSig > 0) {
              Newdata1 <- Newdata[delete.Intercept]
              Newdata1 <- data.frame(y,Newdata1)
              fit1 <- suppressWarnings(lm(y~., data=Newdata1))
              sumfit1 <- suppressWarnings(summary(fit1))
              Number_Beta1 <- as.numeric(nrow(sumfit1$coefficients))

              Decision1 <- 0
              i <- 1
              while( i <= Number_Beta1 ) {
                x <- ifelse(sumfit1$coefficients[i,4]<= Alpha,"Sig","NoSig")
                Decision1[i] <- x
                i = i+1
              }

              view1 <- data.frame(sumfit1$coefficients)
              view1 <- mutate(view1, Decision1)
              Sig <- subset(view1,view1$Decision == "Sig")

              delete.Intercept <- row.names(Sig)
              delete.Intercept <- delete.Intercept[-1]
              NoSig <- subset(view1,view1$Decision == "NoSig")
              RowNoSig <- as.numeric(nrow(NoSig))
              RowNoSig = RowNoSig+0
            }

            error <- fit1$residuals
            error.Group <- factor(error<=median(error))
            Normal <- shapiro.test(error)
            variance <- leveneTest(error,group = error.Group)
            variance_p <- variance$`Pr(>F)`[1]

            plot(fit1,1)
            plot(fit1,2)

            print(suppressWarnings(summary(fit1)))

            writeLines(c ("----------------------",
                          "Checking  Assumptions",
                          "---------------------"))


            ifelse(Normal$p.value <= Alpha,
                   print.noquote (" The errors do not follow a normal distribution."),
                   print.noquote (" The errors follow a normal distribution."))

            ifelse( variance_p <= Alpha,
                    print.noquote(" The variance of the errors is not constant (Heteroscedastic)."),
                    print.noquote(" The variance of the errors is constant (Homoscedastic)."))

            list(
              coefficients = fit1$coefficients,
              residuals = fit1$residuals,
              fitted.values = fit1$fitted.values,
              rank = fit1$rank,
              df.residual = fit1$df.residual,
              call = fit1$call,
              terms = fit1$terms,
              model = fit1$model,
              lambda = lambda
            )
          }
          else if (beta0 <= Alpha && RowNoSig <= 0) {
            error_end <- fit_end$residuals
            error.Group <- factor(error_end<=median(error_end))
            Normal <- shapiro.test(error_end)
            variance <- leveneTest(error_end,group = error.Group)
            variance_p <- variance$`Pr(>F)`[1]
            plot(fit_end,1)
            plot(fit_end,2)

            print(suppressWarnings(summary(fit_end)))

            writeLines(c ("----------------------",
                          "Checking  Assumptions",
                          "---------------------"))



            ifelse(Normal$p.value <= Alpha,
                   print.noquote (" The errors do not follow a normal distribution."),
                   print.noquote (" The errors follow a normal distribution."))

            ifelse( variance_p <= Alpha,
                    print.noquote(" The variance of the errors is not constant (Heteroscedastic)."),
                    print.noquote(" The variance of the errors is constant (Homoscedastic)."))

            list(
              coefficients = fit_end$coefficients,
              residuals = fit_end$residuals,
              fitted.values = fit_end$fitted.values,
              rank = fit_end$rank,
              df.residual = fit_end$df.residual,
              call = fit_end$call,
              terms = fit_end$terms,
              model = fit_end$model,
              lambda = lambda
            )
          }
          else if (beta0 > Alpha && RowNoSig > 0) {
            delete.Intercept <- row.names(Sig)
            Newdata1 <- Newdata[delete.Intercept]

            while(RowNoSig > 0) {
              Newdata1 <- Newdata[delete.Intercept]
              Newdata1 <- data.frame(y,Newdata1)
              fit1 <- suppressWarnings(lm(y~0+.,data=Newdata1))
              sumfit1 <- suppressWarnings(summary(fit1))
              Number_Beta1 <- as.numeric(nrow(sumfit1$coefficients))

              Decision1 <- 0
              i <- 1
              while( i <= Number_Beta1 ) {
                x <- ifelse(sumfit1$coefficients[i,4]<= Alpha,"Sig","NoSig")
                Decision1[i] <- x
                i = i+1
              }

              view1 <- data.frame(sumfit1$coefficients)
              view1 <- mutate(view1, Decision1)
              Sig <- subset(view1,view1$Decision == "Sig")

              delete.Intercept <- row.names(Sig)
              NoSig <- subset(view1,view1$Decision == "NoSig")
              RowNoSig <- as.numeric(nrow(NoSig))
              RowNoSig = RowNoSig+0
            }

            error <- fit1$residuals
            error.Group <- factor(error<=median(error))
            Normal <- shapiro.test(error)
            variance <- leveneTest(error,group = error.Group)
            variance_p <- variance$`Pr(>F)`[1]

            plot(fit1,1)
            plot(fit1,2)

            print(suppressWarnings(summary(fit1)))

            writeLines(c ("----------------------",
                          "Checking  Assumptions",
                          "---------------------"))


            ifelse(Normal$p.value <= Alpha,
                   print.noquote (" The errors do not follow a normal distribution."),
                   print.noquote (" The errors follow a normal distribution."))

            ifelse( variance_p <= Alpha,
                    print.noquote(" The variance of the errors is not constant (Heteroscedastic)."),
                    print.noquote(" The variance of the errors is constant (Homoscedastic)."))

            list(
              coefficients = fit1$coefficients,
              residuals = fit1$residuals,
              fitted.values = fit1$fitted.values,
              rank = fit1$rank,
              df.residual = fit1$df.residual,
              call = fit1$call,
              terms = fit1$terms,
              model = fit1$model,
              lambda = lambda
            )


          }
          else {
            delete.Intercept <- row.names(Sig)
            Newdata1 <- Newdata[delete.Intercept]
            Newdata1 <- data.frame(y,Newdata1)
            fit1 <- suppressWarnings(lm(y~0+.,data=Newdata1))
            sumfit1 <- suppressWarnings(summary(fit1))

            error <- sumfit1$residuals
            error.Group <- factor(error<=median(error))
            Normal <- shapiro.test(error)
            variance <- leveneTest(error,group = error.Group)
            variance_p <- variance$`Pr(>F)`[1]

            plot(fit1,1)
            plot(fit1,2)

            print(suppressWarnings(summary(fit1)))

            writeLines(c ("----------------------",
                          "Checking  Assumptions",
                          "---------------------"))



            ifelse(Normal$p.value <= Alpha,
                   print.noquote (" The errors do not follow a normal distribution."),
                   print.noquote (" The errors follow a normal distribution."))

            ifelse( variance_p <= Alpha,
                    print.noquote(" The variance of the errors is not constant (Heteroscedastic)."),
                    print.noquote(" The variance of the errors is constant (Homoscedastic)."))

            list(
              coefficients = fit1$coefficients,
              residuals = fit1$residuals,
              fitted.values = fit1$fitted.values,
              rank = fit1$rank,
              df.residual = fit1$df.residual,
              call = fit1$call,
              terms = fit1$terms,
              model = fit1$model,
              lambda = lambda
            )


          }


        }




      }
      else if (Normal$p.value >= Alpha && variance_p >= Alpha) {
        plot(fit,1)
        plot(fit,2)

        writeLines(c("-------------------------",
                     "Stepwise regression model ",
                     "-------------------------"))

        print(suppressWarnings(summary(fit)))


        writeLines(c ("----------------------",
                      "Checking  Assumptions",
                      "---------------------"))

        writeLines("[1] The errors follow a normal distribution.")
        writeLines("[1] The variance of the errors is constant (Homoscedastic).")


        list(
          coefficients = fit$coefficients,
          residuals = fit$residuals,
          fitted.values = fit$fitted.values,
          rank = fit$rank,
          df.residual = fit$df.residual,
          call = fit$call,
          terms = fit$terms,
          model = fit$model
        )
      }

    }
    else if (beta0 > Alpha && RowNoSig > 0) {
      delete.Intercept <- row.names(Sig)
      Newdata1 <- Newdata[delete.Intercept]

      while(RowNoSig > 0) {
        Newdata1 <- Newdata[delete.Intercept]
        Newdata1 <- data.frame(y,Newdata1)
        fit1 <- suppressWarnings(lm(y~0+.,data=Newdata1))
        sumfit1 <- suppressWarnings(summary(fit1))
        Number_Beta1 <- as.numeric(nrow(sumfit1$coefficients))

        Decision1 <- 0
        i <- 1
        while( i <= Number_Beta1 ) {
          x <- ifelse(sumfit1$coefficients[i,4]<= Alpha,"Sig","NoSig")
          Decision1[i] <- x
          i = i+1
        }

        view1 <- data.frame(sumfit1$coefficients)
        view1 <- mutate(view1, Decision1)
        Sig <- subset(view1,view1$Decision == "Sig")

        delete.Intercept <- row.names(Sig)
        NoSig <- subset(view1,view1$Decision == "NoSig")
        RowNoSig <- as.numeric(nrow(NoSig))
        RowNoSig = RowNoSig+0
      }

      error <- fit1$residuals
      error.Group <- factor(error<=median(error))
      Normal <- shapiro.test(error)
      variance <- leveneTest(error,group = error.Group)
      variance_p <- variance$`Pr(>F)`[1]

      if (Normal$p.value <= Alpha || variance_p <= Alpha) {

        Find.Lambda <- (boxcox(y~.,data=Newdata1))
        Find.Lambda2 <- data.frame(Find.Lambda$x,Find.Lambda$y)
        Find.Lambda3 <- subset(Find.Lambda2,Find.Lambda$y==max(Find.Lambda$y))
        Optimal.lambda <- Find.Lambda3$Find.Lambda.x
        Optimal.lambda.true  <- Find.Lambda3$Find.Lambda.x

        ifelse (Optimal.lambda >= 0.75 && Optimal.lambda < 1.5,Optimal.lambda <- 1,
                Optimal.lambda <- Optimal.lambda)

        if (Optimal.lambda == 1 ) {

          writeLines(c("-------------------------",
                       "Stepwise regression model ",
                       "-------------------------"))

          print(suppressWarnings(summary(fit1)))

          writeLines(c("----------------------",
                       "Checking  Assumptions",
                       "---------------------"))

          writeLines("[1] The errors do not follow a normal distribution.")
          writeLines("[1] The variance of the errors is not constant (Heteroscedastic).")

          writeLines("------------------------\n Box cox transformation \n------------------------")
          writeLines("Optimal lambda approximate to 1.")

          message("The multiple linear regression may not be appropriate for this data.")
          plot(fit1,1)
          plot(fit1,2)
          list(
            coefficients = fit1$coefficients,
            residuals = fit1$residuals,
            fitted.values = fit1$fitted.values,
            rank = fit1$rank,
            df.residual = fit1$df.residual,
            call = fit1$call,
            terms = fit1$terms,
            model = fit1$model
          )

        } else if (Optimal.lambda!= 1) {

          writeLines(c ("---------------------------------------------------",
                        "Regression model derived from data transformations",
                        "---------------------------------------------------"))

          if (Optimal.lambda> -0.25 && Optimal.lambda< 0.25){
            lambda <- 0

            writeLines(c("The lambda value utilized in the data conversion is 0",
                         "and transformation the dependent variable in the following : y = log(y)"))

            y.prime <- log(y)
            Newdata1$y <- NULL
            Newdata1 <- data.frame(y.prime,Newdata1)
            fit_end <- suppressWarnings(lm(y.prime~0+.,data=Newdata1))
            sum_fit_end <- suppressWarnings(summary(fit_end))

          }
          else if (Optimal.lambda >= 0.25 && Optimal.lambda < 0.75 ) {
            lambda <- 0.5

            writeLines(c("The lambda value utilized in the data conversion is 0.5 ",
                         "and transformation the dependent variable in the following : y=sqrt(y)"))

            y.prime <- sqrt(y)
            Newdata1$y <- NULL
            Newdata1 <- data.frame(y.prime,Newdata1)
            fit_end <- suppressWarnings(lm(y.prime~0+.,data= Newdata1))
            sum_fit_end <- suppressWarnings(summary(fit_end))


          }
          else if (Optimal.lambda >= 1.5 && Optimal.lambda <= 2 ) {
            lambda <- 2

            writeLines(c("The lambda value utilized in the data conversion is 2",
                         "and transformation the dependent variable in the following : y=y^2"))

            y.prime <- (y^2)
            Newdata1$y <- NULL
            Newdata1 <- data.frame(y.prime,Newdata1)
            fit_end <- suppressWarnings(lm(y.prime~0+.,data=Newdata1))
            sum_fit_end <- suppressWarnings(summary(fit_end))

          }
          else if (Optimal.lambda <= -0.25 && Optimal.lambda > -0.75 ) {
            lambda <- -0.5

            writeLines(c("The lambda value utilized in the data conversion is -0.5",
                         "and transformation the dependent variable in the following : y=1/sqrt(y)"))

            y.prime <- 1/sqrt(y)
            Newdata1$y <- NULL
            Newdata1 <- data.frame(y.prime,Newdata1)
            fit_end <- suppressWarnings(lm(y.prime~0+.,data=Newdata1))
            sum_fit_end <- suppressWarnings(summary(fit_end))

          }
          else if (Optimal.lambda <= -0.75 && Optimal.lambda > -1.5 ) {
            lambda <- -1

            writeLines(c("The lambda value utilized in the data conversion is -1",
                         "and transformation the dependent variable in the following : y=1/y"))

            y.prime <- 1/y
            Newdata1$y <- NULL
            Newdata1 <- data.frame(y.prime,Newdata1)
            fit_end <- suppressWarnings(lm(y.prime~0+.,data =  Newdata1))
            sum_fit_end <- suppressWarnings(summary(fit_end))

          }
          else if (Optimal.lambda <= -1.5 && Optimal.lambda >= -2 ) {
            lambda <- -2

            writeLines(c("\n","The lambda value utilized in the data conversion is -2 ",
                         " and transformation the dependent variable in the following : y=1/y^2"))

            y.prime <- 1/(y^2)
            Newdata1$y <- NULL
            Newdata1 <- data.frame(y.prime,Newdata1)
            fit_end <- suppressWarnings(lm(y.prime~0+.,data =  Newdata1))
            sum_fit_end <- suppressWarnings(summary(fit_end))

          }

          Number_Beta <- as.numeric(nrow(sum_fit_end$coefficients))
          Decision <- 0
          i <- 1
          while (i <= Number_Beta) {
            x <- ifelse(sum_fit_end$coefficients[i,4]<=Alpha,"Sig","NoSig")
            Decision[i] <- x
            i = i+1
          }

          view <- data.frame(sum_fit_end$coefficients)
          view <- mutate(view, Decision )
          NoSig <- subset(view,view$Decision == "NoSig")
          RowNoSig <- as.numeric(nrow(NoSig))

          Sig <-  subset(view,view$Decision == "Sig")

          beta0 <- sum_fit_end$coefficients[1,4]

          if (beta0 <= Alpha && RowNoSig > 0) {
            delete.Intercept <- row.names(Sig)
            Newdata1 <- Newdata[delete.Intercept]

            while(RowNoSig > 0) {
              Newdata1 <- Newdata[delete.Intercept]
              Newdata1 <- data.frame(y,Newdata1)
              fit1 <- suppressWarnings(lm(y~0+., data=Newdata1))
              sumfit1 <- suppressWarnings(summary(fit1))
              Number_Beta1 <- as.numeric(nrow(sumfit1$coefficients))

              Decision1 <- 0
              i <- 1
              while( i <= Number_Beta1 ) {
                x <- ifelse(sumfit1$coefficients[i,4]<= Alpha,"Sig","NoSig")
                Decision1[i] <- x
                i = i+1
              }

              view1 <- data.frame(sumfit1$coefficients)
              view1 <- mutate(view1, Decision1)
              Sig <- subset(view1,view1$Decision == "Sig")

              delete.Intercept <- row.names(Sig)
              delete.Intercept <- delete.Intercept[-1]
              NoSig <- subset(view1,view1$Decision == "NoSig")
              RowNoSig <- as.numeric(nrow(NoSig))
              RowNoSig = RowNoSig+0
            }

            error <- fit1$residuals
            error.Group <- factor(error<=median(error))
            Normal <- shapiro.test(error)
            variance <- leveneTest(error,group = error.Group)
            variance_p <- variance$`Pr(>F)`[1]

            plot(fit1,1)
            plot(fit1,2)

            print(suppressWarnings(summary(fit1)))

            writeLines(c ("----------------------",
                          "Checking  Assumptions",
                          "---------------------"))


            ifelse(Normal$p.value <= Alpha,
                   print.noquote (" The errors do not follow a normal distribution."),
                   print.noquote (" The errors follow a normal distribution."))

            ifelse( variance_p <= Alpha,
                    print.noquote(" The variance of the errors is not constant (Heteroscedastic)."),
                    print.noquote(" The variance of the errors is constant (Homoscedastic)."))

            list(
              coefficients = fit1$coefficients,
              residuals = fit1$residuals,
              fitted.values = fit1$fitted.values,
              rank = fit1$rank,
              df.residual = fit1$df.residual,
              call = fit1$call,
              terms = fit1$terms,
              model = fit1$model,
              lambda = lambda
            )
          }
          else if (beta0 <= Alpha && RowNoSig <= 0) {
            error_end <- fit_end$residuals
            error.Group <- factor(error_end<=median(error_end))
            Normal <- shapiro.test(error_end)
            variance <- leveneTest(error_end,group = error.Group)
            variance_p <- variance$`Pr(>F)`[1]
            plot(fit_end,1)
            plot(fit_end,2)

            print(suppressWarnings(summary(fit_end)))

            writeLines(c ("----------------------",
                          "Checking  Assumptions",
                          "---------------------"))



            ifelse(Normal$p.value <= Alpha,
                   print.noquote (" The errors do not follow a normal distribution."),
                   print.noquote (" The errors follow a normal distribution."))

            ifelse( variance_p <= Alpha,
                    print.noquote(" The variance of the errors is not constant (Heteroscedastic)."),
                    print.noquote(" The variance of the errors is constant (Homoscedastic)."))

            list(
              coefficients = fit_end$coefficients,
              residuals = fit_end$residuals,
              fitted.values = fit_end$fitted.values,
              rank = fit_end$rank,
              df.residual = fit_end$df.residual,
              call = fit_end$call,
              terms = fit_end$terms,
              model = fit_end$model,
              lambda = lambda
            )
          }
          else if (beta0 > Alpha && RowNoSig > 0) {
            delete.Intercept <- row.names(Sig)
            Newdata1 <- Newdata[delete.Intercept]

            while(RowNoSig > 0) {
              Newdata1 <- Newdata[delete.Intercept]
              Newdata1 <- data.frame(y,Newdata1)
              fit1 <- suppressWarnings(lm(y~0+.,data=Newdata1))
              sumfit1 <- suppressWarnings(summary(fit1))
              Number_Beta1 <- as.numeric(nrow(sumfit1$coefficients))

              Decision1 <- 0
              i <- 1
              while( i <= Number_Beta1 ) {
                x <- ifelse(sumfit1$coefficients[i,4]<= Alpha,"Sig","NoSig")
                Decision1[i] <- x
                i = i+1
              }

              view1 <- data.frame(sumfit1$coefficients)
              view1 <- mutate(view1, Decision1)
              Sig <- subset(view1,view1$Decision == "Sig")

              delete.Intercept <- row.names(Sig)
              NoSig <- subset(view1,view1$Decision == "NoSig")
              RowNoSig <- as.numeric(nrow(NoSig))
              RowNoSig = RowNoSig+0
            }

            error <- fit1$residuals
            error.Group <- factor(error<=median(error))
            Normal <- shapiro.test(error)
            variance <- leveneTest(error,group = error.Group)
            variance_p <- variance$`Pr(>F)`[1]

            plot(fit1,1)
            plot(fit1,2)

            print(suppressWarnings(summary(fit1)))

            writeLines(c ("----------------------",
                          "Checking  Assumptions",
                          "---------------------"))


            ifelse(Normal$p.value <= Alpha,
                   print.noquote (" The errors do not follow a normal distribution."),
                   print.noquote (" The errors follow a normal distribution."))

            ifelse( variance_p <= Alpha,
                    print.noquote(" The variance of the errors is not constant (Heteroscedastic)."),
                    print.noquote(" The variance of the errors is constant (Homoscedastic)."))

            list(
              coefficients = fit1$coefficients,
              residuals = fit1$residuals,
              fitted.values = fit1$fitted.values,
              rank = fit1$rank,
              df.residual = fit1$df.residual,
              call = fit1$call,
              terms = fit1$terms,
              model = fit1$model,
              lambda = lambda
            )


          }
          else {
            delete.Intercept <- row.names(Sig)
            Newdata1 <- Newdata[delete.Intercept]
            Newdata1 <- data.frame(y,Newdata1)
            fit1 <- suppressWarnings(lm(y~0+.,data=Newdata1))
            sumfit1 <- suppressWarnings(summary(fit1))

            error <- sumfit1$residuals
            error.Group <- factor(error<=median(error))
            Normal <- shapiro.test(error)
            variance <- leveneTest(error,group = error.Group)
            variance_p <- variance$`Pr(>F)`[1]

            plot(fit1,1)
            plot(fit1,2)

            print(suppressWarnings(summary(fit1)))

            writeLines(c ("----------------------",
                          "Checking  Assumptions",
                          "---------------------"))



            ifelse(Normal$p.value <= Alpha,
                   print.noquote (" The errors do not follow a normal distribution."),
                   print.noquote (" The errors follow a normal distribution."))

            ifelse( variance_p <= Alpha,
                    print.noquote(" The variance of the errors is not constant (Heteroscedastic)."),
                    print.noquote(" The variance of the errors is constant (Homoscedastic)."))

            list(
              coefficients = fit1$coefficients,
              residuals = fit1$residuals,
              fitted.values = fit1$fitted.values,
              rank = fit1$rank,
              df.residual = fit1$df.residual,
              call = fit1$call,
              terms = fit1$terms,
              model = fit1$model,
              lambda = lambda
            )


          }

        }




      }
      else if (Normal$p.value >= Alpha && variance_p >= Alpha) {


        writeLines(c ("-------------------------",
                      "Stepwise regression model ",
                      "-------------------------"))

        print(suppressWarnings(summary(fit1)))
        plot(fit1,1)
        plot(fit1,2)

        writeLines(c("----------------------",
                     "Checking  Assumptions",
                     "---------------------"))

        writeLines ("[1] The errors follow a normal distribution.")
        writeLines ("[1] The variance of the errors is constant (Homoscedastic).")
        list(
          coefficients = fit1$coefficients,
          residuals = fit1$residuals,
          fitted.values = fit1$fitted.values,
          rank = fit1$rank,
          df.residual = fit1$df.residual,
          call = fit1$call,
          terms = fit1$terms,
          model = fit1$model

        )
      }



    }
    else  {
      delete.Intercept <- row.names(Sig)
      Newdata1 <- Newdata[delete.Intercept]
      Newdata1 <- data.frame(y,Newdata1)
      fit1 <- suppressWarnings(lm(y~0+.,data=Newdata1))
      sumfit1 <- suppressWarnings(summary(fit1))

      error <- sumfit1$residuals
      error.Group <- factor(error<=median(error))
      Normal <- shapiro.test(error)
      variance <- leveneTest(error,group = error.Group)
      variance_p <- variance$`Pr(>F)`[1]

      if (Normal$p.value <= Alpha || variance_p <= Alpha) {
        Newdata1 <- data.frame(y,Newdata1)
        Find.Lambda <- (boxcox(y~., data = Newdata1))
        Find.Lambda2 <- data.frame(Find.Lambda$x,Find.Lambda$y)
        Find.Lambda3 <- subset(Find.Lambda2,Find.Lambda$y==max(Find.Lambda$y))
        Optimal.lambda <- Find.Lambda3$Find.Lambda.x
        Optimal.lambda.true <- Find.Lambda3$Find.Lambda.x

        ifelse (Optimal.lambda >= 0.75 && Optimal.lambda < 1.5,Optimal.lambda <- 1,
                Optimal.lambda <- Optimal.lambda)

        if (Optimal.lambda == 1 ) {

          writeLines(c("-------------------------",
                       "Stepwise regression model ",
                       "-------------------------"))

          print(suppressWarnings(summary(fit1)))

          writeLines(c("----------------------",
                       "Checking  Assumptions",
                       "---------------------"))

          writeLines ("[1] The errors do not follow a normal distribution.")
          writeLines ("[1] The variance of the errors is not constant (Heteroscedastic).")

          writeLines("------------------------\n Box cox transformation \n------------------------")
          writeLines("Optimal lambda approximate to 1.")

          message("The multiple linear regression may not be appropriate for this data.")

          plot(fit1,1)
          plot(fit1,2)
          list(
            coefficients = fit1$coefficients,
            residuals = fit1$residuals,
            fitted.values = fit1$fitted.values,
            rank = fit1$rank,
            df.residual = fit1$df.residual,
            call = fit1$call,
            terms = fit1$terms,
            model = fit1$model
          )

        }
        if (Optimal.lambda != 1) {

          writeLines(c("--------------------------------------------------",
                       "Regression model derived from data transformations",
                       "--------------------------------------------------"))
          if (Optimal.lambda> -0.25 && Optimal.lambda< 0.25){
            lambda <- 0

            writeLines(c("The lambda value utilized in the data conversion is 0",
                         "and transformation the dependent variable in the following : y = log(y)"))

            y.prime <- log(y)
            Newdata1$y <- NULL
            Newdata1 <- data.frame(y.prime,Newdata1)
            fit_end <- suppressWarnings(lm(y.prime~0+.,data=Newdata1))
            sum_fit_end <- suppressWarnings(summary(fit_end))

          }
          else if (Optimal.lambda >= 0.25 && Optimal.lambda < 0.75 ) {
            lambda <- 0.5

            writeLines(c("The lambda value utilized in the data conversion is 0.5 ",
                         "and transformation the dependent variable in the following : y=sqrt(y)"))

            y.prime <- sqrt(y)
            Newdata1$y <- NULL
            Newdata1 <- data.frame(y.prime,Newdata1)
            fit_end <- suppressWarnings(lm(y.prime~0+.,data= Newdata1))
            sum_fit_end <- suppressWarnings(summary(fit_end))


          }
          else if (Optimal.lambda >= 1.5 && Optimal.lambda <= 2 ) {
            lambda <- 2

            writeLines(c("The lambda value utilized in the data conversion is 2",
                         "and transformation the dependent variable in the following : y=y^2"))

            y.prime <- (y^2)
            Newdata1$y <- NULL
            Newdata1 <- data.frame(y.prime,Newdata1)
            fit_end <- suppressWarnings(lm(y.prime~0+.,data=Newdata1))
            sum_fit_end <- suppressWarnings(summary(fit_end))

          }
          else if (Optimal.lambda <= -0.25 && Optimal.lambda > -0.75 ) {
            lambda <- -0.5

            writeLines(c("The lambda value utilized in the data conversion is -0.5",
                         "and transformation the dependent variable in the following : y=1/sqrt(y)"))

            y.prime <- 1/sqrt(y)
            Newdata1$y <- NULL
            Newdata1 <- data.frame(y.prime,Newdata1)
            fit_end <- suppressWarnings(lm(y.prime~0+.,data=Newdata1))
            sum_fit_end <- suppressWarnings(summary(fit_end))

          }
          else if (Optimal.lambda <= -0.75 && Optimal.lambda > -1.5 ) {
            lambda <- -1

            writeLines(c("The lambda value utilized in the data conversion is -1",
                         "and transformation the dependent variable in the following : y=1/y"))

            y.prime <- 1/y
            Newdata1$y <- NULL
            Newdata1 <- data.frame(y.prime,Newdata1)
            fit_end <- suppressWarnings(lm(y.prime~0+.,data =  Newdata1))
            sum_fit_end <- suppressWarnings(summary(fit_end))

          }
          else if (Optimal.lambda <= -1.5 && Optimal.lambda >= -2 ) {
            lambda <- -2

            writeLines(c("\n","The lambda value utilized in the data conversion is -2 ",
                         " and transformation the dependent variable in the following : y=1/y^2"))

            y.prime <- 1/(y^2)
            Newdata1$y <- NULL
            Newdata1 <- data.frame(y.prime,Newdata1)
            fit_end <- suppressWarnings(lm(y.prime~0+.,data =  Newdata1))
            sum_fit_end <- suppressWarnings(summary(fit_end))

          }

          Number_Beta <- as.numeric(nrow(sum_fit_end$coefficients))
          Decision <- 0
          i <- 1
          while (i <= Number_Beta) {
            x <- ifelse(sum_fit_end$coefficients[i,4]<=Alpha,"Sig","NoSig")
            Decision[i] <- x
            i = i+1
          }

          view <- data.frame(sum_fit_end$coefficients)
          view <- mutate(view, Decision )
          NoSig <- subset(view,view$Decision == "NoSig")
          RowNoSig <- as.numeric(nrow(NoSig))

          Sig <-  subset(view,view$Decision == "Sig")

          beta0 <- sum_fit_end$coefficients[1,4]

          if (beta0 <= Alpha && RowNoSig > 0) {
            delete.Intercept <- row.names(Sig)
            delete.Intercept <- delete.Intercept[-1]
            Newdata1 <- Newdata[delete.Intercept]

            while(RowNoSig > 0) {
              Newdata1 <- Newdata[delete.Intercept]
              Newdata1 <- data.frame(y,Newdata1)
              fit1 <- suppressWarnings(lm(y~., data=Newdata1))
              sumfit1 <- suppressWarnings(summary(fit1))
              Number_Beta1 <- as.numeric(nrow(sumfit1$coefficients))

              Decision1 <- 0
              i <- 1
              while( i <= Number_Beta1 ) {
                x <- ifelse(sumfit1$coefficients[i,4]<= Alpha,"Sig","NoSig")
                Decision1[i] <- x
                i = i+1
              }

              view1 <- data.frame(sumfit1$coefficients)
              view1 <- mutate(view1, Decision1)
              Sig <- subset(view1,view1$Decision == "Sig")

              delete.Intercept <- row.names(Sig)
              delete.Intercept <- delete.Intercept[-1]
              NoSig <- subset(view1,view1$Decision == "NoSig")
              RowNoSig <- as.numeric(nrow(NoSig))
              RowNoSig = RowNoSig+0
            }

            error <- fit1$residuals
            error.Group <- factor(error<=median(error))
            Normal <- shapiro.test(error)
            variance <- leveneTest(error,group = error.Group)
            variance_p <- variance$`Pr(>F)`[1]

            plot(fit1,1)
            plot(fit1,2)

            print(suppressWarnings(summary(fit1)))

            writeLines(c ("----------------------",
                          "Checking  Assumptions",
                          "---------------------"))


            ifelse(Normal$p.value <= Alpha,
                   print.noquote (" The errors do not follow a normal distribution."),
                   print.noquote (" The errors follow a normal distribution."))

            ifelse( variance_p <= Alpha,
                    print.noquote(" The variance of the errors is not constant (Heteroscedastic)."),
                    print.noquote(" The variance of the errors is constant (Homoscedastic)."))

            list(
              coefficients = fit1$coefficients,
              residuals = fit1$residuals,
              fitted.values = fit1$fitted.values,
              rank = fit1$rank,
              df.residual = fit1$df.residual,
              call = fit1$call,
              terms = fit1$terms,
              model = fit1$model,
              lambda = lambda
            )
          }
          else if (beta0 <= Alpha && RowNoSig <= 0) {
            error_end <- fit_end$residuals
            error.Group <- factor(error_end<=median(error_end))
            Normal <- shapiro.test(error_end)
            variance <- leveneTest(error_end,group = error.Group)
            variance_p <- variance$`Pr(>F)`[1]
            plot(fit_end,1)
            plot(fit_end,2)

            print(suppressWarnings(summary(fit_end)))

            writeLines(c ("----------------------",
                          "Checking  Assumptions",
                          "---------------------"))



            ifelse(Normal$p.value <= Alpha,
                   print.noquote (" The errors do not follow a normal distribution."),
                   print.noquote (" The errors follow a normal distribution."))

            ifelse( variance_p <= Alpha,
                    print.noquote(" The variance of the errors is not constant (Heteroscedastic)."),
                    print.noquote(" The variance of the errors is constant (Homoscedastic)."))

            list(
              coefficients = fit_end$coefficients,
              residuals = fit_end$residuals,
              fitted.values = fit_end$fitted.values,
              rank = fit_end$rank,
              df.residual = fit_end$df.residual,
              call = fit_end$call,
              terms = fit_end$terms,
              model = fit_end$model,
              lambda = lambda
            )
          }
          else if (beta0 > Alpha && RowNoSig > 0) {
            delete.Intercept <- row.names(Sig)
            Newdata1 <- Newdata[delete.Intercept]

            while(RowNoSig > 0) {
              Newdata1 <- Newdata[delete.Intercept]
              Newdata1 <- data.frame(y,Newdata1)
              fit1 <- suppressWarnings(lm(y~0+.,data=Newdata1))
              sumfit1 <- suppressWarnings(summary(fit1))
              Number_Beta1 <- as.numeric(nrow(sumfit1$coefficients))

              Decision1 <- 0
              i <- 1
              while( i <= Number_Beta1 ) {
                x <- ifelse(sumfit1$coefficients[i,4]<= Alpha,"Sig","NoSig")
                Decision1[i] <- x
                i = i+1
              }

              view1 <- data.frame(sumfit1$coefficients)
              view1 <- mutate(view1, Decision1)
              Sig <- subset(view1,view1$Decision == "Sig")

              delete.Intercept <- row.names(Sig)
              NoSig <- subset(view1,view1$Decision == "NoSig")
              RowNoSig <- as.numeric(nrow(NoSig))
              RowNoSig = RowNoSig+0
            }

            error <- fit1$residuals
            error.Group <- factor(error<=median(error))
            Normal <- shapiro.test(error)
            variance <- leveneTest(error,group = error.Group)
            variance_p <- variance$`Pr(>F)`[1]

            plot(fit1,1)
            plot(fit1,2)

            print(suppressWarnings(summary(fit1)))

            writeLines(c ("----------------------",
                          "Checking  Assumptions",
                          "---------------------"))


            ifelse(Normal$p.value <= Alpha,
                   print.noquote (" The errors do not follow a normal distribution."),
                   print.noquote (" The errors follow a normal distribution."))

            ifelse( variance_p <= Alpha,
                    print.noquote(" The variance of the errors is not constant (Heteroscedastic)."),
                    print.noquote(" The variance of the errors is constant (Homoscedastic)."))

            list(
              coefficients = fit1$coefficients,
              residuals = fit1$residuals,
              fitted.values = fit1$fitted.values,
              rank = fit1$rank,
              df.residual = fit1$df.residual,
              call = fit1$call,
              terms = fit1$terms,
              model = fit1$model,
              lambda = lambda
            )


          }
          else {
            delete.Intercept <- row.names(Sig)
            Newdata1 <- Newdata[delete.Intercept]
            Newdata1 <- data.frame(y,Newdata1)
            fit1 <- suppressWarnings(lm(y~0+.,data=Newdata1))
            sumfit1 <- suppressWarnings(summary(fit1))

            error <- sumfit1$residuals
            error.Group <- factor(error<=median(error))
            Normal <- shapiro.test(error)
            variance <- leveneTest(error,group = error.Group)
            variance_p <- variance$`Pr(>F)`[1]

            plot(fit1,1)
            plot(fit1,2)

            print(suppressWarnings(summary(fit1)))

            writeLines(c ("----------------------",
                          "Checking  Assumptions",
                          "---------------------"))



            ifelse(Normal$p.value <= Alpha,
                   print.noquote (" The errors do not follow a normal distribution."),
                   print.noquote (" The errors follow a normal distribution."))

            ifelse( variance_p <= Alpha,
                    print.noquote(" The variance of the errors is not constant (Heteroscedastic)."),
                    print.noquote(" The variance of the errors is constant (Homoscedastic)."))

            list(
              coefficients = fit1$coefficients,
              residuals = fit1$residuals,
              fitted.values = fit1$fitted.values,
              rank = fit1$rank,
              df.residual = fit1$df.residual,
              call = fit1$call,
              terms = fit1$terms,
              model = fit1$model,
              lambda = lambda
            )


          }


        }




      }
      else if (Normal$p.value >= Alpha && variance_p >= Alpha) {


        writeLines(c("-------------------------",
                     "Stepwise regression model",
                     "-------------------------"))

        print(suppressWarnings(summary(fit1)))
        plot(fit1,1)
        plot(fit1,2)

        writeLines(c ("----------------------",
                      "Checking  Assumptions",
                      "---------------------"))

        writeLines("[1] The errors follow a normal distribution")
        writeLines("[1] The variance of the errors is constant (Homoscedastic).")
        list(
          coefficients = fit1$coefficients,
          residuals = fit1$residuals,
          fitted.values = fit1$fitted.values,
          rank = fit1$rank,
          df.residual = fit1$df.residual,
          call = fit1$call,
          terms = fit1$terms,
          model = fit1$model
        )
      }


    }
  }


}


