score=c(278,272,276,281,279,276,281,289,280,276,272,274,278,279,280,270,283,279,273)
temp=c(70,70,75,61,75,74,79,72,75,71,82,79,80,81,73,83,72,74,78)
plot(temp,score)
lm(score~temp)
## Temperature and Score Robust Regression Models

## Data and line without outliers

score2=c(278,276,281,279,276,281,280,276,272,274,278,279,280,270,283,279,273)
temp2=c(70,75,61,75,74,79,75,71,82,79,80,81,73,83,72,74,78)
plot(temp2, score2)

lm(score2~temp2)

## Robust-resistant line

rrline1 <- function(x,y) {
  n <- length(x); nmod3 <- n%%3;
  if (nmod3 == 0) n3 <- n/3;
  if (nmod3 == 1) n3 <- (n - 1)/3;
  if (nmod3 == 2) n3 <- (n + 1)/3;
  x.order <- order(x)
  medxL <- median(x[x.order][1:n3])
  medxR <- median(rev(x[x.order])[1:n3])
  medyL <- median(y[x.order][1:n3])
  medyR <- median(rev(y[x.order])[1:n3])
  slope1 <- (medyR - medyL)/(medxR - medxL)
  int1 <- median(y - slope1 * x)
  newy <- y - slope1 * x - int1
  sumres <- sum(abs(newy))
  list(a = int1, b = slope1, sumres = sumres, res = newy) }

run.rrline <- function (xx, yy, iter = 10) {
  out.coef <- matrix(0, iter, 3)
  ll <- (1:length(xx))[!is.na(xx) & !is.na(yy)]
  n <- length(ll); x <- xx[ll]; y <- yy[ll];
  newy <- y
  for (i in 1:iter) {
    rr <- rrline1(x, newy)
    out.coef[i, ] <- c(rr$a, rr$b, rr$sumres)
    newy <- rr$res }
  dimnames(out.coef) <- list(format(1:iter),c("a","b","|res|"))
  aa <- sum(out.coef[, 1]); bb <- sum(out.coef[, 2]);
  cc <- sum(abs(y - aa - bb * x)^2)
  res <- (yy - aa - bb * xx)
  out.coef <- rbind(out.coef, c(aa, bb, cc))
  print(round(out.coef, 5))
  list(a = aa, b = bb, res = res, coef = out.coef) }

run.rrline(temp, score)
run.rrline(temp, score, iter=10)

## M-estimation

library("L1pack")

abs.fit=lad(score~temp)
abs.fit$coefficients

## Bi-Square

library("MASS")

bisq.fit<-rlm(score~temp, method="MM")
summary(bisq.fit)

## Huber

huber.fit<-rlm(score~temp, method="M")
summary(huber.fit)

# Plot all lines

plot(temp, score)
abline(a=lm(score~temp)$coefficient[1], b=lm(score~temp)$coefficient[2], col="black")
abline(a=run.rrline(temp, score, iter=3)$a, b=run.rrline(temp, score, iter=3)$b, col="blue")
abline(a=abs.fit$coefficients[1], b=abs.fit$coefficients[2], col="red", lwd=2)
abline(a=bisq.fit$coefficients[1], b=bisq.fit$coefficients[2], col="green", lwd=2)
abline(a=huber.fit$coefficients[1], b=huber.fit$coefficients[2], col="orange", lwd=2)
legend("bottomleft", c("OLS", "RR Line", "M-Estimation", "Bi-square", "Huber"), lwd=2, col=c("black", "blue", "red", "green", "orange"), cex = 0.75)