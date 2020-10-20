
library(deltaGLMtree)
options(na.action = "na.fail")

context("test deltaGLMtree")
test_that("test deltaGLMtree", {
  data(sim_dat)
  dat=sim_dat

  area <- rep(1, nrow(dat))
  area[dat$Lon>167.5] <- 2
  area[dat$Lon<=167.5 & dat$Lat>40] <- 3
  area[dat$Lon<=167.5 & dat$Lat<=40 & dat$Lon>160] <- 4
  area[dat$Lon<=167.5 & dat$Lat>40 & dat$Lon>147.5] <- 5

  dat$Area <- area
  dat$CPUE01 <- ifelse(dat$CPUE==0,0,1)

  dat_p = dat %>% filter(CPUE>0)

  binom1 = glm(CPUE01 ~ 0 + Year*Area + SST + I(SST^2), data = dat, family = binomial)
  gamma1 = glm(CPUE ~ 0 + Year+Area + SST + I(SST^2), data = dat_p, family = Gamma(log), control=glm.control(maxit=10000))
  gamma2 = glm2::glm2(CPUE ~ 0 + Year+Area + SST + I(SST^2), data = dat_p, family = Gamma(log), control=glm.control(maxit=10000))
  lognorm1 = glm(log(CPUE) ~ 0 + Year+Area + SST + I(SST^2), data = dat_p, family = gaussian(), control=glm.control(maxit=10000))
  lognorm2 = glm2::glm2(log(CPUE) ~ 0 + Year+Area + SST + I(SST^2), data = dat_p, family = gaussian(), control=glm.control(maxit=10000))

  # testthat::expect_equal(as.numeric((gamma1$coefficients - gamma2$coefficients)/gamma1$coefficients),rep(0,length(gamma1$coefficients)),tolerance=1.0e-6)
  # testthat::expect_equal(as.numeric((lognorm1$coefficients - lognorm2$coefficients)/lognorm1$coefficients),rep(0,length(lognorm1$coefficients)),tolerance=1.0e-6)

  testthat::expect_equivalent(gamma1$coefficients, gamma2$coefficients,tolerance=1.0e-6)
  testthat::expect_equivalent(lognorm1$coefficients, lognorm2$coefficients,tolerance=1.0e-6)

  range(dat$Lon)
  range(dat$Lat)


  testres = get(load(system.file("extdata","res_gamma.rda",package = "deltaGLMtree")))
  res = deltaGLMtree(binom1, gamma2, criteria = "BIC", delta = 10, trend = TRUE, cat = TRUE, numeric = c("SST"),
                     min.lon = 125,max.lon = 197.5, min.lat = 22.5,max.lat = 60)
  # warnings()
  # save(res,file="inst/extdata/res_gamma.rda")
  expect_equivalent(res$sep,testres$sep)
  expect_equivalent(res$objective,testres$objective)
  expect_equivalent(res$LonLat,testres$LonLat)
  expect_equal(res$best$n.area,2)
  expect_equivalent(res$best$p.trend0,testres$best$p.trend0)
  expect_equivalent(res$best$m.trend0,testres$best$m.trend0)
  expect_equivalent(res$best$y.trend0,testres$best$y.trend0)
  expect_equivalent(res$best$p.trend,res$best$p.trend0)
  expect_equivalent(res$best$m.trend,res$best$m.trend0/mean(res$best$m.trend0))
  expect_equivalent(res$best$y.trend,res$best$y.trend0/mean(res$best$y.trend0))

  testres2 = get(load(system.file("extdata","res_lognorm.rda",package = "deltaGLMtree")))
  res2 = deltaGLMtree(binom1, lognorm2, criteria = "BIC", delta = 10, trend = TRUE, cat = TRUE, numeric = c("SST"),
                     min.lon = 125,max.lon = 197.5, min.lat = 22.5,max.lat = 60)
  # save(res2,file="inst/extdata/res_lognorm.rda")
  expect_equivalent(res2$sep,testres2$sep)
  expect_equivalent(res2$objective,testres2$objective)
  expect_equivalent(res2$LonLat,testres2$LonLat)
  expect_equal(res2$best$n.area,2)
  expect_equivalent(res2$best$p.trend0,testres2$best$p.trend0)
  expect_equivalent(res2$best$m.trend0,testres2$best$m.trend0)
  expect_equivalent(res2$best$y.trend0,testres2$best$y.trend0)
  expect_equivalent(res2$best$p.trend,res2$best$p.trend0)
  expect_equivalent(res2$best$m.trend,res2$best$m.trend0/mean(res2$best$m.trend0))
  expect_equivalent(res2$best$y.trend,res2$best$y.trend0/mean(res2$best$y.trend0))

  expect_equal(res$best$loglik.binom,res2$best$loglik.binom)
  expect_equal(res$best$loglik.posi > res2$best$loglik.posi,TRUE)
  expect_equal(res$best$loglik.sum > res2$best$loglik.sum,TRUE)

  expect_equal(res$best$binom$rank, res2$best$binom$rank)
  expect_equal(res$best$posi$rank, res2$best$posi$rank)
  })

