set.seed(12345)

test_that("Consistency with hclust on random data", {
  x <- rnorm(10)
  ward_univar <- univarclust:::ward_1d(x)
  reference   <- hclust(dist(x), method = "ward.D2")

  expect_equal(ward_univar$height, reference$height)  
  
  cl_univar    <- as.list(as.data.frame(cutree(ward_univar, 1:9)))
  cl_reference <- as.list(as.data.frame(cutree(reference, 1:9)))
  res <- mapply(aricode::ARI, cl_univar, cl_reference)
  
  expect_true(all(res == 1))
})

test_that("Consistency with hclust on aves data", {
  
  data(aves)
  ward_univar <- univarclust:::ward_1d(aves$weight)
  reference   <- hclust(dist(aves$weight), method = "ward.D2")

  expect_lt(sum((ward_univar$height- reference$height)^2)/sum(reference$height^2), 1e-6)  
  
  cl_univar    <- as.list(as.data.frame(cutree(ward_univar, 1:9)))
  cl_reference <- as.list(as.data.frame(cutree(reference, 1:9)))
  res <- mapply(aricode::ARI, cl_univar, cl_reference)
  
  expect_true(all(res == 1))
})


test_that("Consistency with hclust on aves data with labels", {
  
  data(aves)
  x <- tapply(aves$weight, aves$family, mean)
  ward_univar <- univarclust:::ward_1d(x)
  reference   <- hclust(dist(x), method = "ward.D2")

  expect_lt(sum((ward_univar$height - reference$height)^2)/sum(reference$height^2), 1e-6)  
  
  cl_univar    <- as.list(as.data.frame(cutree(ward_univar, 1:9)))
  cl_reference <- as.list(as.data.frame(cutree(reference, 1:9)))
  res <- mapply(aricode::ARI, cl_univar, cl_reference)
  
  expect_true(all(res == 1))
})


test_that("Consistency with hclust with repeated data", {
  x <- sample(rnorm(10), replace = TRUE)
  ward_univar <- univarclust:::ward_1d(x)
  reference   <- hclust(dist(x), method = "ward.D2")

  expect_equal(ward_univar$height, reference$height)  
  
  cl_univar    <- as.list(as.data.frame(cutree(ward_univar, 1:9)))
  cl_reference <- as.list(as.data.frame(cutree(reference, 1:9)))
  res <- mapply(aricode::ARI, cl_univar, cl_reference)
  
  expect_true(all(res == 1))
})
