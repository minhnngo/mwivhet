set.seed(2024)
devtools::load_all()
# Judge Data (No Covariates)
df_judge <- GenData_nocov()

# Keep only observable variables for dataset with no covariates (dnc)
dnc <- data.frame(
  group = df_judge$group,
  X = df_judge$X,
  Y = df_judge$Y,
  e = df_judge$e,
  MX = df_judge$MX,
  Me = df_judge$Me,
  MY = df_judge$MY
)

# Interacted QOB Data (With Covariates)
df_cov <- GenData_cov()

# Keep only observable variables for dataset with covariates (dc)
dc <- data.frame(
  group = df_cov$group,
  groupW = df_cov$groupW,
  X = df_cov$X,
  Y = df_cov$Y,
  e = df_cov$e,
  MX = df_cov$MX,
  Me = df_cov$Me,
  MY = df_cov$MY
)

# Saving data
usethis::use_data(dnc, dc, overwrite = TRUE)
