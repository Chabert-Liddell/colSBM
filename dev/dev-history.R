### CREATE PACKAGE
use_package_doc()
use_roxygen_md()

### ADD DESCRIPTION FILES
edit_file("DESCRIPTION")
use_mit_license("Saint-Clair Chabert-Liddell")

### ADD DEPENDENCIES (IF APPLICABLE)
use_package("sbm")
use_package("ape")
use_package("stats")
use_package("utils")
use_package("ggplot2")
use_package("purrr")
use_package("dplyr")
use_package("tidyr")
use_package("patchwork")
use_pipe()
use_package("cluster")

### ADD R FUNCTIONS AND THEIR DOCUMENTATION
# use_r("name_of_r_file")
# Code > Insert Roxygen Skeleton (Ctrl+Alt+Shift+R)
use_data_raw()
use_r("data")
use_r("generic-function")
use_vignette("tutorial", "Tutorial on food webs")


### ADD VIGNETTE-LIKE RMD TO YOUR PACKAGE
use_readme_rmd()

### ADD UNITARY TESTS
use_testthat()
use_test("name_of_test_file")

test()


### LINK TO TRAVIS-CI
use_travis()
# Go to travis-ci.org and link your repo

### ADD PAT TO R
# Create PAT: https://github.com/settings/tokens with "public_repo" scope.
# Add the PAT to .Renviron
edit_r_environ()
# GITHUB_PAT=<GENERATED PAT>
git_vaccinate()

# [travis::travis_set_pat() | BUT NEEDS MORE PERMISSIONS]

### SETTING A WEBPAGE FOR YOUR PACKAGE
use_pkgdown()
use_pkgdown_travis()

### CREATE SSH KEYPAIR FOR TRAVIS DEPLOY
# Run use_travis_deploy_manual.R

#----- 'use_travis_deploy', but manual

# Generate SSH keypair
key <- openssl::rsa_keygen()


# Go to the GitHub repository page - then Settings > Deploy keys,
# add the public key there (check the box for write-access).
title <- "travis+tic"

# Go to the Travis repository page - then More Options > Settings.
# Add an environment-variable named id_rsa. Paste your clipboard-contents into its value.
# Make sure you are not displaying the value in the log. Save.
name <- "id_rsa"
private_key

#----- End of 'use_travis_deploy', but manual

# Go to the GitHub repository page - then Settings > Deploy keys, add the public key there (check the box for write-access).
# Go to the Travis repository page - then More Options > Settings. Add an environment-variable named id_rsa. Paste your clipboard-contents into its value. Make sure you are not displaying the value in the log. Save.

# [OR travis::use_travis_deploy() | TOO MANY PERMISSIONS]

### ADD TEST COVERAGE REPORTS
use_coverage(type = "codecov")
# Go to the codecov repository page - then Add a repository > 'repo_name' > Copy Token. Then go to the Travis repository page - then More Options > Settings. Add an environment-variable named CODECOV_TOKEN. Paste your clipboard-contents into its value.


use_version(which = "minor")
usethis::use_version(which = "minor")
