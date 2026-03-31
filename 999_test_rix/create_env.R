# using rix to make a reproducible nix environment for the given libraries
# https://docs.ropensci.org/rix/articles/project-environments.html

# run this R script first to generate default.nix and start-rstudio.sh

library(rix)

rix(r_ver = "4.5.0",
    r_pkgs = c("Rsubread", "tidyverse", "plyranges", "DESeq2", "AnnotationHub", "Rsamtools"),
    system_pkgs = NULL,
    git_pkgs = NULL,
    ide = "rstudio",
    project_path = ".",
    overwrite = TRUE,
    print = TRUE)

rix::make_launcher("rstudio", project_path = ".")

# now run $ nix-build, then launch rstudio with $ sudo chmod +x start-rstudio.sh, then $ ./start-rstudio.sh