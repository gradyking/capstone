findProximalBinding.R is my R script to try to find whether or not Musashi-1 binds directly to RNA proximally (very closely) to AGO2 to suppress binding. this would mean a fun interaction where Musashi-1 is competitively inhibiting binding or something like that

MSI1-with_input.regions.ucsc.bed has the binding sites of Musashi-1 in mice retinal cells for comparison to AGO2
diff_chimeric_clip.zip contains a tab-separated table of AGO2 binding sites, and how they change in wild-type versus Musashi-1 knockout
forReference--CSD_CLIP.r has an excerpt of R code given by Dr. Stoilov that helps convert from binding sites to actual DNA 
