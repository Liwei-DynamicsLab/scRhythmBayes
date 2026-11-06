dir.create("outputs", showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), "outputs/sessionInfo_current.txt")
packages <- as.data.frame(utils::installed.packages()[, c("Package", "Version")])
utils::write.table(packages, "outputs/packages_current.tsv", sep = "\t",
                   quote = FALSE, row.names = FALSE)
