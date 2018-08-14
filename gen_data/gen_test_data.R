library(EigenH5)
library(testthat)
N <- 50
p <- 40000
sample_ids <- as.character(sample(1:(N*p), N, replace=F))
ntsnpf <- file.path(getwd(), paste0("test_data.txt.gz"))
tsnp_mat <- matrix(sprintf("%.3f", runif(min=0, max=2, N*p)), nrow = N, byrow=T)
wtsnp_mat <- cbind(cbind(sample_ids, rep("DOSE", N)), tsnp_mat)
readr::write_delim(tibble::as_data_frame(wtsnp_mat[sample(1:N), ]), path = ntsnpf, delim = "\t", col_names = F, append = F)
test_input_f <- file.path(getwd(), "test_inp.h5")
test_output_f <- file.path(getwd(), "test_out.h5")
if(file.exists(test_input_f)){
  file.remove(test_input_f)
}
if(file.exists(test_output_f)){
  file.remove(test_output_f)
}
atsnp_mat <-wtsnp_mat[, -c(1, 2)]
class(atsnp_mat) <- "numeric"
sample_names <- wtsnp_mat[, 1]
p <- ncol(atsnp_mat)
t_idx <- sort(sample(1:p, min(100, as.integer(p/2)), replace=F))
ttsnp_mat <- t(atsnp_mat[, t_idx])
attr(ttsnp_mat, "dimnames") <- NULL
write_vector_h5(sample_names, filename = test_input_f, "sample_id")
write_vector_h5(t_idx-1, test_input_f, "snp_id")
systemcmd <- paste0("../cmake-build-debug/mach2h5 ",
                    "--machfile=", ntsnpf,
                    " --hdf5=", test_output_f,
                    " --samplenames=", test_input_f,
                    " --SNPlist=", test_input_f)
#"../cmake-build-debug/mach2h5 --machfile=/home/nwknoblauch/Dropbox/Repos/mach2h5/gen_data/test_data.txt.gz  --hdf5=/home/nwknoblauch/Dropbox/Repos/mach2h5/gen_data/test_out.h5 --samplenames=/home/nwknoblauch/Dropbox/Repos/mach2h5/gen_data/test_inp.h5 --SNPlist=/home/nwknoblauch/Dropbox/Repos/mach2h5/gen_data/test_inp.h5"
system(systemcmd)
ret_mat <- read_matrix_h5(test_output_f, "dosage")
expect_equal(t(ret_mat), ttsnp_mat)
