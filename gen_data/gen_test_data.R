library(EigenH5)
N <- 50
p <- 40000
sample_ids <- as.character(sample(1:(N*p),N,replace=F))
ntsnpf <- paste0("test_data.txt.gz")
tsnp_mat <- matrix(sprintf("%.3f",runif(min=0,max=2,N*p)),nrow = N,byrow=T)
# expect_true(all(nchar(tsnp_mat)==5))
wtsnp_mat <- cbind(cbind(sample_ids,rep("DOSE",N)),tsnp_mat)
readr::write_delim(tibble::as_data_frame(wtsnp_mat[sample(1:N),]),path = ntsnpf,delim = "\t",col_names = F,append = F)
test_input_f <- "test_inp.h5"
if(file.exists("test_inp.h5")){
  file.remove("test_inp.h5")
}
if(file.exists("test_out.h5")){
  file.remove("test_out.h5")
}
atsnp_mat <-wtsnp_mat[,-c(1,2)]
class(atsnp_mat) <- "numeric"
sample_names <- wtsnp_mat[,1]
p <- ncol(atsnp_mat)
t_idx <- sort(sample(1:p,min(100,as.integer(p/2)),replace=F))
ttsnp_mat <- t(atsnp_mat[,t_idx])
attr(ttsnp_mat,"dimnames") <- NULL
write_vector_h5(sample_names,filename = test_input_f,"sample_id")
write_vector_h5(t_idx-1,test_input_f,"snp_id")
system("../cmake-build-debug/mach2h5 --machfile=/home/nwknoblauch/Dropbox/Repos/mach2h5/gen_data/test_data.txt.gz  --hdf5=/home/nwknoblauch/Dropbox/Repos/mach2h5/gen_data/test_out.h5 --samplenames=/home/nwknoblauch/Dropbox/Repos/mach2h5/gen_data/test_inp.h5 --SNPlist=/home/nwknoblauch/Dropbox/Repos/mach2h5/gen_data/test_inp.h5")
ret_mat <- read_matrix_h5("test_out.h5","dosage")
expect_equal(t(ret_mat),ttsnp_mat)
