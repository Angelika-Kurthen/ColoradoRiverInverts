# code for exporting into numpy 
install.packages("reticulate")
library(reticulate)
np = import("numpy")
a <- np$array(outarray)
pa <- r_to_py(a)
py_save_object(pa, "test.npy")
np$load
np$load("test.npy", allow_pickle=T)
