install.packages("meta")
library(meta)

logOR=c(-0.4741,-0.5605)###beta值填入

selogOR=c(0.26358107,0.18854499) ##se值填入

meta=metagen(logOR, selogOR, sm = "OR")
meta
