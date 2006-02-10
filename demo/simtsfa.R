#       Various (important working) pieces have been moved to tests.
#  Note - This demo file is the "working experiment" 


require("tsfa")  
# tfsa requires "tframe", "dse1"(simulate), "setRNG", and "dse2"(EstEval)
data("CanadianMoneyData.asof.6Feb2004", package="CDNmoney")


# This is not really a valid seed for this rng kind, but generates one.
rngValue10 <- list(seed=10, kind="Mersenne-Twister", normal.kind="Inversion")
rngValue20 <- list(seed=20, kind="Mersenne-Twister", normal.kind="Inversion")
rngValue30 <- list(seed=30, kind="Mersenne-Twister", normal.kind="Inversion")


######################################################################

#Generate simulated data using factors based on real per capita  total M1 
# as transaction factor and real per capita M2++ as savings, and then estimate.

######################################################################

#  Move stuff from Calculations.R   !!!!!!!!!!!!!!!!!!
