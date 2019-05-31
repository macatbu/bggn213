# This is Margaret Burns' submission for Class 6 Homework of BGGN 213


#install relevant packages needed and load the data for my function 

install.packages("bio3d", dependencies=TRUE)

library(bio3d)

# input of my function is three PDB protein codes for the proteins that you want to compare 

# output of my function is a plot of each of the input proteins 

# defining my function below

my_function <- function(protein_name, chain_name) {
  
  # first the function calls in the PDB database for the protein of interest
  s1 <- read.pdb(protein_name) 

  # then we trim down the large PDB file for just the chain of the protein that
  # we are interested in analyzing 
  s1.chainx <- trim.pdb(s1, chain = chain_name, elety="CA")
  
  # from this chain we pull out the atoms we want
  s1.b <- s1.chainx$atom$b

  # plot our region of interest
  return(plotb3(s1.b, sse=s1.chainx, typ="l", ylab="Bfactor"))
}

# practice by callin the function for "4AKE"

my_function("4AKE", "A")
my_function("1AKE", "A")
my_function("1E4Y", "A")
