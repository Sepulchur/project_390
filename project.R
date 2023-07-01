ms_obs_lines = readLines("ms_obs_final.out")
pars_final =as.numeric(readLines("pars_final.txt"))

referenceDataset <- matrix(nrow = 50, ncol = 132)
for (i in 1:50) {
  for (j in 1:132) {
    referenceDataset[i, j] <- substr(ms_obs_lines[i], j, j)
  }
}

# Read the text file into a character vector
data <- readLines("ms_sim_final.out")

# Create an empty list to store the individual matrices
simulatedDatasets <- list()

# Split the data into individual datasets
datasets <- split(data, cumsum(nchar(data) == 0))

# Loop through each dataset and convert them into matrices
for (dataset in datasets) {
  if (length(dataset) > 0) {
    # Remove empty lines from the dataset
    dataset <- dataset[nchar(dataset) > 0]
    
    # Create a matrix for the dataset
    matrix_data <- do.call(rbind, strsplit(dataset, ""))
    simulatedDatasets[[length(simulatedDatasets) + 1]] <- matrix_data
  }
}

# Function to calculate pairwise differences
pairwise_differences <- function(seq1, seq2) {
  differences <- sum(seq1 != seq2)
  return(differences)
}

#calculate the number of differences between the observed dataset and every simulated dataset

pairdiffobserve = function(dataset ){
  n = 50
  # Create a matrix of all pairwise comparisons
  pairwise_diff <- matrix(0, n, n)
  
  # Calculate the pairwise differences
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      pairwise_diff[i, j] <- sum(dataset[i, ] != dataset[j, ])
      pairwise_diff[j, i] <- pairwise_diff[i, j] # Fill in the symmetric value
    }
  }
  
  # Sum up all the pairwise differences
  sum_differences <- sum(pairwise_diff)
  return (sum_differences)
}

#calculate k=the k statistic for the observed dataset

differnces = pairdiffobserve(referenceDataset)
diffsum = sum(differnces)


paronomastis = choose(50,2)
k0 = diffsum / paronomastis
k = c()
for(i in 1 : 10000){
  dataset = simulatedDatasets[i]
  dataset = unlist(dataset)
  dataset = matrix(dataset , nrow=50)
  k = c(k , pairdiffobserve(dataset) )
}
k = k / paronomastis
#calculate the w statistic for every simulated dataset 
calculatew = function(dataset_list){
  a = c()
  w = c()
  num_datasets = length(dataset_list)
  for(i in 1:49){
    x=1/i
    a = c(a,x)
  }
  a = sum(a)
  
  for(i in 1:num_datasets){
    dataset = dataset_list[[i]]
    S = ncol(dataset)
    x = S/a
    w=c(w,x)
  }
  
  return(w)
}

w = calculatew(simulatedDatasets)

#calculate the w statistic for the observed dataset and the a1 variable
a1=c()
for(i in 1:49){
  x=1/i
  a1 = c(a1,x) 
}
a1 = sum(a1)

w0 = ncol(referenceDataset)/a1

#calculate the a2 variable 
a2 = c()
for(i in 1:49){
  x=1/(i^2)
  a2 = c(a2,x)
}
a2 = sum(a2)

#calculate the b1 and b2 variable 
b1 = (50+1)/(3*(50-1))
b2 = 2*((50^2) + 50 + 3 ) /((9*50)*(50-1))

#calculate the c1 and c2 variables
c1 = b1 - (1/a1)
c2 = b2 - ((50+1) / (a1*50)) + (a2/(a1^2))

#calculate the e1 and e2 variables 
e1 = c1/a1
e2 = c2 / ((a1^2) +a2)


#calculate the D statistic for all the simulated datasets
calculateD = function(k,w,e1,e2,dataset_list){
  D = c()
  num_datasets = length(dataset_list)
  for(i in 1:num_datasets){
    dataset = dataset_list[[i]]
    S = ncol(dataset)
    paronomastis = (e1*S) + ((e2*S)*(S-1))
    paronomastis = sqrt(paronomastis)
    arithmitis = k[i] - w[i]
    x = arithmitis / paronomastis
    D=c(D,x)
  }
  return(D)
}

D = calculateD(k,w,e1,e2,simulatedDatasets)

#calculate D0
S = ncol(referenceDataset)
arithmitis = k0 - w0
paronomastis = (e1*S) + ((e2*S)*(S-1))
paronomastis = sqrt(paronomastis)
D0 = arithmitis / paronomastis 

#normalize each vector 
meank = mean(k)
meanw = mean(w)
meanD = mean(D)
vark = var(k)
varw = var(w)
varD = var(D)

knorm = (k - meank) / vark
k0norm = (k0 - meank) / vark
wnorm = (w - meanw) /varw
w0norm = (w0 - meanw)/varw
Dnorm = (D - meanD)/varD
D0norm = (D0 - meanD) / varD

#calculate the Euklidean distances 
d = c()
for(i in 1:10000){
  x = (D0norm - Dnorm[i])^2
  y = (w0norm - wnorm[i])^2
  z = (k0norm - knorm[i])^2
  h = x+y+z
  h = sqrt(h)
  d=c(d,h)
}

#get the 500 smallest distances
indexes = order(d)[1:500]
dsmallest = d[indexes]

#match the distances with the values in pars_final.txt
corresponding = pars_final[indexes]

#calculate the mean and median for the 500 smallest distances 

meancorr = mean(corresponding)
mediancorr = median(corresponding)

#make the histogram and density plot 

hist(corresponding)
plot(density(corresponding))

