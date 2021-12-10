install.packages("dplyr")
install.packages("rjson")
install.packages("tidyverse")
install.packages("lsa")
install.packages("sjmisc")
library("sjmisc")
library("rjson")
library("dplyr")
library("tidyverse")
library("lsa")

json_data <- fromJSON(file = "TVs-all-merged 2.json")
#----------------------------------------------------------------------------
#Here a matrix is created with 1 for duplicates and zero otherwise
g<- matrix(0, nrow = 1624, ncol = 1624)
x <- c(1:1262)
j <- 1

for (p in x){
  i <- length(json_data[[p]])
  if (i == 1) {
    j <- j +1
  } else if (i == 2) {
    g[j,j+1] <- 1
    j <- j +2
  } else if (i == 3) {
    g[j,j+1] <- 1
    g[j,j+2] <- 1
    g[j+1,j+2] <- 1
    j <- j +3
  } else if (i == 4) {
    g[j,j+1] <- 1
    g[j,j+2] <- 1
    g[j,j+3] <- 1
    g[j+1,j+2] <- 1
    g[j+1,j+3] <- 1
    g[j+2,j+3] <- 1
    j <- j +4
  }
}

#-----------------------------------------------------------------------------
#Here all the different titles from the data are extracted
x <- c(1:1262)
s <- 0
for (p in x) {
  i <- length(json_data[[p]])
  for (k in 1:i) {
    y <- json_data[[p]][[k]][["title"]]
    if (str_contains(json_data[[p]][[k]][["featuresMap"]][["Brand"]], "")) {
      y <- paste(y, json_data[[p]][[k]][["featuresMap"]][["Brand"]])
    }
    s <- append(s,y)
  }}
s <- s[-1]
titles <- c(s)
titles <- tolower(titles)
titles <- gsub("\"", "-inch", titles)
titles <- gsub("''", "-inch", titles)
titles <- gsub(" hz", "hz", titles)
titles <- gsub(" -inch", "-inch", titles)

#create a vector "brand" with all different brands and create a vector brand.list
# with the brand per product and 0 otherwise 
x <- c(1:1262)
s <- 0
brand <- c()
for (p in x) {
  i <- length(json_data[[p]])
  for (k in 1:i) {
    if (str_contains(json_data[[p]][[k]][["featuresMap"]][["Brand"]], "")) {
      brand <- append(brand, json_data[[p]][[k]][["featuresMap"]][["Brand"]])
    } 
  }}
brand <- tolower(brand)
brand <- unique(brand)
brand.list <- c()
for (p in x) {
  i <- length(json_data[[p]])
  for (k in 1:i) {
    brad <- c(NA)
    for (l in 1:length(brand)) {
      if (str_contains(tolower(json_data[[p]][[k]][["title"]]), brand[l])) {
        brad <- brand[l]
      }}
    if (is.na(brad)) {
      if (str_contains(json_data[[p]][[k]][["featuresMap"]][["Brand"]], "")) {
        brad <- json_data[[p]][[k]][["featuresMap"]][["Brand"]]
      }}
    brand.list <- append(brand.list,brad)
  }}
brand.list <- tolower(brand.list)
for (i in 1:1624) {
  if (is.na(brand.list[i])){
    brand.list[i] <- 0
  }
}

# The start of the function which is used in the bootstrap
begin <- function(usedset, alpha) {
  titles1 <- 0
  brand.list1 <- 0
  for (i in usedset) {
    titles1 <- append(titles1,titles[i])
    brand.list1<- append(brand.list1,brand.list[i])
  }
  brand.list1 <- brand.list1 [-1]
  titles1 <- titles1[-1]
  woord <- c()
  x1 <- c(1:25)
  # Find all model words and add tv brands to the model words
  for (i in x1) {
    woord <-c(woord, word(titles1, i))
  }
  woord <- unique(woord)
  woord <- na.omit (woord)
  m.woord <- c()
  m.woord <- str_extract(woord, "[a-zA-Z0-9]*(([0-9]+[^0-9, ]+)|([^0-9, ]+[0-9]+))[a-zA-Z0-9]*")
  m.woord <- na.omit(m.woord)
  m.woord <- c(m.woord)
  m.woord <- unique(m.woord)
  m.woord <- append(m.woord, brand)
  
  #Create a binary matrix which has a 1 if a model word is present in the specific title   
  b.matrix <- matrix(0,length(usedset),length(m.woord))
  for (p in 1:length(usedset)) {
    titels <- titles1[p]
    end <- FALSE
    i <- 1
    while (end == FALSE) {
      mw <- word(titels, i)
      if(is.na(mw)) {
        end = TRUE
      } else for (j in 1:length(m.woord)) {
        if(mw == m.woord[j]) {
          b.matrix[p,j] <- 1 
        }
      }
      i = i+1
    }
  }
  b.matrix <- t(b.matrix)
  
  #Min hashing. Here different permutations are created to eventually form the signature matrix
  per.matrix <- t(matrix(0, nrow(b.matrix),600))
  for (i in 1:600) {
    set.seed(i)
    per.matrix[i,] <- matrix(sample(seq(1:nrow(b.matrix))))
  }
  per.matrix <- t(per.matrix) 
  signature <- matrix(1400, nrow = ncol(per.matrix), ncol = ncol(b.matrix))
  for (r in 1:nrow(b.matrix)) {
    for (c in 1:ncol(b.matrix)) {
      if (b.matrix[r,c] == 1) {
        for (i in 1:ncol(per.matrix))
          if (per.matrix[r,i] < signature[i,c]) {
            signature[i,c] <- per.matrix[r,i]
          }
      }
    }}
  
  # set bands and row for LSH
  band <- 500 #set brands
  row <- 500/band
  
  # Start the LSH and obtain the candidate pairs
  neighbores <- matrix(0, nrow = ncol(signature), ncol = ncol(signature))
  end <- FALSE
  i <- 1 
  j <- row #number of rows per band
  while ( end == FALSE) {
    band.matrix <- matrix(signature[i:j,], nrow = row, ncol = ncol(signature))
    for (t in 1:ncol(signature)){
      for (s in t:ncol(signature)) {
        if (all(band.matrix[,t] == band.matrix[,s])) {
          neighbores[t,s] <- 1
        } 
      } 
      neighbores[t,t] <- 0
    }
    i <- i + row
    j <- j + row
    if (i > band*row) {    
      end <- TRUE
    }
  }
  
  #--------------------------------------------------------------------------------
  # start by removing all products that do not have the same brand. Second:  
  #By using the cosine similarity between titles of products we find which pairs are real
  # duplicates and which are not
  neighbores1 <- neighbores
  cos <- matrix(0, nrow = ncol(signature), ncol = ncol(signature))
  for (t in 1:ncol(neighbores)) {
    for (s in t:ncol(neighbores)) {
      if (neighbores1[t,s] == 1 & brand.list1[t] != brand.list1[s]) {
        neighbores1[t,s] <- 0
      }
      
      if (neighbores1[t,s] == 1) {
        cos[t,s] <- cosine(signature[,t],signature[,s])
        if(cos[t,s] < alpha) {
          neighbores1[t,s] <- 0
        }
      }
      
    }
  }
  
  
  #------------------------------------------------------------------------------
  
  #Create the matrix containing the real duplicates. with this matrix the algoritm
  #can see how many duplicates that were found are indeed real duplicates
  dupl <- matrix(0,nrow= length(usedset), ncol=length(usedset))
  for (j in 1:length(usedset)) {
    for (i in 1:length(usedset)) {
      dupl[j,i] <- g[usedset[j],usedset[i]]
    }}
  
  
  #Calculate the f1 star and the f1 and return the answers
  new <- dupl + neighbores
  new1 <- dupl + neighbores1
  duplicates.found <- sum(new == 2)
  duplicates.found1 <- sum(new1 == 2) 
  comparissons <- (sum(neighbores))
  fraction.compare <- comparissons/(((length(usedset)^2)-length(usedset))/2)
  pair.completness <- duplicates.found/sum(dupl)
  pair.quality <- duplicates.found/ comparissons
  f1.star <- 2/ ((1/pair.completness) + (1/pair.quality))
  precision <- duplicates.found1/sum(neighbores1)
  recall <- duplicates.found1/sum(dupl)
  f1 <- 2/ ((1/precision) + (1/recall))
  answer <- c(fraction.compare, pair.completness, pair.quality, f1.star, precision, recall, f1)
  return(answer)
}


#Here we let the algoritm do a bootstrap five times
bootstraps <- 5
f1 <- 0
train.answer <- c()
test.answer <- c()
best.para <- c()
for(b in 1:bootstraps){
  train <- c()
  test <- c()
  set.seed(b)
  draws <- sample(1:length(titles),replace=TRUE)
  for(i in 1:length(titles)){
    if(i %in% draws){
      train <- append(train,i)
    } else {
      test <- append(test,i)
    }
  }
  #Here a grid search is performed for different values 
  allpha <- 0
  #  grid <- c(0.9,0.8,0.99)
  grid <- c(0.99)
  null <- 0
  for (i in grid) {
    train.answer <- append(train.answer,begin(train,i))
    if (train.answer[6] > null) {
      null <- train.answer[6]
      allpha <- i
      
    }
    train.answer <- 0 
  }
  best.para <- append(best.para,allpha)
  test.answer <- append(test.answer,begin(test,allpha))
  test.answer
}