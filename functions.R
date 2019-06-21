preprocessing <- function(mutations ,mut_one, mut_two, answer, nucleotide1, nucleotide2){
  data_need_mut <-  rbind(mutations[mutations['name'] == mut_one,],   #mut_one = mutation, which we want to edit
                          mutations[mutations['name'] == mut_two,])   #mut_two = complement mutation, which we want to edit
  
  answer = c(answer,nrow(data_need_mut))
  
  cic <-c()
  i <- 1
  while (i <= nrow(data_need_mut)){
    if ( grepl("athogenic", data_need_mut$clinSign[i]) == T   &&     grepl("Conflict", data_need_mut$clinSign[i]) == F  ){
      cic = c(cic, i)
    }
    i <- i + 1
  }
  data_need_mut = data_need_mut[cic, ]
  
  
  
  i <- 1
  while ( i <= nrow(data_need_mut)){
    data_need_mut$nucleotide[i] <- seq[[1]][data_need_mut$chromEnd[i]:data_need_mut$chromEnd[i]]
    i <- i + 1
  }
  
  data_need_mut["MedGen"] <- NA
  i <- 1
  while(i<=nrow(data_need_mut)) {
    place_medgen <- gregexpr(pattern ='MedGen:',data_need_mut$phenotype[i])
    list_medgen <- c(substr(data_need_mut$phenotype[i], place_medgen[[1]][1]+7, place_medgen[[1]][1]+14))
    g <- 2
    while(g<=length(place_medgen[[1]])) {
      list_medgen <- append(list_medgen, substr(data_need_mut$phenotype[i], place_medgen[[1]][g]+7, place_medgen[[1]][g]+14) )
      g <- g + 1
    }
    data_need_mut$MedGen[i] <- paste(list_medgen, collapse=" ")
    i <- i + 1
  }
  
  
  answer = c(answer, nrow(data_need_mut[data_need_mut['nucleotide'] == nucleotide1,]),nrow(data_need_mut[data_need_mut['nucleotide'] == nucleotide2,]) )
  return(list(data_need_mut[data_need_mut['nucleotide'] == nucleotide1,],data_need_mut[data_need_mut['nucleotide'] == nucleotide2,],answer))

  }

