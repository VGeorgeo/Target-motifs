find_PAM_G <- function(data, start, finish, PAM, type_PAM, len_pam, nucleotide, length_RNA ,answer){
  data['type_PAM'] <- type_PAM
  data['PAM_search'] <- NA
  data['PAM'] <- NA
  data['sequence'] <- NA
  data['frequency'] <- NA
  data['mut_location'] <- NA
  data['RNA'] <- NA
  data['PAM_seq'] <- NA
  i <- 1
  while ( i <= nrow(data)){
    finish_wind_PAM = data$chromEnd[i] + length_RNA - start
    start_wind_PAM = finish_wind_PAM - len_pam - (finish - start - 1)
    data$PAMsearch[i] <- gsub(", ","",toString(seq[[1]][start_wind_PAM:finish_wind_PAM]))
    data$PAM[i] <- regexpr( PAM, data$PAMsearch[i])
    i = i + 1
  }
  data <- data[data$PAM != -1,]
  answer = c(answer,nrow(data))
  data$PAM <- sapply(data$PAM, as.numeric)
  paam = data
  
  i <- 1
  while (i <= nrow(data)){
    zum <- 0
    #seque - if we find PAM, delete left part of seq to find more PAMs
    seque <- substr(data$PAMsearch[i], (data$PAM[i] + 1), nchar(data$PAMsearch[i]))
    if( regexpr(PAM, seque) != -1 ){
      while (nchar(seque) >= len_pam && regexpr(PAM, seque) != -1 ){
        zum <- zum + regexpr(PAM, seque)
        data$PAM[i] <- data$PAM[i] + regexpr(PAM, seque)
        seque <- substr(data$PAMsearch[i], (data$PAM[i] + 1), nchar(data$PAMsearch[i]))
        data[nrow(data) + 1,] = data[c(i), ]
      }
    }
    data$PAM[i] <- data$PAM[i] - zum
    i <-i+1
  }
  data = unique(data)
  answer = c(answer,nrow(data))
  
  i <- 1
  while ( i <= nrow(data)){
    loc_start = (data$chromEnd[i] + length_RNA - start) - ( nchar(data$PAMsearch[i]) - data$PAM[i] - (len_pam - 1) ) 
    open = loc_start -  length_RNA + start
    clothe = open + (finish - start)
    data$sequence[i] <- gsub(", ","",toString(seq[[1]][open:clothe]))
    data$frequency[i] <- length(grep(nucleotide, data$sequence[i]))
    tolower(data$sequence[i])
    data$mut_location[i] <- (finish + 1) - data$PAM[i]
    i <- i + 1
  }
  if (nrow(data)){
    data = data[data['frequency'] == 0,]
  }
  else{
    data = data[c(0),]
  }
  answer = c(answer,length(unique(data$chromEnd)), nrow(data))
  
  i <- 1
  while (i <= nrow(data)){
    clothe <- (data$chromEnd[i] + length_RNA - start) - ( nchar(data$PAMsearch[i]) - data$PAM[i] - (len_pam - 1) ) 
    open <- clothe - length_RNA + 1
    data$RNA[i] <- gsub(", ","",toString(seq[[1]][open:clothe]))
    data$PAM_seq[i] <-  toupper(gsub(", ","",toString(seq[[1]][ (clothe - (len_pam - 1)) :clothe])))
    i <- i + 1
  }
  data <- data[c(1,3,4,14,17,18,27,35,36,39,41,42,43)]
  return(list(data,answer,paam))
}