a <- scan("shakespeare.txt",what="character",skip=83,nlines=196043-83,
          fileEncoding="UTF-8")

# STEP 4.a =====================================================================

loc_1 <- grep("[",a,fixed=TRUE) #location/index of [
length_a <- length(a)

all_loc <- NULL
for (i in loc_1){
  loc_2 <- grep("]",a[i:min((i+100),length_a)],fixed=TRUE)
  loc_3 <- grep(".",a[i:min((i+100),length_a)],fixed=TRUE)
  
  if(length(loc_2)>0){
    loc_exc <- c(i:(loc_2[1]+(i)-(1)))
  } else{loc_exc <- c(i:(loc_3[1]+(i)-(1)))}
  
  ifelse (i==loc_1[1], 
          all_loc <- loc_exc, 
          all_loc <- append(all_loc,loc_exc))
}

a <- a[-(all_loc)]