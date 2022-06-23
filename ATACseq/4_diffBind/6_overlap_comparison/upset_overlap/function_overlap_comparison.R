bedToVennlist <- function(filename){
  #filename <- in.ls2[[2]]
  #df <- try(read.table(filename, sep='\t', header=F))
  
  df <- tryCatch(read.table(filename, sep='\t', header=F),
                 error=function(e){
                   print(paste0("empty file ", filename))
                   return(NULL)
                 },
                 finally = print(paste0("read ", filename))
  )
  
  if(!is.null(df)){
    df$V4 <- paste0(df$V1,":",df$V2,"-", df$V3)
    return(df$V4)
    
    #df=NULL
  } 
  #df <- read.table(filename, sep='\t', header=F)
  #df$V4 <- paste(df$V1, df$V2, df$V3, sep = '_')
  #return(df$V4)
}

makeUpsetPlots <- function(pa="u.bed",path, w=12, h=5, text.scale = 1.3, 
                           f.n="upset.pdf", sets=NULL) {
  print(pa)
  in.ls <- list.files(path = path, pattern = pa , full.names = T)
  n.ls <- list.files(path = path, pattern = pa , full.names = F)
  names(in.ls) <- lapply(n.ls, function(x) {
    v <- unlist(strsplit(x, '[.]'))
    n <- paste0(v[1],".",v[6])
    return(n)
  })
  
  #pa <- unlist(strsplit(pa, '[.]'))[1]
  #dir.create(output, showWarnings = F)
  f.name <- paste0(path, "/", f.n)
  
  in_bed.ls <- lapply(in.ls, bedToVennlist)
  names(in_bed.ls) <- names(in.ls)
  in_bed.ls <- Filter(Negate(is.null),in_bed.ls)
  print("list finished")
  print(names(in_bed.ls))
  
  if (!is.null(sets)){
    sets <- names(in_bed.ls)
  }
  
  if (is.null(in_bed.ls)) {
    print("all lists are empty")
  } else {
    u.p <- upset(fromList(in_bed.ls), 
                 #nintersects = NA, 
                 nsets = length(in_bed.ls), 
                 sets=sets, order.by = "freq", number.angles = 45, 
                 text.scale = text.scale, keep.order = T)
    
    
    pdf(f.name, w,h, useDingbats = F)
    print(u.p)
    dev.off()
  }
  
  
  
  
  #print(in.ls)
}

