columns.to.coMut = function(maf, variant.categories=tcga.mapping(), type=c("mutation", "annotation"), sample.list=NULL, gene.list=NULL, gene.order.by="frequency", gene.order=NULL, sample.order.by="mutation.type", sample.order=NULL, sample.column="Tumor_Sample_Barcode", variant.column="Variant_Classification", gene.column="Hugo_Symbol", custom.sample.grouping=NULL, na.samples=NULL, na.genes=NULL, verbose=TRUE) {
 
  type = match.arg(type)
 
  ## Pull relevant columns from maf file
  samples = as.character(maf[,sample.column])
  genes = as.character(maf[,gene.column])
  variants = as.character(maf[,variant.column])
  
  ## Perform error check on gene/sample list
  sample.list = check.samples(samples, sample.list, verbose=verbose)
  maf = subset(maf, maf[,sample.column] %in% sample.list)
  gene.list = check.genes(genes, gene.list, verbose=verbose)
  
  ## Generate summary counts for mutations by sample and by gene from entire maf
  samples = factor(as.character(maf[,sample.column]), levels = sample.list)
  genes = factor(as.character(maf[,gene.column]), levels = gene.list)
  variants = as.character(maf[,variant.column])
    
  total.variants.samples.genes.full = xtabs(~ variants + samples + genes)
  
  ## Perform category name substitution  							
  variant.id = rep(0, length(variants))
  for(i in 1:nrow(variant.categories)) {
    ind = grepl(variant.categories[i,"maf.classification"], variants)
    variant.id[ind] = variant.categories[i,"id"]
  }
  variant.id = as.integer(variant.id)
  
  ## If there are other categories in the maf that were not in the variant.categories variable, they will be listed
  variants.extra = unique(variants[variant.id == 0])
  if(length(variants.extra) > 0) {
    if(verbose == TRUE) { cat("The following variant categories were in the maf but not in the variant.categories variable:", "\n", paste(variants.extra, collapse="\n"), "\nThese will not be considered in the final matrix\n") }
  }

  ## Subset to selected genes, samples, and variants
  ind = maf[,sample.column] %in% sample.list & maf[,gene.column] %in% gene.list & variant.id != 0
  samples = factor(as.character(maf[ind,sample.column]), levels = sample.list)
  genes = factor(as.character(maf[ind,gene.column]), levels = gene.list)
  variants = as.character(maf[ind,variant.column])
  variant.id = variant.id[ind]
  
  total.variants.samples.genes = xtabs(~ variants + samples + genes)
  
  variant.agg <- aggregate(formula = variant.id ~ genes + samples, FUN = max, na.rm=TRUE) 
  variant.matrix <- as.matrix(xtabs(formula = variant.id ~ genes + samples, data=variant.agg, drop.unused.levels = FALSE))
  class(variant.matrix) = "matrix"

  ## Fill in samples with missing data (i.e. NA samples)
  if(!is.null(na.samples)) {
    variant.matrix = cbind(variant.matrix, matrix(NA, nrow = nrow(variant.matrix), ncol = length(na.samples), dimnames = list(rownames(variant.matrix), na.samples)))
  }                      
  if(!is.null(na.genes)) {
    variant.matrix = rbind(variant.matrix, matrix(NA, nrow = length(na.genes), ncol = ncol(variant.matrix), dimnames = list(na.genes, colnames(variant.matrix))))
  }                      

  ## Order rows 
  if(nrow(variant.matrix) > 1) { 
    row.o = sort.coMut(variant.matrix, ordering=gene.order.by, names=gene.order, index=gene.order, margin="rows")
  } else {
    row.o = 1
  }

  ## Order columns
  if(ncol(variant.matrix) > 1) {
    col.o = sort.coMut(variant.matrix, ordering=sample.order.by, names=sample.order, index=sample.order, margin="columns")
  } else {
    col.o = 1
  }

  variant.matrix = variant.matrix[row.o,col.o,drop=FALSE]


  ## Add flag for plotting each row. Can be set to FALSE afterwards and the plotting function will ignore
  keep.in.plot = rep(TRUE, nrow(variant.matrix))
  names(keep.in.plot) = rownames(variant.matrix)

  ## Add flag for using each row in matrix ordering. Can be set to FALSE afterwards and the ordering function will ignore
  keep.in.ordering = rep(TRUE, nrow(variant.matrix))
  names(keep.in.ordering) = rownames(variant.matrix)

  return(list(matrix=variant.matrix, type=type, mapping=variant.categories, counts=total.variants.samples.genes, full.counts=total.variants.samples.genes.full, keep=keep.in.ordering, keep.in.plot))
}



matrix.to.coMut = function(mat, type=c("mutation", "annotation"), variant.categories=tcga.mapping(), sample.list=NULL, gene.list=NULL, gene.order.by="frequency", gene.order=NULL, sample.order.by="mutation.type", sample.order=NULL, na.samples=NULL, na.genes=NULL, verbose=TRUE) {
 
  type = match.arg(type)
  
  ## Pull relevant columns from maf file
  samples = colnames(mat)
  genes = rownames(mat)
  
  ## Set up sample list
  s1 =  setdiff(sample.list, samples)
  s2 =  setdiff(samples, sample.list)
  if(!is.null(sample.list) & length(s1) > 0) {
    if(verbose == TRUE) { cat("The following samples were not columns in the matrix:", "\n", paste(s1, collapse="\n"), "\n\nThey will be given zero counts across all genes\n") }
    new.mat = matrix(0, nrow=nrow(mat), ncol=length(s1), dimnames=list(rownames(mat), s1))
    mat = cbind(mat, new.mat)
  } else if (!is.null(sample.list) & length(s2) > 0) {
    if(verbose == TRUE) { cat("The following samples were columns in the matrix but were not in sample.list:", "\n", paste(s2, collapse="\n"), "\n\nThey will be excluded from the matrix and variant counts.\n") }
    mat = mat[,colnames(mat) %in% sample.list]
  } else {
    sample.list = unique(samples)
  }

  ## Set up gene list
  s1 =  setdiff(gene.list, genes)
  s2 =  setdiff(genes, gene.list)
  if(!is.null(gene.list) & length(s1) > 0) {
    if(verbose == TRUE) { cat("The following genes were not columns in the matrix:", "\n", paste(s1, collapse="\n"), "\n\nThey will be given zero counts across all samples\n") }
    new.mat = matrix(0, nrow=length(s1), ncol=ncol(mat), dimnames=list(s1, colnames(mat)))
    mat = rbind(mat, new.mat)
  } else if (!is.null(gene.list) & length(s2) > 0) {
    if(verbose == TRUE) { cat("The following genes were rows in the matrix but were not in gene.list:", "\n", paste(s2, collapse="\n"), "\n\nThey will be excluded from the matrix and variant counts.\n") }
    mat = mat[rownames(mat) %in% gene.list,]
  } else {
    gene.list = unique(genes)
  }
  
  
  ## Get list of all unique values in matrix
  mat.values = setdiff(unique(c(mat)), variant.categories[,2])
  ## If there are other categories in the maf that were not in the variant.categories variable, they will be listed
  variants.extra = setdiff(mat.values, variant.categories[,2])
  if(length(variants.extra) > 0) {
    if(verbose == TRUE) { cat("The following variant categories were in the matrix but not in the variant.categories variable:", "\n", paste(variants.extra, collapse="\n"), "\nThese will set to 0 in final matrix\n") }
    mat[mat %in% variants.extra] = 0
  }

  ## Perform category name substitution  							
  variant.matrix = matrix(0, nrow=nrow(mat), ncol=ncol(mat), dimnames=list(gene.list, sample.list))
  for(i in 1:nrow(variant.categories)) {
    ind = grepl(variant.categories[i,"maf.classification"], mat)
    variant.matrix[ind] = variant.categories[i,"id"]
  }

  
  ## Generate summary counts for mutations by sample and by gene from entire maf
  variant.matrix.melt = melt(variant.matrix)
  colnames(variant.matrix.melt) = c("genes", "samples", "variants")
  total.variants.samples.genes = xtabs(~ variants + samples + genes, data=variant.matrix.melt) 
  
  ## Fill in samples with missing data (i.e. NA samples)
  if(!is.null(na.samples)) {
    variant.matrix = cbind(variant.matrix, matrix(NA, nrow = nrow(variant.matrix), ncol = length(na.samples), dimnames = list(rownames(variant.matrix), na.samples)))
  }                      
  if(!is.null(na.genes)) {
    variant.matrix = rbind(variant.matrix, matrix(NA, nrow = length(na.genes), ncol = ncol(variant.matrix), dimnames = list(na.genes, colnames(variant.matrix))))
  }                      

  ## Order rows 
  if(nrow(variant.matrix) > 1) { 
    row.o = sort.coMut(variant.matrix, ordering=gene.order.by, names=gene.order, index=gene.order, margin="rows")
  } else {
    row.o = 1
  } 

  ## Order columns
  if(ncol(variant.matrix) > 1) {
    col.o = sort.coMut(variant.matrix, ordering=sample.order.by, names=sample.order, index=sample.order, margin="columns")
  } else {
    col.o = 1
  }

  variant.matrix = variant.matrix[row.o, col.o, drop=FALSE]
  
  ## Add flag for plotting each row. Can be set to FALSE afterwards and the plotting function will ignore
  keep.in.plot = rep(TRUE, nrow(variant.matrix))
  names(keep.in.plot) = rownames(variant.matrix)

  ## Add flag for using each row in matrix ordering. Can be set to FALSE afterwards and the ordering function will ignore
  keep.in.ordering = rep(TRUE, nrow(variant.matrix))
  names(keep.in.ordering) = rownames(variant.matrix)

  return(list(matrix=variant.matrix, type=type, mapping=variant.categories, counts=total.variants.samples.genes, keep=keep.in.ordering, keep.in.plot))
}








tabler = function(tabler.obj, plot.samples="all", plot.samples.index=NULL,
  xlab = "", ylab = "", 
  horizontal.border = TRUE,
  vertical.border = TRUE,
  horizontal.border.args = list(lwd=0.25, col="grey90"),
  vertical.border.args = list(lwd=0.25, col="grey90"),
  box = TRUE,
  box.args = list(col="black", lwd=1),
  row.label = TRUE,
  column.label = TRUE,
  na.col = "grey90",
  bg.col = "white",
  row.freq = TRUE,
  row.label.args = list(),
  column.label.args = list(),
  row.freq.args = list(),
  image.args = list(),
  column.label.side = 'bottom', row.label.side = "left", row.freq.side = "right",
  y.adj = 0, # adjustment for row labels
  x.adj = 0, # adjustment for col labels
  bar.right.args= list(),
  bar.left.args = list(),
  cexCol=NULL,
  cexRow=NULL,
  ...)
  {

    mat = tabler.obj$matrix
    tabler.mapping = tabler.obj$mapping

    if(nrow(mat) > 1) {
      mat = mat[nrow(mat):1,]
    }
    row.freq.percent = round((rowSums(mat>0, na.rm=TRUE) / rowSums(!is.na(mat)))*100)
    row.freq.label = paste(row.freq.percent, "%", sep="")
    row.freq.label[row.freq.percent > 0 & row.freq.percent < 1] = "<1%"

    ## Remove unaltered samples if requested
    if(plot.samples == "mutated") {
      one.mut.present = colSums(mat == 0 | is.na(mat) | is.infinite(mat), na.rm=TRUE) != nrow(mat)
      mat = mat[,one.mut.present]

    } else if (plot.samples == "selected" & !is.null(plot.samples.index)) {
      mat = mat[,plot.samples.index]
    }


    tabler.mapping.id = as.integer(tabler.mapping$id)
    tabler.mapping.color = as.character(tabler.mapping$bg.col)
 
    ## Add color for zero if it isn't already in there
    if(sum(tabler.mapping.id == 0) == 0) {
      tabler.mapping.id = c(tabler.mapping.id, 0)
      tabler.mapping.color = c(tabler.mapping.color, bg.col)
    }
    ## Add color for NA samples if they exist in the matrix
    if(sum(is.na(mat)) > 0) {
      na.id = min(tabler.mapping.id)-1
     
      tabler.mapping.id = c(tabler.mapping.id, na.id)
      tabler.mapping.color = c(tabler.mapping.color, na.col)
      
      mat[is.na(mat)] = na.id
    }

    ## Order mapping ids and colors for use in image()
    ind = order(tabler.mapping.id)
    tabler.mapping.id = tabler.mapping.id[ind]
    tabler.mapping.color = tabler.mapping.color[ind]
    
    ## Remove values that were in the mapping scheme but not in the actual matrix
    id.in.matrix = tabler.mapping.id %in% c(mat)
    tabler.mapping.id = tabler.mapping.id[id.in.matrix]
    tabler.mapping.color = tabler.mapping.color[id.in.matrix]
    
    ## Rescale IDs such that values in matrix are increasing by one    
    ## Needed for image(), otherwise it behaves strange when trying to color
    new.mat = mat
    new.ids = 1:length(tabler.mapping.id)
    for(i in 1:length(new.ids)) {
      new.mat[mat == tabler.mapping.id[i]] = new.ids[i] 
    }
 

    ## Make matrix plot    
    image.required.args = list(1:ncol(new.mat), 1:nrow(new.mat), t(new.mat), col=tabler.mapping.color, axes=FALSE, xlab=xlab, ylab=ylab, zlim=c(min(new.mat, na.rm=TRUE)-1, max(new.mat, na.rm=TRUE)+1))
    image.args = c(image.required.args, image.args)
    do.call.args("image", image.args)
    
    ## Add Vertical borders
    if(vertical.border == TRUE) {
      vertical.segs = list(.5 + 1:(ncol(new.mat)-1), .5, .5 + 1:(ncol(new.mat)-1), .5 + nrow(new.mat))
      do.call.args("segments", vertical.segs, vertical.border.args)
    }
    
    ## Add Horizontal borders
    if(horizontal.border == TRUE) {
      horizontal.segs = list(.5, .5 + 1:(nrow(new.mat)-1), .5 + ncol(new.mat), .5 + 1:(nrow(new.mat)-1))
      do.call.args("segments", horizontal.segs, horizontal.border.args)
    }
 

    
    ## Add column labels
    if(column.label == TRUE) {
      if (column.label.side == "top") {
        x.side = 3
      } else {
        x.side = 1
      } 
      
      ## Get cex to fit colnames in margin
      if(is.null(cexCol)) {
        cexCol <- compute.cex(width = (par("mar")[1] - par("mgp")[2])* mar.to.mai() * 0.8, height = par("pin")[1] / ncol(new.mat), label = colnames(new.mat))            
      }
      column.label.required.args = list(side=x.side, line=0.2, at = 1:ncol(new.mat), labels = colnames(new.mat), tick = FALSE, las = 2, cex.axis = cexCol)
      do.call.args("axis", column.label.required.args, column.label.args)
    }  
    
    ## Add row labels
    if(row.label == TRUE) {
      if (row.label.side == "left") {
        y.side = 2
      } else {
        y.side = 4
      } 
      
      ## Get cex to fit rownames in margin
      if(is.null(cexRow)) {
        cexRow <- compute.cex(width = (par("mar")[2] - par("mgp")[2]) * mar.to.mai() * 0.8, height = par("pin")[2] / nrow(new.mat), label = rownames(new.mat))
      }    
      row.label.required.args = list(side=y.side, line=0.2, at = 1:nrow(new.mat), labels = rownames(mat), tick = FALSE, las = 1, cex.axis = cexRow)
      do.call.args("axis", c(row.label.required.args, row.label.args))
    }
    
    ## Add row frequecies
    if(row.freq == TRUE) {
          
      if (row.freq.side == "left") {
        y.side = 2
      } else {
        y.side = 4
      } 

      row.freq.required.args = list(side=y.side, line=0.2, at = 1:nrow(new.mat), labels = row.freq.label, tick = FALSE, las = 1, cex.axis = cexRow)
      do.call.args("axis", c(row.freq.required.args, row.freq.args)) 
    }  
    
    ## Plot box around matrix    
    if (box == TRUE) {
      do.call.args("box", c(list(which="plot"), box.args))
    }
   
}






coMut.with.bars = function(coMut.obj, left.bar=c("counts", "frequency", "none", "custom"), left.bar.args=list(), right.bar=c("frequency","counts","none","custom"), right.bar.args=list(), tabler.args=list(row.freq=FALSE), mar=par("mar"), close.screens=TRUE) {

  left.bar = match.arg(left.bar)
  right.bar = match.arg(right.bar)
  oma = c(0,0,0,0)
  mgp = c(0,0,0)
  ## Set up screen dimensions
  screen.dims = rbind(c(0,0.1,0,1), c(0.10,0.70,0,1), c(0.70,0.90,0,1), c(0.90,1,0,1))
  local.screens = split.screen(screen.dims)

  ## close screens on exit if requested
  on.exit( {
    if(close.screens == TRUE) {
      close.screen(local.screens)
  }})
  
  #### Plot left bar plot
  screen(local.screens[1])
  temp.mar = mar
  temp.mar[4] = 0
  par(mar=temp.mar, oma=oma, mgp=mgp)
#  if(left.bar != "none") {
#    left.bar.required.args = list(counts=coMut.obj$counts, horiz=TRUE, gene.order=rev(rownames(coMut.obj$matrix)), draw.gene.labels=FALSE)
#    do.call.args("gene.barplot", left.bar.required.args, left.bar.args) 
#  }  
  
  ### Plot matrix of mutations  
  screen(local.screens[2])
  plot.row.freq = FALSE
  if(!is.null(tabler.args[["row.freq"]])) {
    if(coMut.obj$type == "mutation" & tabler.args[["row.freq"]] == TRUE) {
      plot.row.freq = TRUE  
    }
  }
  
  temp.mar = mar
  temp.mar[2] = 4.1
  if(plot.row.freq == TRUE) {
    temp.mar[4] = 2.1
  } else {
    temp.mar[4] = 0
  }
  par(mar=temp.mar, oma=oma, mgp=mgp)
  do.call.args("tabler", list(tabler.obj=coMut.obj, row.freq=plot.row.freq), tabler.args)

  #### Plot right bar plot
  screen(local.screens[3])
  temp.mar = mar
  temp.mar[2] = 0.5
  
  if(!is.null(right.bar.args[["draw.freq"]])) {
    if(right.bar.args[["draw.freq"]] == TRUE) {
      if(is.null(right.bar.args[["vertical.group"]])) {
        temp.mar[4] = 1
      } else {
        temp.mar[4] = length(unique(right.bar.args[["vertical.group"]]))
      }
    }
  } else {
    temp.mar[4] = 0  
  }  
  temp.mgp = c(0.5,0.1,0.02)

  par(mar=temp.mar, oma=oma, mgp=temp.mgp)    
  if(right.bar != "none" & coMut.obj$type == "mutation") {
    right.bar.required.args = list(counts=coMut.obj$counts, horiz=TRUE, gene.order=rev(rownames(coMut.obj$matrix)), draw.gene.labels=FALSE)
    do.call.args("gene.barplot", right.bar.required.args, right.bar.args) 
  }  
  
  ## Add legend on the last screen
  screen(local.screens[4])
  temp.mar = mar
  temp.mar[2] = 0
  par(mar=temp.mar, oma=oma, mgp=mgp, xpd=NA)

  #cex.legend = compute.cex(width = par("pin")[1]-0.3, height=par("pin")[2] / nrow(coMut.obj$mapping), label = coMut.obj$mapping[,"final.classification"])
  cex.legend = 0.5
  plot(0, ylab="", xlab="", type="n", xaxt="n", yaxt="n", bty="n", yaxs="i", xaxs="i", ann=FALSE)

  legend.required.args = list(x="left", legend=as.character(coMut.obj$mapping[,"final.classification"]), fill=as.character(coMut.obj$mapping[,"bg.col"]), bty="n", pt.cex=1, cex=cex.legend, y.intersp=0.75)
  do.call.args("legend", legend.required.args)
  
}





multi.coMut = function(table.list, sample.order="mutation.type", close.screens=TRUE, between.table.space = 0.025, row.space = 0.025, mar=c(1,3,1,1), bg.col="white", na.col="grey90", labCol=TRUE, plot.samples="all", left.bar.type=c("none", "counts", "frequency", "custom"), right.bar.type=c("counts", "frequency", "custom", "none"), right.bar.beside=NULL, right.bar.stack=NULL, right.bar.args=list(), tabler.args=list()) {
  
  
  left.bar.type = match.arg(left.bar.type)
  right.bar.type = match.arg(right.bar.type)
    
  ntable = length(table.list)
  
  ## Check to verify that column names match for each table in the list  
  table.colnames = colnames(table.list[[1]][["matrix"]])
  if(ntable > 1) {
    for(i in 2:ntable) {
      temp.table = table.list[[i]][["matrix"]]
      if(!is.overlapping(table.colnames, colnames(temp.table))) {
        stop(paste("Column names in table", i, "do not match those from table 1. Cannot align columns tables properly\n", sep=" "))
      }
    }
  }
  
  ## Create meta-matrix for ordering, start with last one
  meta.matrix = table.list[[ntable]][["matrix"]]
  rows.for.ordering = table.list[[ntable]][["keep"]]
  meta.matrix = meta.matrix[,table.colnames,drop=FALSE][rows.for.ordering,]
  
  if(ntable > 1) {
    for(i in (ntable-1):1) {

      ## From the last table to the first, increase or decrease the ids so they are non-overlapping 
      ## with all other tables. Then concatenate into one super table for sorting
      meta.matrix.max = max(meta.matrix, na.rm=TRUE)
      meta.matrix.min = min(meta.matrix, na.rm=TRUE)

      rows.for.ordering = table.list[[i]][["keep"]]
      matrix.to.add = table.list[[i]][["matrix"]][,table.colnames,drop=FALSE][rows.for.ordering,,drop=FALSE]

      
      ## Shift baseline value in matrix to zero if no zeros are present for each row
      ## Mainly for annotation like matrices that do not have a baseline (i.e. unmutated) state
      matrix.to.add.min = apply(matrix.to.add, 1, min, na.rm=TRUE)
      ind = matrix.to.add.min > 0
      matrix.to.add[ind,] = sweep(matrix.to.add[ind,,drop=FALSE], 1, matrix.to.add.min[ind], "-")

      ## Shift values to be greater/less than the values currently in the meta.matrix
      ## This will give priority to tables earlier in the table.list for sorting
      ind = matrix.to.add < 0 & !is.na(matrix.to.add)
      matrix.to.add[ind] = matrix.to.add[ind] + meta.matrix.min
      
      ind = matrix.to.add > 0 & !is.na(matrix.to.add)
      matrix.to.add[ind] = matrix.to.add[ind] + meta.matrix.max
      
      meta.matrix = rbind(matrix.to.add[,table.colnames,drop=FALSE], meta.matrix)
    }
  } 

  ## Determine which samples do and do not have at least one mutation
  if(plot.samples == "all") {
    samples.to.show = rep(TRUE, ncol(meta.matrix))
  } else {
    samples.to.show = rep(FALSE, ncol(meta.matrix))
    for(i in 1:ntable) {
      is.mut = apply(table.list[[i]][["matrix"]][,table.colnames, drop=FALSE] != 0, 2, sum, na.rm=TRUE) > 0
      samples.to.show[is.mut] = TRUE  
    }
  } 

  ## Perform ordering of samples across all tablers
  global.sample.order = sort.coMut(meta.matrix, margin="columns", ordering=sample.order)
  global.sample.names = colnames(meta.matrix)[global.sample.order]

  samples.to.show.order = samples.to.show[global.sample.order]
  
  
  ## Plot tables one by one
  table.nrows = unlist(lapply(1:ntable, function(i) { nrow(table.list[[i]][["matrix"]]) }))
  total.nrow = sum(table.nrows)

  mar.fig.fraction = mar.to.fig.fraction(mar)
  table.spacing = row.space*table.nrows
  total.table.space = sum(table.spacing)
  total.between.table.space = between.table.space*(ntable-1)
  total.space = total.table.space + total.between.table.space + mar.fig.fraction[1] + mar.fig.fraction[3] 
  
  ## Shrink row height if the total number of rows exceed the amount available in plot
  if(total.space > 1) {    
    usable.space = 1 - mar.fig.fraction[1] - mar.fig.fraction[3] 
    inflation.factor = (total.table.space + total.between.table.space) / usable.space
    row.space = row.space/inflation.factor
    between.table.space = between.table.space/inflation.factor
    table.spacing = row.space*table.nrows
    warning("Number of rows exceeded space available in plot. 'row.space' and 'between.table.space' have been lowered to accomodate extra rows.")
  }
  
  
  ## Add space to bottom and top for margins
  table.spacing[1] = table.spacing[1] + mar.fig.fraction[3]
  table.spacing[ntable] = table.spacing[ntable] + mar.fig.fraction[1]

  ## Set up split.screen matrix
  between.table.shift = c(cumsum(rep(between.table.space, ntable))) 
  ss.top = c(1, 1-cumsum(table.spacing) - between.table.shift)
  ss.bottom = c(1-cumsum(table.spacing),0) - c(0, between.table.shift)
  ss.left = 0  ## par(mar) will be used for now to move the sides of the plots in
  ss.right = 1

  ss.matrix = round(cbind(ss.left, ss.right, ss.bottom, ss.top), 3)
  
  ## Remove dummy row at the end and ensure it is a matrix (for the case when only one element in table.list)
  ss.matrix = matrix(ss.matrix[-nrow(ss.matrix),], byrow=FALSE, ncol=4) 
  if(sum(ss.matrix < 0) > 0) { stop("split.screen matrix has negative values") }


  ## Check to see if ss.matrix values are outside [0,1] interval and adjust if necessary
  #ss.matrix = rescale.screen(ss.matrix)
  
  
  ## Get gene name cex to use for all coMuts
  all.names = c()
  for(i in 1:ntable) {
    all.names = c(all.names, rownames(table.list[[i]]$matrix))
  }
  gene.name.cex = compute.cex(width=mar.fig.fraction[2]*par("fin")[1], height=0.8*row.space*par("fin")[2], label=all.names, cex.lim=c(0.01, 1))
  
  ## Get groups for right.bar beside groups 
  right.bar.beside.group = NULL
#  bar.col = NULL
  if(!is.null(right.bar.beside)) {
    ind = as.numeric(right.bar.beside[1])

    temp.matrix = table.list[[ind]]$matrix
    temp.mapping = table.list[[ind]]$mapping

	right.bar.beside.group = temp.matrix[right.bar.beside[2],]
    matrix.values = unique(c(right.bar.beside.group))
    temp.mapping = subset(temp.mapping, temp.mapping[,1] %in% matrix.values)
    right.bar.beside.group = mapvalues(right.bar.beside.group, as.character(temp.mapping[,1]), as.character(temp.mapping[,3]), warn_missing=FALSE)
 #   bar.col = rbind(bar.col, as.character(temp.mapping$bg.col))
    names(right.bar.beside.group) = colnames(temp.matrix)
    right.bar.beside = "samples"
  } 

  ## Get groups for right.bar stack groups
  right.bar.stack.group = NULL  
  if(!is.null(right.bar.stack)) {
    ind = as.numeric(right.bar.stack[1])
    temp.matrix = table.list[[ind]]$matrix
    temp.mapping = table.list[[ind]]$mapping

	right.bar.stack.group = temp.matrix[right.bar.stack[2],]
    matrix.values = unique(c(right.bar.stack.group))
    temp.mapping = subset(temp.mapping, temp.mapping[,1] %in% matrix.values)
        
    right.bar.stack.group = mapvalues(right.bar.stack.group, as.character(temp.mapping[,1]), as.character(temp.mapping[,3]), warn_missing=FALSE)
    
#    if(is.null(bar.col)) {
#      bar.col = cbind(bar.col, as.character(temp.mapping$bg.col))
#    } else {
#      bar.col = rbind(bar.col, as.character(temp.mapping$bg.col))
#    }
    names(right.bar.stack.group) = colnames(temp.matrix)
    right.bar.stack = "samples"
  }  

  ## Get maximum bar to use for right bar
  max.count = NULL
  if(right.bar.type != "none") {
    l = list(type=right.bar.type, horizontal.group.by=right.bar.beside, horizontal.group=right.bar.beside.group, vertical.group.by=right.bar.stack, vertical.group=right.bar.stack.group)
    max.count = get.max.bar.from.table.list(table.list, l)
  } 

  ## Plot each coMut by using split screen
  screen.index = split.screen(ss.matrix)
  
  ## close screens on exit if requested
  on.exit( {
    if(close.screens == TRUE) {
      close.screen(screen.index)
  }})

  plot.bottom.labels = FALSE
  for(i in 1:ntable) {
    
    if(i == ntable & labCol == TRUE) {
      plot.bottom.labels = TRUE
      plot.bottom.axis = TRUE
      temp.mar = mar
      temp.mar[3] = 0
    } else if (i == ntable & labCol == FALSE) {
      plot.bottom.labels = FALSE
      plot.bottom.axis = TRUE
      temp.mar = mar
      temp.mar[3] = 0
    } else if (i == 1) {
      plot.bottom.axis = FALSE
      temp.mar = mar
      temp.mar[1] = 0
    } else {
      plot.bottom.axis = FALSE
      temp.mar[c(1,3)] = 0
    }
  
    screen(screen.index[i])
    new.tabler = table.list[[i]]
    new.tabler$matrix = new.tabler$matrix[,global.sample.names,drop=FALSE]

    ## Plot table
    coMut.with.bars(new.tabler, mar=temp.mar,
    	tabler.args=c(tabler.args, list(cexRow=gene.name.cex, column.label=plot.bottom.labels, plot.samples="selected", plot.samples.index=samples.to.show.order, bg.col=bg.col, na.col=na.col)),
    	right.bar.args=c(right.bar.args, list(type=right.bar.type, draw.axis=plot.bottom.axis, draw.gene.labels=FALSE, draw.axis.title=plot.bottom.axis, axis.lim=c(0,max.count), horizontal.group=right.bar.beside.group, horizontal.group.by=right.bar.beside, vertical.group.by=right.bar.stack, vertical.group=right.bar.stack.group)))

  }  

  invisible(list(meta.matrix, between.table.space, row.space))
}





gene.barplot = function(counts, col=NULL, gene.order=NULL, type=c("counts", "frequency"), vertical.group.by=c("none", "variants", "samples"), vertical.group=NULL, horizontal.group.by=c("none", "variants", "samples"), horizontal.group=NULL, horiz=FALSE, reverse=FALSE, axis.lim=NULL, space=0.1, draw.axis=TRUE, draw.gene.labels=TRUE, draw.axis.title=TRUE, draw.freq=FALSE, plot=TRUE) {

  ## Takes a 3-dimensional array of variant by sample by gene counts and generates a barplot according to the supplied grouping structure
  ## Written by: Josh Campbell
  ## 2-15-2016
  
  ## Parameters: ##########################################################################################
  #
  # counts				3 dimentional array corresponding to variants, samples, gene counts
  # col					A matrix with nrows equal to the number of vertical groups and ncols equal to the number of
  #						horizontal groups. Each entry in the matrix will correspond to the color of that subgroup
  #						within a gene.
  # gene.order			List of gene names to use for ordering the barplot. Otherwise, it stays in the same order as in the array
  # type				If set to 'counts', then the total sum of mutations will be displayed for each group within each gene
  #						If set to 'frequency', then the frequency of patients mutated for each group will be displayed. Only valid
  #						if 'vertical.group.by' or 'horizontal.group.by' is 'none' or 'samples'. Not valid with 'variant' grouping.
  # vertical.group		A vector with values corresponding to the group to re-map each variant/sample. Names of the vector should correspond 
  #						to the actual variant/sample to be re-mapped. The stacked part of the barplot will be grouping according to this variable.
  #	horizontal.group	Same as 'vertical.group.by' but for grouping the side-by-side bars
  # vertical.group.by	One of 'none', 'variants', 'samples'. Variants or samples will be grouped on the stacked bar accordingly
  # horizontal.group.by	One of 'none', 'variants', 'samples'. Variants or samples will be grouped on the side-by-side bar accordingly
  # horiz				Whether to rotate the barplot by 90 degrees. Similar to 'horiz=TRUE' in base barplot. 
  # reverse				Whether to refect the barplot
  # axis.lim			A vector with 2 entries. Defaults to zero and just above the maximum value in the plot
  # space				Space to have between bars for different genes
  # draw.axis			Draw axis. Default is TRUE.
  # draw.gene.labels			Draw gene names. Default is TRUE.
  # draw.axis.title			Draw axis label. Default is TRUE.
  # plot				Whether to draw the plot
  #########################################################################################################
  
  ## Values: ##############################################################################################
  #
  #	bars.by.group		If requested, the new array counts/frequencies per group will be requested
  #	counts				The matrix of counts/frequences used for plotting. Directly correpsonds value in barplot
  #	rect.details		The values used in the 'rect' function to draw the bars
  #
  #########################################################################################################
  require(plyr)
  require(reshape2)

  type = match.arg(type)
  vertical.group.by = match.arg(vertical.group.by)
  horizontal.group.by = match.arg(horizontal.group.by)

  ## Additional error checking
  if(type == "frequency" & (vertical.group.by == "variants" | horizontal.group.by == "variants")) {
    warning("Using type 'frequency' with 'variant' grouping may cause unexpected results as two different variants from the same patient will be considered separately. Therefore this is disabled.")
  }
  
  if(is.null(gene.order)) {
    gene.order = dimnames(counts)$genes
  } else if (!is.overlapping(gene.order, dimnames(counts)$genes)) {
    stop("gene.order needs to be a vector of genes names perfectly overlapping with dimnames(counts)[3]")
  } else {
    counts = counts[,,gene.order,drop=FALSE]
  }

  n.variants = dim(counts)[1]
  n.samples = dim(counts)[2]
  n.gene = dim(counts)[3]

  ## Decompose array into data.frame based on counts, or summing out the different variants to get at frequency
  if(type == "counts" ) {
    counts.melt = melt(counts)  
  } else if (type == "frequency") {
    counts.melt = data.frame(variants="Any", melt(apply(counts, c("genes", "samples"), sum) > 0))
  }

  ## Create a vectors correpsonding to the horiztonal and vertical grouping variables  
  if(vertical.group.by != "none" & !is.null(vertical.group)) {
    vertical.group.ix = as.factor(mapvalues(counts.melt[,vertical.group.by], names(vertical.group), vertical.group, warn_missing=FALSE))
  } else if (vertical.group.by != "none") {
    vertical.group.ix = as.factor(counts.melt[,vertical.group.by])
  } else {
    vertical.group.ix = as.factor(rep(1, nrow(counts.melt)))
  }  

  if(horizontal.group.by != "none" & !is.null(horizontal.group)) {  
    horizontal.group.ix = as.factor(mapvalues(counts.melt[,horizontal.group.by], names(horizontal.group), horizontal.group, warn_missing=FALSE))
  } else if (horizontal.group.by != "none") {
    horizontal.group.ix = as.factor(counts.melt[,horizontal.group.by])
  } else {
    horizontal.group.ix = as.factor(rep(1, nrow(counts.melt)))
  }  
  
  n.horizontal.groups = length(levels(horizontal.group.ix))
  n.vertical.groups = length(levels(vertical.group.ix))
 

  ## Add new group indexes to data frame
  counts.melt = data.frame(counts.melt, horizontal.group.ix, vertical.group.ix, stringsAsFactors=FALSE)

  ## Recount according to the new grouping scheme 
  if(type == "counts") {
  
    counts.by.group = xtabs(value ~ genes + horizontal.group.ix + vertical.group.ix, data=counts.melt)
    axis.lab = "Mutation counts"
    
  } else if (type == "frequency" & vertical.group.by == "samples" & horizontal.group.by == "none") {

    counts.by.group = xtabs(value ~ genes + horizontal.group.ix + vertical.group.ix, data=counts.melt) 
    n.samp = table(vertical.group)
    n.samp = n.samp[dimnames(counts.by.group)$vertical.group.ix]
    counts.by.group = sweep(counts.by.group, 3, n.samp, "/")
    axis.lab = "Mutation frequency"
    
  } else if (type == "frequency" & horizontal.group.by == "samples" & vertical.group.by == "none") {

    counts.by.group = xtabs(value ~ genes + horizontal.group.ix + vertical.group.ix, data=counts.melt) 
    n.samp = table(horizontal.group)
    n.samp = n.samp[dimnames(counts.by.group)$horizontal.group.ix]
    counts.by.group = sweep(counts.by.group, 2, n.samp, "/") 
    axis.lab = "Mutation frequency"

  } else if (type == "frequency" & horizontal.group.by == "none" & vertical.group.by == "none") {

    n.samp = dim(counts)[2]
    counts.by.group = xtabs(value ~ genes + horizontal.group.ix + vertical.group.ix, data=counts.melt) 
    counts.by.group = counts.by.group / n.samp
    axis.lab = "Mutation frequency"

  } else {
    stop("type='frequency' cannot be used with 'variants' or with both horizontal.group.by and vertical.group.by set to 'samples'")
  }

 
  ## Combine counts so that all values for the same gene are next to each other
  final = t(apply(counts.by.group, "vertical.group.ix", cbind))
  final = rbind(0, final)
  final.cs = sapply(as.data.frame(final), FUN=cumsum)

  ## Add names
  rownames(final) = c("Offset", dimnames(counts.by.group)[["vertical.group.ix"]])
  cn.1 = rep(dimnames(counts.by.group)[["genes"]], n.horizontal.groups)
  cn.2 = rep(dimnames(counts.by.group)[["horizontal.group.ix"]], each=n.gene)
  colnames(final) = paste(cn.1, cn.2, sep="_")


  ## Determine axis limits
  if(is.null(axis.lim)) {
    m = max(final.cs)
    if(m < 1) {
      axis.max = round_any(m*100, 10, ceiling)/100
    } else if (m < 10) {
      axis.max = round_any(m, 1, ceiling)
    } else {
      axis.max = round_any(m, 10, ceiling)
    }
    axis.min = 0
  } else {
    axis.min = axis.lim[1]
    axis.max = axis.lim[2]
  }

  ## If the plot needs to be reflected, this can be done by taking the negative of each count
  if(reverse == TRUE) {
    final.cs = -final.cs
    temp = axis.max
    axis.max = -axis.min
    axis.min = -temp
  }
  
  ## Get spacing on y axis
  NR = nrow(final.cs) 
  NC = ncol(final.cs)
  y.start = c(t(final.cs[-NR,]))
  y.end = c(t(final.cs[-1,]))
  y.actual.counts = c(t(final[-1,]))
  
  ## Get spacing on x axis
  x.start.base = space[1L] + 0:(n.gene-1)
  x.end.base = 1:n.gene - space[1L]
  x.adj = sapply(1:n.gene, function(i) seq(x.start.base[i], x.end.base[i], length=n.horizontal.groups+1))
  x.start = c(t(x.adj[-(n.horizontal.groups+1),]))
  x.end = c(t(x.adj[-1,]))

  ## Replicate for each 
  x.start = rep(x.start, n.vertical.groups)
  x.end = rep(x.end, n.vertical.groups)
  
  ## Replicate color for each set of bars within each gene
  if(is.null(col)) {
    ## Default set to gray colors
    temp.col=matrix(gray.colors(n.horizontal.groups*n.vertical.groups), ncol=n.horizontal.groups, nrow=n.vertical.groups, byrow=TRUE)
  } else if(!is.matrix(col) & length(col) == n.horizontal.groups & n.vertical.groups == 1) {
    ## If only horizontal groups are used
    temp.col = matrix(col, nrow=1, ncol=n.horizontal.groups, byrow=TRUE)
  } else if(!is.matrix(col) & length(col) == n.vertical.groups & n.horizontal.groups == 1) {
    ## If only vertical groups are used
    temp.col = matrix(col, ncol=1, nrow=n.vertical.groups, byrow=TRUE)
  } else if(is.matrix(col) & (ncol(col) == n.horizontal.groups & nrow(col) == n.vertical.groups)) {
    temp.col = col
  } else if(is.matrix(col) & (ncol(col) != n.horizontal.groups | nrow(col) != n.vertical.groups)) {
    stop("'col' must be a matrix with ncol equal to the number of unique horizontal groups and nrow equal to the number of unique vertical groups.")
  } else if(is.vector(col)) {
    stop("'col' must be a vector equal in length to the number of groups.")
  }
  new.col = rep(c(t(temp.col)), each=n.gene)

  ## Perform the plotting, different if horiz=TRUE or FALSE  
  plotting.matrix = NULL
  if(plot == TRUE) {
  if(horiz==TRUE) {

    plot(0, xlim=c(axis.min,axis.max), ylim=c(0,n.gene), ylab="", xlab="", type="n", xaxt="n", yaxt="n", bty="n", yaxs="i", xaxs="i", ann=FALSE)
    rect(y.start, x.start, y.end, x.end, xaxt="n", yaxt="n", col=new.col)
    plotting.matrix = data.frame(x1=y.start, y1=x.start, x2=y.end, y2=x.end, color=new.col, y.actual.counts, stringsAsFactors=FALSE)

    axis.side = 1
    fin = 2
    axp = "xaxp"    
    if(reverse == TRUE) {
      label.side=4
      freq.side=2
      freq.las=2
    } else {
      label.side=2
      freq.side=4
      freq.las=2
    }  

  } else {

    plot(0, xlim=c(0,n.gene), ylim=c(axis.min,axis.max), ylab="", xlab="", type="n", xaxt="n", yaxt="n", bty="n", xaxs="i", yaxs="i", ann=FALSE)
    rect(x.start, y.start, x.end, y.end, xaxt="n", yaxt="n", col=new.col)
    plotting.matrix = data.frame(x1=x.start, y1=y.start, x2=x.end, y2=y.end, color=new.col, y.actual.counts, stringsAsFactors=FALSE)

    axis.side = 2
    fin = 1
    axp = "yaxp"
    if(reverse == TRUE) {
      label.side=3
      freq.side=1
      freq.las=1
    } else {
      label.side=1
      freq.side=3
      freq.las=1    
    }  
  } 

  
  
  ## Plot addtional features/labels if requested
  mar = par("mar")
  mai = par("mai")
  
  if(draw.axis == TRUE) {
    cex = 0.8 * compute.cex(width = par("pin")[axis.side] , height = mai[axis.side], label = axis.lab)
    xaxp = par(axp)
    xaxp.lab = seq(xaxp[1], xaxp[2], length=xaxp[3]+1)
    s1 = seq(1, length(xaxp.lab), by=2)  ## Trick R into plotting all axis labels
    s2 = seq(2, length(xaxp.lab), by=2)
    axis(axis.side, las=1, at=xaxp.lab[s1], cex.axis=cex, labels=abs(xaxp.lab)[s1])
    axis(axis.side, las=1, at=xaxp.lab[s2], cex.axis=cex, labels=abs(xaxp.lab)[s2])
  }
  if(draw.axis.title == TRUE) {
    cex = 0.8 * compute.cex(width = par("pin")[axis.side], height = mai[axis.side], label = axis.lab)  
    mtext(axis.lab, side=axis.side, cex=cex, line=par("mgp")[1])
  }
  if(draw.gene.labels == TRUE) { 
    cex = compute.cex(width = (mar[label.side] - par("mgp")[2] - 0.5)* mar.to.mai(), height = par("pin")[fin] / n.gene, label = dimnames(counts.by.group)[[1]])
    axis(label.side, las=2, at=1:n.gene - 0.5, cex.axis=cex, labels=dimnames(counts.by.group)[[1]])
  }
  if(draw.freq == TRUE) {
    ## Make labels differently for frequency vs. counts
    if(type == "frequency") {
      y.round = as.character(round(y.actual.counts*100))
      y.lab = as.character(y.round)
      y.lab[y.round > 0 & y.round < 1] = "<1"
    } else {
      y.lab = as.character(y.actual.counts)
    }
    
    ## Get the mid point for each horizontal.bar
	mar.loc = x.start + ((x.end - x.start)/2)
    mar.line = rep(0:(n.vertical.groups-1), each=n.horizontal.groups*n.gene) + 0.5
    cex = 0.8 * compute.cex(width = ((mar[freq.side] - par("mgp")[2]) / n.vertical.groups) * mar.to.mai(), height = par("pin")[fin] / n.gene / n.horizontal.groups, label = y.lab)
    
    ## Plot for each line (vertical.group)
    for(i in unique(mar.line)) {
      ind = mar.line == i
      s1 = seq(1, length(ind), by=2)  ## Trick R into plotting all labels even if they overlap by plotting every other one 
      s2 = seq(2, length(ind), by=2)
      axis(at=mar.loc[ind][s1], labels=y.lab[ind][s1], side=freq.side, line=i, tick=FALSE, las=freq.las, hadj=0.5, padj=0.5, cex.axis=cex)
      axis(at=mar.loc[ind][s2], labels=y.lab[ind][s2], side=freq.side, line=i, tick=FALSE, las=freq.las, hadj=0.5, padj=0.5, cex.axis=cex)
    }
  }
  }

  invisible(list(bars.by.group=counts.by.group, counts=final, rect.details=plotting.matrix, lim=c(axis.min, axis.max)))
}   





sort.coMut = function(	M,
						ordering=c("frequency", "mutation.type", "frequency.mutation.type", "names", "index", "none"),
						margin=c("rows", "columns"),
						names=NULL,
						index=NULL,
						mutation.cutoff=0,
						decreasing=TRUE) {

  ## Takes a mutation matrix M and performs sorting on rows or columns
  ## based off of types of mutations, frequencies, or a combination of both
  ## Written by: Josh Campbell
  ## 2-1-2016
  
  ## Parameters: ##########################################################################################
  #
  # M				Integer matrix of mutations. Genes/features are rows. Samples are columns. Different types
  #					of mutations should be indicated by different integers.
  # ordering		One of "frequency", "frequency.mutation.type", "mutation.type", "names", or "index".
  #					frequency: order by most frequently mutated. A gene is considered mutated in a sample if it is greater than "mutation.cutoff"			
  #					mutation.type: order by mutation type. Mutations represented by higher integers are given priority
  #					names: order the margin according to the names given in the vector "names"
  #					index: order the margin according to the indicies given in the vector "index"
  #					none: no ordering is performed 
  # margin			One of "rows" or "columns"
  # names			Character vector with names to order the rows or columns. Useful if custom sorting has
  #					already been performed. Should be the same length as number of rows or columns and 
  #					overlap with the rownames or colnames of the matrix.
  # index			Integer vector with indicies to order the rows or columns. E.g. Useful if custom sorting has
  #					already been performed. Should be the same length as number of rows or columns and 
  #					should not be out of range of the rows or columns.
  #	mutation.cutoff	Used as a binary cutoff to determine whether a gene is mutated in a sample or not. Default 0.
  # decreasing		Direction of ordering. Default: decreasing (higher frequency or mutation type first)
  #########################################################################################################
  
  ## Values: ###########################################################################
  #
  #	index			An integer vector the same length as the rows or columns of the matrix
  #					that was sorted. This vector can be used to sort the matrix. E.g.
  #					M[index,] or M[,index] for rows or columns, respectively.
  #
  ######################################################################################

  
  margin = match.arg(margin)
  
  if(length(ordering) == 1) {
    ordering = match.arg(ordering)
  }
  
  
  ## The rest of the function assumes ordering will be done by row. Transpose matrix if columns need sorting
  if(margin == "columns") {
    M = t(M)
  }
  
  ## Check to see if names or index was properly filled if selected
  if(ordering == "names" & !is.overlapping(names, rownames(M))) {
    stop(paste("names do not match matrix names for margin: ", margin, sep=""))
  }
  if(ordering == "index" & !is.overlapping(index, 1:nrow(M))) {
    stop(paste("index is out of range for matrix margin: ", margin, sep=""))
  }
  
  ## Handle NAs as lowested possible kind of mutation
  M.min = min(M, na.rm=TRUE)
  if(M.min > -1) { M.min = -1 }
  M[is.na(M) | is.infinite(M) | is.null(M)] = M.min
  
  ## Generate different matricies 
  index = 1:nrow(M)
  M.pos = M>mutation.cutoff  	# Used to sort by frequency
  M.any = M != 0           		# Used when sorting by type, will still sor NA's or other mutations < mutation.cutoff after the mutations > mutation.cutoff
  M.row.freq = rowSums(M.pos, na.rm=TRUE)
  M.column.freq = colSums(M.pos, na.rm=TRUE)
  
  ## Order samples by overall frequency, mutation type, or a given set of indices or labels
  if (ordering == "frequency") {
  
    M.combined = cbind(M.row.freq, M.pos, M.any, M)
    index = do.call("order", c(lapply(1:ncol(M.combined), function(x) M.combined[,x]), decreasing=decreasing))
    
  } else if (ordering == "frequency.mutation.type") {
    
    ## Sorts columns first by whether or not they are mutated and then by the type of mutation. Columns are ordered by frequency first
    temp.o = order(M.column.freq, decreasing=decreasing)    
    M.combined = cbind(M.pos[,temp.o], M.any[,temp.o], M[,temp.o])
    index = do.call("order", c(lapply(1:ncol(M.combined), function(x) M.combined[,x]), decreasing=decreasing))

  } else if (class(ordering) == "character" & is.overlapping(ordering, rownames(M))) {

    index = match(rownames(M), ordering) 

  } else if (ordering == "mutation.type") {
    
    ## Sorts columns first by whether or not they are mutated and then by the type of mutation.
    M.combined = cbind(M.pos, M.any, M)
    index = do.call("order", c(lapply(1:ncol(M.combined), function(x) M.combined[,x]), decreasing=decreasing))

  } else if (class(ordering) == "integer" & is.overlapping(ordering, 1:nrow(M))) {

    index = ordering

  } else if (ordering == "none") {

    index = 1:nrow(M)

  }

  return(index)
}






my.mapping = function() {
    require(RColorBrewer)
    variant.categories = data.frame(rbind(c(9, "Nonsense_Mutation", "Nonsense"),
								c(8, "Frame_Shift.*", "Frame Shift"),
								c(7, "[iI]n_[fF]rame.*", "In Frame Indel"),
								c(6, "Splice_Site", "Splice Site"),
								c(5, "Missense_Mutation", "Missense"),
								c(4, "Nonstop_.*|De_novo_Start.*|De_novo_Stop.*|Start_.*|Stop_.*|Translation_.*|Read\\-through.*", "Other"),
								c(3, "3'\\-?UTR|5'\\-?UTR|3'\\-?Flank|5'\\-?Flank|Intron", "UTR/Intron/Flanking"),
								c(2, "IGR", "IGR"),
								c(1, "RNA|lincRNA", "RNA")), stringsAsFactors=FALSE)
								
	return(variant.categories)	
}


tcga.mapping = function() {
  
  require(RColorBrewer)
  
  ## These 7 Categories are used in many of the TCGA papers
  id = c(-1, 1:6)
  maf.classification = c("Silent",
  						"Nonstop_.*|De_novo_Start.*|De_novo_Stop.*|Start_.*|Stop_.*|Translation_.*|Read\\-through.*",
  						"Missense_.*",
  						"Splice_Site.*",
  						"[iI]n_[fF]rame.*",
  						"Frame_Shift.*",
  						"Nonsense_.*")
  
  final.classification = c("Syn.", "Other non syn.", "Missense", "Splice site", "In frame indel", "Frame shift", "Nonsense")						
  
  brewer.pal.set1 = brewer.pal(7, "Set1")
  bg.col = c(brewer.pal(7, "Set1")[c(3, 1, 2, 4, 6, 5, 1)]) 
  
  return(data.frame(id, maf.classification, final.classification, bg.col, stringsAsFactors=FALSE))
}



do.call.args = function(FUN, required.args, additional.args=NULL) {
  
  ## Identifies overlap in the required args and user supplied args
  ## Overwrites required args using teh user supplied args.
  if(!is.null(additional.args) & !is.null(names(required.args)) & !is.null(names(additional.args))) {
    ind = !(names(required.args) %in% names(additional.args))
    do.call(FUN, c(required.args[ind], additional.args))
  } else if (!is.null(additional.args)) {
    do.call(FUN, c(required.args, additional.args))
  } else {
    do.call(FUN, required.args)
  }  
}


compute.cex <- function (label, width, height, units="inches", cex.lim=c(0.1, 1)) {

    #####################################################################################
	# Adam Gower, 2009
	# Modified by Josh Campbell 2016
	# Computes the cex factor needed to place a label (or the largest in a vector of labels) in a given plot region
	
	# label:  string of text to be drawn, or a character vector of text strings to be drawn, the largest of which will be used for computations
	# width:  width of the plot region in which to draw the text
	# height: height of the plot region in which to draw the text
	# units:  units in which width and height are reported
		
	# Get the width, height, and aspect ratio of the (largest) label
	label.width <- max(strwidth(label, units=units));
	label.height <- max(strheight(label, units=units));

    label.width.ratio =  label.width / width
    label.height.ratio = label.height / height
    
    cex.label = 1 / max(c(label.width.ratio, label.height.ratio))
    
    if(cex.label > cex.lim[2]) { cex.label = cex.lim[2] }
    if(cex.label < cex.lim[1]) { cex.label = cex.lim[1] }

	return(cex.label);
}


mar.to.mai = function() {
  return(max(par("mai") / par("mar"), na.rm=TRUE))
}

mar.to.fig.fraction = function(mar) {
  mar.fig.fraction = (mar * mar.to.mai()) / rep(rev(par("fin")), 2)  
  return(mar.fig.fraction)
}



is.overlapping = function(set1, set2) {
  ## Determines if two vectors contain the same values (not necessarily in the same order)
  flag = FALSE
  l1 = length(set1)
  l2 = length(set2)
  if(l1 == l2 & length(intersect(set1, set2)) == l1) {
    flag = TRUE
  }
  return(flag)
}


read.maf = function(filename, stringsAsFactors=FALSE) {
  maf = read.table(filename, header=TRUE, stringsAsFactors=stringsAsFactors, sep="\t", quote="")
  return(maf)
}


check.samples = function(samples, sample.list=NULL, verbose=FALSE) {
  ## Set up sample list
  s1 =  setdiff(sample.list, samples)
  s2 =  setdiff(samples, sample.list)
  if(!is.null(sample.list) & length(s1) > 0) {
    if(verbose == TRUE) { cat("The following samples were not in the maf or matrix:", "\n", paste(s1, collapse="\n"), "\n\nThey will be given zero counts across all genes or annotations\n") }
    samples.list = unique(c(sample.list, samples))
  } else if (!is.null(sample.list) & length(s2) > 0) {
    if(verbose == TRUE) { cat("The following samples in the maf or matrix were not in sample.list:", "\n", paste(s2, collapse="\n"), "\n\nThey will be excluded from the matrix and variant counts.\n") }
  } else {
    sample.list = unique(samples)
  }
  return(sample.list)
}

check.genes = function(genes, gene.list=NULL, verbose=TRUE)  {
  ## Set up gene list
  if(!is.null(gene.list)) {
    g1 = setdiff(gene.list, genes)
    if(length(g1) > 0) {
    if(verbose == TRUE) { cat("The following genes were not in the maf or matrix:", "\n", paste(g1, collapse="\n"), "\n\nThey will be given zero counts across all samples\n") }
    gene.list = unique(c(gene.list, genes))
    }
  } else {
    gene.list = unique(genes)
  }
  return(gene.list)
}


rescale.screen = function(screen.matrix) {
  if(ncol(screen.matrix) != 4) {
    stop("Matrix for split.screen must be four columns")
  }

  ## Rescale left and right boundries if they excede [0,1]
  if(min(screen.matrix[,1]) < 0) {
    screen.matrix[,1:2] = screen.matrix[,1:2] + abs(min(screen.matrix[,1]))
  }
  if(max(screen.matrix[,4]) > 1) {
    screen.matrix[,1:2] = screen.matrix[,1:2] / max(screen.matrix[,2])
  }
  
  ## Rescale top and bottom boundries if they excede [0,1]
  if(min(screen.matrix[,3]) < 0) {
    screen.matrix[,3:4] = screen.matrix[,3:4] + abs(min(screen.matrix[,3]))
  }
  if(max(screen.matrix[,4]) > 1) {
    screen.matrix[,3:4] = screen.matrix[,3:4] / max(screen.matrix[,4])
  }

  return(screen.matrix)
}


get.max.bar.from.table.list = function(table.list, params) {
  max.count = 0
  ntable = length(table.list)
  for(i in 1:ntable) {
    if(table.list[[i]]$type == "mutation") {
      params.required = list(counts=table.list[[i]]$counts, plot=FALSE) 
      b = do.call.args("gene.barplot", params.required, params)
      max.count = max(max.count, max(abs(b$lim)))  
    }  
  }
  return(max.count)
}




