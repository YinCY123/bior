myGOChord <- function (data, title, space, gene.order, gene.size, gene.space, 
                       nlfc = 1, lfc.col, lfc.min, lfc.max, ribbon.col, border.size, 
                       process.label, limit) 
{
  y <- id <- xpro <- ypro <- xgen <- ygen <- lx <- ly <- ID <- logFC <- NULL
  Ncol <- dim(data)[2]
  if (missing(title)) 
    title <- ""
  if (missing(space)) 
    space = 0
  if (missing(gene.order)) 
    gene.order <- "none"
  if (missing(gene.size)) 
    gene.size <- 3
  if (missing(gene.space)) 
    gene.space <- 0.2
  if (missing(lfc.col)) 
    lfc.col <- c("brown1", "azure", "cornflowerblue")
  if (missing(lfc.min)) 
    lfc.min <- -3
  if (missing(lfc.max)) 
    lfc.max <- 3
  if (missing(border.size)) 
    border.size <- 0.5
  if (missing(process.label)) 
    process.label <- 11
  if (missing(limit)) 
    limit <- c(0, 0)
  if (gene.order == "logFC") 
    data <- data[order(data[, Ncol], decreasing = T), ]
  if (gene.order == "alphabetical") 
    data <- data[order(rownames(data)), ]
  if (sum(!is.na(match(colnames(data), "logFC"))) > 0) {
    if (nlfc == 1) {
      cdata <- check_chord(data[, 1:(Ncol - 1)], limit)
      lfc <- sapply(rownames(cdata), function(x) data[match(x, 
                                                            rownames(data)), Ncol])
    }
    else {
      cdata <- check_chord(data[, 1:(Ncol - nlfc)], limit)
      lfc <- sapply(rownames(cdata), function(x) data[, 
                                                      (Ncol - nlfc + 1)])
    }
  }
  else {
    cdata <- check_chord(data, limit)
    lfc <- 0
  }
  #使基因数少于15个
  # while(dim(cdata)[1]>15){
  #   limit[1] = limit[1]+1
  #   if (sum(!is.na(match(colnames(data), "logFC"))) > 0) {
  #     if (nlfc == 1) {
  #       cdata <- check_chord(data[, 1:(Ncol - 1)], limit)
  #       lfc <- sapply(rownames(cdata), function(x) data[match(x, 
  #                                                             rownames(data)), Ncol])
  #     }
  #     else {
  #       cdata <- check_chord(data[, 1:(Ncol - nlfc)], limit)
  #       lfc <- sapply(rownames(cdata), function(x) data[, 
  #                                                       (Ncol - nlfc + 1)])
  #     }
  #   }
  #   else {
  #     cdata <- check_chord(data, limit)
  #     lfc <- 0
  #   }
  # }
  
  if (missing(ribbon.col)) 
    colRib <- grDevices::rainbow(dim(cdata)[2])
  else colRib <- ribbon.col
  nrib <- colSums(cdata)
  ngen <- rowSums(cdata)
  Ncol <- dim(cdata)[2]
  Nrow <- dim(cdata)[1]
  colRibb <- c()
  for (b in 1:length(nrib)) colRibb <- c(colRibb, rep(colRib[b], 
                                                      202 * nrib[b]))
  r1 <- 1
  r2 <- r1 + 0.1
  xmax <- c()
  x <- 0
  for (r in 1:length(nrib)) {
    perc <- nrib[r]/sum(nrib)
    xmax <- c(xmax, (pi * perc) - space)
    if (length(x) <= Ncol - 1) 
      x <- c(x, x[r] + pi * perc)
  }
  xp <- c()
  yp <- c()
  l <- 50
  for (s in 1:Ncol) {
    xh <- seq(x[s], x[s] + xmax[s], length = l)
    xp <- c(xp, r1 * sin(x[s]), r1 * sin(xh), r1 * sin(x[s] + 
                                                         xmax[s]), r2 * sin(x[s] + xmax[s]), r2 * sin(rev(xh)), 
            r2 * sin(x[s]))
    yp <- c(yp, r1 * cos(x[s]), r1 * cos(xh), r1 * cos(x[s] + 
                                                         xmax[s]), r2 * cos(x[s] + xmax[s]), r2 * cos(rev(xh)), 
            r2 * cos(x[s]))
  }
  df_process <- data.frame(x = xp, y = yp, id = rep(c(1:Ncol), 
                                                    each = 4 + 2 * l))
  xp <- c()
  yp <- c()
  logs <- NULL
  x2 <- seq(0 - space, -pi - (-pi/Nrow) - space, length = Nrow)
  xmax2 <- rep(-pi/Nrow + space, length = Nrow)
  for (s in 1:Nrow) {
    xh <- seq(x2[s], x2[s] + xmax2[s], length = l)
    if (nlfc <= 1) {
      xp <- c(xp, (r1 + 0.05) * sin(x2[s]), (r1 + 0.05) * 
                sin(xh), (r1 + 0.05) * sin(x2[s] + xmax2[s]), 
              r2 * sin(x2[s] + xmax2[s]), r2 * sin(rev(xh)), 
              r2 * sin(x2[s]))
      yp <- c(yp, (r1 + 0.05) * cos(x2[s]), (r1 + 0.05) * 
                cos(xh), (r1 + 0.05) * cos(x2[s] + xmax2[s]), 
              r2 * cos(x2[s] + xmax2[s]), r2 * cos(rev(xh)), 
              r2 * cos(x2[s]))
    }
    else {
      tmp <- seq(r1, r2, length = nlfc + 1)
      for (t in 1:nlfc) {
        logs <- c(logs, data[s, (dim(data)[2] + 1 - t)])
        xp <- c(xp, (tmp[t]) * sin(x2[s]), (tmp[t]) * 
                  sin(xh), (tmp[t]) * sin(x2[s] + xmax2[s]), 
                tmp[t + 1] * sin(x2[s] + xmax2[s]), tmp[t + 
                                                          1] * sin(rev(xh)), tmp[t + 1] * sin(x2[s]))
        yp <- c(yp, (tmp[t]) * cos(x2[s]), (tmp[t]) * 
                  cos(xh), (tmp[t]) * cos(x2[s] + xmax2[s]), 
                tmp[t + 1] * cos(x2[s] + xmax2[s]), tmp[t + 
                                                          1] * cos(rev(xh)), tmp[t + 1] * cos(x2[s]))
      }
    }
  }
  if (lfc[1] != 0) {
    if (nlfc == 1) {
      df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:Nrow), 
                                                      each = 4 + 2 * l), logFC = rep(lfc, each = 4 + 
                                                                                       2 * l))
    }
    else {
      df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:(nlfc * 
                                                             Nrow)), each = 4 + 2 * l), logFC = rep(logs, 
                                                                                                    each = 4 + 2 * l))
    }
  }
  else {
    df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:Nrow), 
                                                    each = 4 + 2 * l))
  }
  aseq <- seq(0, 180, length = length(x2))
  angle <- c()
  for (o in aseq) if ((o + 270) <= 360) 
    angle <- c(angle, o + 270)
  else angle <- c(angle, o - 90)
  #设置Term角度
  p1<-nrib[1]/sum(nrib)
  for (r in 2:length(nrib)) {
    perc <- nrib[r]/sum(nrib)
    p1 <- c(p1, p1[length(p1)]+perc)
  }
  p2<-p1[1]/2
  for(k in 2:length(p1)){
    p2[k]=p1[(k-1)]+nrib[k]/(sum(nrib)*2)
  }
  aseq1 <- p2*180
  angle1 <- c()
  for (o in aseq1) if ((o + 270) <= 360) 
    angle1 <- c(angle1, 90-o)
  else angle1 <- c(angle1, 90-o)
  
  df_texg <- data.frame(xgen = (r1 + gene.space) * sin(x2 + 
                                                         xmax2/2), ygen = (r1 + gene.space) * cos(x2 + xmax2/2), 
                        labels = rownames(cdata), angle = angle)
  df_texp <- data.frame(xpro = (r1 + 0.25) * sin(x + xmax/2), 
                        ypro = (r1 + 0.25) * cos(x + xmax/2), labels = colnames(cdata), 
                        stringsAsFactors = FALSE)
  cols <- rep(colRib, each = 4 + 2 * l)
  x.end <- c()
  y.end <- c()
  processID <- c()
  for (gs in 1:length(x2)) {
    val <- seq(x2[gs], x2[gs] + xmax2[gs], length = ngen[gs] + 
                 1)
    pros <- which((cdata[gs, ] != 0) == T)
    for (v in 1:(length(val) - 1)) {
      x.end <- c(x.end, sin(val[v]), sin(val[v + 1]))
      y.end <- c(y.end, cos(val[v]), cos(val[v + 1]))
      processID <- c(processID, rep(pros[v], 2))
    }
  }
  df_bezier <- data.frame(x.end = x.end, y.end = y.end, processID = processID)
  df_bezier <- df_bezier[order(df_bezier$processID, -df_bezier$y.end), 
  ]
  x.start <- c()
  y.start <- c()
  for (rs in 1:length(x)) {
    val <- seq(x[rs], x[rs] + xmax[rs], length = nrib[rs] + 
                 1)
    for (v in 1:(length(val) - 1)) {
      x.start <- c(x.start, sin(val[v]), sin(val[v + 1]))
      y.start <- c(y.start, cos(val[v]), cos(val[v + 1]))
    }
  }
  df_bezier$x.start <- x.start
  df_bezier$y.start <- y.start
  df_path <- bezier(df_bezier, colRib)
  if (length(df_genes$logFC) != 0) {
    tmp <- sapply(df_genes$logFC, function(x) ifelse(x > 
                                                       lfc.max, lfc.max, x))
    logFC <- sapply(tmp, function(x) ifelse(x < lfc.min, 
                                            lfc.min, x))
    df_genes$logFC <- logFC
  }
  g <- ggplot() + geom_polygon(data = df_process, aes(x, y,group = id), 
                               fill = 'white', inherit.aes = F, color = 'white') + #通路方块底
    geom_polygon(data = df_process, aes(x, y, group = id), 
                 fill = cols, inherit.aes = F, alpha = 1, color = cols) + #通路方块
    # geom_point(aes(x = xpro, y = ypro, size = factor(labels,levels = labels), 
    #                shape = NA), data = df_texp) + 
    #guides(size = guide_legend("Term", ncol = 2, byrow = T,override.aes = list(shape = 22, fill = unique(cols),size = 8))) + 
    #theme(legend.text = element_text(size = process.label)) + 
    geom_text(aes(xgen, ygen, label = labels, angle = angle),
              data = df_texg, size = gene.size) +#基因名
    geom_polygon(aes(x = lx,
                     y = ly, group = ID), data = df_path, fill = colRibb, size = border.size, inherit.aes = F, alpha = 0.6) +#弦
    geom_text(aes(xpro, ypro, label = labels),data = df_texp,angle = angle1,size=gene.size) +#通路名
    labs(title = title) + theme_blank
  if (nlfc >= 1) {
    g=g + geom_polygon(data = df_genes, aes(x, y, group = id,
                                          fill = logFC), inherit.aes = F, color = 'white') +
      scale_fill_gradient2("logFC", space = "Lab", low = lfc.col[3],
                           mid = lfc.col[2], high = lfc.col[1],
                           guide = guide_colorbar(title.position = "top",
                                                  title.hjust = 0.5),
                           breaks = c(min(df_genes$logFC),
                                      max(df_genes$logFC)), labels = c(round(min(df_genes$logFC)),
                                                                       round(max(df_genes$logFC)))) +
      theme(legend.position = "bottom",
            legend.background = element_rect(fill = "transparent"),
            legend.box = "horizontal", legend.direction = "horizontal",
            legend.key.size = unit(12,'pt'),
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 8),#legend.spacing = unit(1,"cm"),
            plot.margin = unit(c(1,1,1,1),"lines"))
  }
  else {
    g=g + geom_polygon(data = df_genes, aes(x, y, group = id),
                     fill = logFC, inherit.aes = F, color = "black") +
      theme(legend.position = "bottom", legend.background = element_rect(fill = "transparent"),
            legend.box = "horizontal", legend.direction = "horizontal")#vertical
  }
  return(g)
}
environment(myGOChord) <- environment(GOChord)


mychord_dat=function (data, genes, process){
  id <- term <- logFC <- BPprocess <- NULL
  colnames(data) <- tolower(colnames(data))
  if (missing(genes)) {
    if (is.null(data$logFC)) {
      genes <- as.character(unique(data$genes))
    }
    else {
      genes <- subset(data, !duplicated(genes), c(genes, 
                                                  logFC))
    }
  }
  else {
    if (is.vector(genes)) {
      genes <- as.character(genes)
    }
    else {
      if (class(genes[, 2]) != "numeric") 
        genes[, 2] <- as.numeric(levels(genes[, 2]))[genes[, 
                                                           2]]
      genes[, 1] <- as.character(genes[, 1])
      colnames(genes) <- c("genes", "logFC")
    }
  }
  if (missing(process)) {
    process <- as.character(unique(data$term))
  }
  else {
    if (class(process) != "character") 
      process <- as.character(process)
  }
  if (strsplit(process[1], ":")[[1]][1] == "GO"|substr(process[1],1,3)%in% c("hsa","mmu")) {
    subData <- subset(data, id %in% process)
    colnames(subData)[which(colnames(subData) == "id")] <- "BPprocess"
  }
  else {
    subData <- subset(data, term %in% process)
    colnames(subData)[which(colnames(subData) == "term")] <- "BPprocess"
  }
  if (is.vector(genes)) {
    M <- genes[genes %in% unique(subData$genes)]
    mat <- matrix(0, ncol = length(process), nrow = length(M))
    rownames(mat) <- M
    colnames(mat) <- process
    for (p in 1:length(process)) {
      sub2 <- subset(subData, BPprocess == process[p])
      for (g in 1:length(M)) mat[g, p] <- ifelse(M[g] %in% 
                                                   sub2$genes, 1, 0)
    }
  }
  else {
    genes <- subset(genes, genes %in% unique(subData$genes))
    N <- length(process) + 1
    M <- genes[, 1]
    mat <- matrix(0, ncol = N, nrow = length(M))
    rownames(mat) <- M
    colnames(mat) <- c(process, "logFC")
    mat[, N] <- genes[, 2]
    for (p in 1:(N - 1)) {
      sub2 <- subset(subData, BPprocess == process[p])
      for (g in 1:length(M)) mat[g, p] <- ifelse(M[g] %in% 
                                                   sub2$genes, 1, 0)
    }
  }
  return(mat)
}

environment(mychord_dat) <- environment(chord_dat)