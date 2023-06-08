##########
# Import #
##########

# library("gdata", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
library('readxl')
library('gsubfn')
library('tidyr')
library("ChemmineR")
library('ggplot2')
library('nplr')
library('xlsx')
library('dplyr')
library('ggthemes')
library('magick')

#############
# Functions #
#############

# Function to parse the name of the compound
parse_cmpd <- function(x) {
  out <- paste(gsubfn(".", list(" " = "", "-" = "", "+" = "", "(" = "/("),
                      x), '/', sep = '')
  return(out)
}

# Function to parse the name of the compound from ambinter
substramb <- function(x, n = 4) {
  paste(substr(x, 0, 3),
        substr(x, nchar(x) - n + 1, nchar(x)), sep = '')
}

# Function to create a dataframe with the names of all compounds registered in the project
gx.cmpd.name <- function(synth.name,
                         df.list.cmpd = data.frame(
                          read_excel("/home/dreano/Desktop/mppase/Shared/20200220_mPPaseCompoundRegister_updated.xlsx", 2))) {
  # src_path = "/home/dreano/Desktop/mppase/Shared/mPPaseCompoundRegister_updated 23022018.xlsx"
  # df.list.cmpd = data.frame(read_excel(src_path,2))
  parse.synth.name <- gsubfn(".", list(" " = "", "-" = "", "+" = ""), synth.name)
  l.gx.names <- apply(df.list.cmpd['Synthetic'], 2, parse_cmpd)
  df.gx.names <- cbind(df.list.cmpd, l.gx.names)
  colnames(df.gx.names) <- c(colnames(df.list.cmpd), 'l.synth.name')
  if (substr(synth.name, 1, 4) == 'mPP-') {
    synth.id <- which(df.list.cmpd$Generic == synth.name)
    gx.name <- df.list.cmpd[synth.id, c('SMILES', 'Generic', 'Paper')]
  }else if (substr(parse.synth.name, 0, 3) == 'Amb') {
    df.amb <- data.frame(l.synth.name =
                          substramb(
                            as.character(
                              df.gx.names$l.synth.name)))
    gx.name <- df.list.cmpd[grep(paste(parse.synth.name,
                                       '/',
                                       sep = ''),
                                 df.amb$l.synth.name),
                            c(2, 3, 4)]
  }else {
    gx.name <- df.list.cmpd[grep(paste(parse.synth.name,
                                       '/',
                                       sep = ''),
                                 df.gx.names$l.synth.name),
                            c(2, 3, 4)]
  }
  if (nrow(gx.name) == 0) {
    gx.name <- data.frame(SMILES = NA, Generic = synth.name, Paper = NA)
  }
  return(gx.name)
}

# Function to read the IC50s from the excel file
parsing_IC50 <- function(src_path = "~/Desktop/mppase/preliminary_IC50s-2018.xlsx") {
  src_path_1 <- "/home/dreano/Desktop/mppase/Shared/20200220_mPPaseCompoundRegister_updated.xlsx"
  df.list.cmpd <- data.frame(read_excel(src_path_1, 2))
  l_sheets <- excel_sheets(src_path)
  estimate.IC50 <- lapply(l_sheets, read_excel,
                          path = src_path,
                          col_names = FALSE)
  df.out <- data.frame()
  nb_noname <- 0
  nb_id <- 0
  for (i in 1:(length(l_sheets))) {
    summary.IC50 <- data.frame(estimate.IC50[[i]][, seq(1, 5)])
    not.summary.IC50 <- data.frame(estimate.IC50[[i]][, -seq(1, 5)])
    id.mean <- as.data.frame(which(not.summary.IC50 == 'mean', T))
    if (nrow(id.mean) == 0) {
      next
    }
    id.mean <- id.mean[order(id.mean$row),]
    for (j in 1:nrow(id.mean)) {
      name.cmpd <- not.summary.IC50[id.mean$row[j] + 2, id.mean$col[j] - 8]
      if (is.na(name.cmpd) |
        name.cmpd == 1000 |
        name.cmpd == 100) {
        name.cmpd <- not.summary.IC50[id.mean$row[j] + 2, id.mean$col[j] - 9]
      }
      if (is.na(name.cmpd) |
        name.cmpd == 1000 |
        name.cmpd == 100) {
        name.cmpd <- not.summary.IC50[id.mean$row[j] + 2, id.mean$col[j] - 10]
      }
      if (length(name.cmpd) == 0) {
        nb_noname <- nb_noname + 1
        next
      }
      if (is.na(name.cmpd) | name.cmpd == 0) {
        nb_noname <- nb_noname + 1
        next
      }
      id.range <- seq(max(id.mean$col[j] - 10, 0)
        , id.mean$col[j])
      id.col <- which(not.summary.IC50[id.mean$row[j] + 1,] == 0 |
                       not.summary.IC50[id.mean$row[j] + 1,] == 0.1)
      id.rm.1 <- which(id.col > id.mean$col[j] |
                        id.col < id.mean$col[j] - 10)
      if (length(id.rm.1) > 0) {
        id.col <- id.col[-id.rm.1]
      }
      id.na <- which(is.na(not.summary.IC50[, id.mean$col[j]]))
      id.row <- id.na[which(id.na > id.mean$row[j])]
      id.seq.1 <- seq(id.mean$row[j] + 1, id.na[which(id.na > id.mean$row[j])][1] - 1)
      x.axis <- not.summary.IC50[id.seq.1, id.col]
      x.axis <- as.numeric(replace(x.axis, x.axis == 0, 0.1))
      id.100 <- head(which(not.summary.IC50[id.mean$row[j] + 1,] == 100), -1)
      id.rm <- which(is.na(not.summary.IC50[id.mean$row[j] + 2, id.100]) |
                      id.100 >= id.mean$col[j] |
                      id.100 < id.mean$col[j] - 10)
      if (length(id.rm) > 0) {
        id.100 <- id.100[-id.rm]
      }
      x.axis <- rep(x.axis, length(id.100))
      y.axis <- as.numeric(gather(not.summary.IC50[id.seq.1, id.100])$value)
      nb_id <- nb_id + 1
      df.translate <- gx.cmpd.name(as.character(name.cmpd), df.list.cmpd)
      df.tmp <- data.frame(df.translate$Generic,
                           df.translate$Paper,
                           name.cmpd,
                           x.axis,
                           y.axis,
                           nb_id,
                           df.translate$SMILES)
      colnames(df.tmp) <- c('gx.name', 'paper.name', 'synth.name', 'x', 'y', 'id', 'SMILES')
      df.out <- rbind(df.out, df.tmp)
    }
  }
  print(paste(nb_noname, 'compound whith no name'))
  return(df.out)
}

# Function to plot the IC50s based on nplr package
plot.IC50 <- function(df.result = parsing_IC50(), conf.level = 0.9, showSDerr = F,
                      showInfl = F, showPoints = T, showCI = T, plot_curve = F,
                      showGOF = F, showEstim = F, save_pic = F, up.lim = 5e3, n = 5,
                      target = 'mPPase') {
  set.seed(123) # for reproducibility
  setwd('~/Desktop/mppase/') # set working directory
  df.return <- data.frame()
  l.result <- split(df.result, df.result$id) # split by id
  # loop over the different compounds
  for (i in 1:length(l.result)) {
    df.cmpd <- l.result[[i]]
    df.cmpd$x <- df.cmpd$x * 10^-9
    np.cmpd <- nplr(df.cmpd$x, df.cmpd$y / 100) # Y axis scale from 0 to 1
    name.cmpd <- df.cmpd$gx.name[1]
    paper.name <- df.cmpd$paper.name[1]
    synt.name <- df.cmpd$synth.name[1]
    smiles.cmpd <- df.cmpd$SMILES[1]
    id.cmpd <- df.cmpd$id
    # plot the IC50 curve
    if (!grepl('SI-', paper.name)) { # Color the curve for compound in the paper
      linecol <- 'skyblue1'
      legcol <- 'skyblue3'
    }else {                                # Color the curve for compound in the supplementary
      linecol <- 'black'
      legcol <- 'black'
    }
    if (is.na(paper.name)) {
      p.title <- paste(target, "response to", name.cmpd)
      x.title <- paste('[', name.cmpd, '] (M)', sep = '')
      r.path <- "./Picture_260521/"
      curve_name <- paste(r.path, "curve/", 'Curve_IC50_', name.cmpd, '_', synt.name, '_', id.cmpd, '.png', sep = '')
    }else {
      # p.title = paste("mPPase response to §",paper.name," ( ",name.cmpd, " )",
      p.title <- paste(target, "response to", paper.name)
      # x.title = paste('[ §',paper.name,' ] (M)',sep='')
      x.title <- paste('[', paper.name, '] (M)', sep = '')
      r.path <- "./Picture_260521/Paper/"
      curve_name <- paste(r.path, "curve/", paper.name, '_Curve_IC50_', name.cmpd, '_', synt.name, '_', id.cmpd, '.png', sep = '')
    }
    IC.50 <- getEstimates(np.cmpd, targets = 0.5, conf.level = conf.level) * 10^6
    if (IC.50$y != 5e+05 |
      IC.50$x == Inf |
      IC.50$x < 1e-5 |
      IC.50$x > up.lim) {
      IC.50 <- data.frame(y = 0, x.05 = NA, x = NA, x.90 = NA) # if no IC50 found
      lgd.ic <- bquote(
        bold(
          atop(NA,
               atop(IC[50] * ' = NA',
                    '[NA–NA]'))))
    }else {
      if (IC.50[[3]] > 10) {
        r.ic50 <- round(IC.50[[3]])
      }else {
        r.ic50 <- format(round(IC.50[[3]], digits = 1), nsmall = 1) }
      if (IC.50[[2]] > 10) {
        rinf.ic50 <- round(IC.50[[2]])
      }else {
        rinf.ic50 <- format(round(IC.50[[2]], digits = 1), nsmall = 1) }
      if (IC.50[[4]] > 10) {
        rsup.ic50 <- round(IC.50[[4]])
      }else {
        rsup.ic50 <- format(round(IC.50[[4]], digits = 1), nsmall = 1) }
      lgd.ic <- bquote(
        bold(
          atop(NA,
               atop(IC[50] * .(paste(' = ', r.ic50, ' µM', sep = '')),
                    .(paste('[',
                            rinf.ic50,
                            '–',
                            rsup.ic50,
                            ']', sep = ''))))))
    }
    df.IC50 <- data.frame(name.cmpd, paper.name, synt.name, IC.50$x, IC.50[2], IC.50[4], smiles.cmpd)
    colnames(df.IC50) <- c("name.cmpd", 'paper.name', 'synth.name', "IC50", "IC50.low", "IC50.up", 'SMILES')
    if (plot_curve) {
      plot(np.cmpd, pcol = "grey40", lcol = linecol,
           showEstim = showEstim, showInfl = showInfl,
           showPoints = showPoints, showCI = showCI,
           showGOF = showGOF, cex.main = 1.5,
           showSDerr = showSDerr, lwd = 3,
           main = p.title,
           conf.level = conf.level,
           unit <- ' µM',
           xlab = '',
           ylab = paste(target, 'activity'), ylim = c(0, 1.2), xlim = c(-10, -4), cex.lab = 1.2, las = 1)
      text(labels = lgd.ic, x = -9, y = 0.15, font = 2, col = legcol, cex = 1.5)
      title(xlab = bquote(Log[10] * .(x.title)))
      if (save_pic) {
        png(curve_name,
            width = 595.3684 * n, height = 400 * n, res = 100 * n, family = "serif")
        plot(np.cmpd, pcol = "grey40", lcol = linecol,
             showEstim = showEstim, showInfl = showInfl,
             showPoints = showPoints, showCI = showCI,
             showGOF = showGOF, cex.main = 1.5,
             showSDerr = showSDerr, lwd = 3,
             main = p.title,
             conf.level = conf.level,
             unit = ' µM',
             xlab = '',
             ylab = 'mPPase activity', ylim = c(0, 1.2), xlim = c(-10, -4), cex.lab = 1.2, las = 1)
        text(labels = lgd.ic, x = -9, y = 0.15, font = 2, col = legcol, cex = 1.5)
        title(xlab = bquote(Log[10] * .(x.title)))
        dev.off()
      }
      split.df.cmpd <- split(df.cmpd[4:5], df.cmpd$x)
    }
    df.return <- rbind(df.return, df.IC50)
  }
  return(df.return[order(df.return$IC50),])
}

# Read new batch of data
read.cmpd.test <- function(src_path = "~/Desktop/mppase/New compounds_170119.xlsx",
                           nb_replic = 4, nb_id = 0, start_sheet = 1, end_sheet = 0,
                           skip_col = 1, skip_row = 8) {
  src_path_1 <- "/home/dreano/Desktop/mppase/Shared/20200220_mPPaseCompoundRegister_updated.xlsx"
  df.list.cmpd <- data.frame(read_excel(src_path_1, 2))
  l_sheets <- excel_sheets(src_path)
  cmpd.data <- lapply(l_sheets, read_excel,
                      path = src_path,
                      col_names = FALSE)
  df.out <- data.frame()
  df.ini <- data.frame(gx.name = as.character(),
                       paper.name = as.character(),
                       synth.name = as.character(),
                       x = as.numeric(),
                       y = as.numeric(),
                       SMILES = as.character())
  nb_noname <- 0
  if (end_sheet == 0) end_sheet <- length(l_sheets)
  if (length(skip_row) == 1) skip_row <- rep(skip_row, end_sheet - start_sheet + 1)
  for (i in start_sheet:end_sheet) {
    df.full <- data.frame(cmpd.data[[i]])
    act.id <- data.frame(which(df.full == '% Activity', T))
    if (nrow(act.id) == 0) next
    row.id <- seq(act.id$row + 2, act.id$row + 9)
    plate.scheme <- df.full[row.id, seq(2, 12, nb_replic)]
    col.id <- seq(act.id$col + 1, act.id$col + 12)
    inib.data <- df.full[row.id, col.id]
    for (r in seq(1, nrow(inib.data))) {
      for (c in seq(1, ncol(inib.data))) {
        if (r <= skip_row[i - start_sheet + 1] & c <= skip_col * nb_replic) next()
        exp.name <- plate.scheme[r, ceiling(c / nb_replic)]
        if (is.na(exp.name)) next()
        exp.info <- strsplit(exp.name, ' ')[[1]]
        if (is.na(exp.info[1])) next()
        name.cmpd <- paste(exp.info[1:length(exp.info) - 1], collapse = ' ')
        y <- inib.data[r, c]
        x <- as.numeric(gsubfn(".", list("u" = "", "M" = "", "n" = ""),
                               exp.info[length(exp.info)])) * 1000 # convert to nM #
        df.translate <- gx.cmpd.name(as.character(name.cmpd), df.list.cmpd)
        df.tmp <- data.frame(df.translate$Generic,
                             df.translate$Paper,
                             name.cmpd,
                             x,
                             as.numeric(y),
                             df.translate$SMILES)
        colnames(df.tmp) <- c('gx.name', 'paper.name', 'synth.name', 'x', 'y', 'SMILES')
        df.ini <- rbind(df.ini, df.tmp)
      }
    }
  }
  df.split <- split(df.ini, df.ini$synth.name)
  for (l in seq(1, length(df.split))) {
    df.split[[l]] <- rbind(data.frame(df.split[[l]][1, c(1, 2, 3, 6)], x = 0.1, y = 100),
                           df.split[[l]])
    df.split[[l]] <- cbind(df.split[[l]], id = nb_id + l)
    df.out <- rbind(df.out, df.split[[l]][, c(1, 2, 3, 5, 6, 7, 4)])

  }
  return(df.out)
}

# cluster the IC50 dataframe summary
clust.df.out <- function(df.no.clust) {
  rm.name <- c()
  df.no.clust <- df.no.clust[order(df.no.clust$IC50),]
  df.out <- data.frame()
  for (name in df.no.clust$name.cmpd) {
    if (name %in% rm.name) {
      next()
    }
    rm.name <- cbind(rm.name, name)
    df.out <- rbind(df.out, df.no.clust[which(df.no.clust$name.cmpd == name),])
  }
  return(df.out)
}

# update the summary of IC50 file
update.summary.IC50 <- function(file_src, file_out, df.new) {
  df.old <- read_excel(file_src, 1)
  df.update <- rbind(df.old, df.new)
  df.out <- as.data.frame(clust.df.out(df.update))
  write.xlsx(df.out, file_out, col.names = F)
  return(df.out)
}

# update/add new compound to the register of compounds file
update.cmpd.register <- function(ic50_sum, cpd_reg, file_out) {
  df.sum <- read_excel(ic50_sum, 1)
  df.register <- read_excel(cpd_reg, 2)
  instr <- read_excel(cpd_reg, 1)
  # df.1 = data.frame(Generic = df.sum$name.cmpd, new_Ic50 = round(df.sum$IC50*(10^-3),3))
  df.1 <- data.frame(Generic = df.sum$name.cmpd, new_Ic50 = paste(round(df.sum$IC50 * (10^-3), 3), ' [', round(df.sum$IC50.low * (10^-3), 3), ';', new_Ic50 = round(df.sum$IC50.up * (10^-3), 3), ']'))
  df.conc1 <- mutate(group_by(df.1, Generic), new_Ic50 = paste0(new_Ic50, collapse = ",\n"))
  df.1.clean <- df.conc1[!duplicated.data.frame(df.conc1),]
  l_col <- c(colnames(df.register), 'new_Ic50')
  df.2 <- merge(df.register, df.1.clean, all.x = T)
  df.3 <- df.2[, l_col]
  write.xlsx(instr, file_out, sheetName = 'Instructions')
  write.xlsx(df.3, file_out, row.names = F, append = T, sheetName = 'CompoundRegister')
  # write.xlsx(instr,file_out, sheetName = 'Instructions')
  # write.xlsx(df.3,file_out,row.names = F,append = T, sheetName = 'CompoundRegister')
  return(df.3)
}

setwd('~/Desktop/mppase/')

# crop the curve image for the paper
crop_folder <- function(folder_name = 'Picture_050719/Paper/curve/') {
  for (file_name in list.files()) {
    curve <- image_read(file_name)
    curve_crop <- image_crop(curve, '2976x1900')
    image_write(curve_crop, paste('../crop/', file_name, sep = ''))
  }
}

# add new data to the file compatible with Canvas
update_IC50_csv <- function(ic50_tab, write_csv = F, csv_out = '') {
  fold_out <- 'Canvas_file/update_ic50/'
  `Generic code` <- ic50_tab$name.cmpd
  IC50 <- round(ic50_tab$IC50, 3)
  `new IC50 [CI95%] (uM)` <- paste(round(ic50_tab$IC50, 3), ' [ ', round(ic50_tab$IC50.low, 3), ' - ', round(ic50_tab$IC50.up, 3), ' ]')
  tab.out <- data_frame(`Generic code`, IC50, `new IC50 [CI95%] (uM)`)
  if (write_csv) {
    write.csv(tab.out, paste(fold_out, csv_out, sep = ''), quote = F, row.names = F, na = '', sep = '')
  }
  return(tab.out)
}

########
# Main #
########

# Initialisation of the different reports files
# write.table(res.IC50,"Summary_IC50_240519.csv",col.names = T,row.names = F)
#res.IC50.clust = clust.df.out(res.IC50)
#write.xlsx(res.IC50.clust,"Summary_IC50_240519.xlsx",row.names = F)
# write.table(res.IC50.clust,"Summary_IC50_clust.csv",col.names = T,row.names = F)

# Parameters for the run (show/saving pictures, saving dataframes)
pics <- T
save_p <- T
df.out <- read.cmpd.test(nb_id = 215)

# Read of experimental data
df.ki <- parsing_IC50()
df.out <- read.cmpd.test('IC50_determination_of_compounds.xlsx', nb_replic = 3, nb_id = 240)
df.cmpd2205 <- read.cmpd.test('KV_inhibition assay (Sep-Nov2017).xlsx', nb_replic = 4, nb_id = 0)
df.cmpd1701 <- read.cmpd.test('New compounds_170119.xlsx', nb_replic = 4, nb_id = 0)
df.cmpd0601.1 <- read.cmpd.test('Aaron_Inhibition test and IC50_310118.xlsx', nb_replic = 4, nb_id = 0, end_sheet = 1)
df.cmpd0601.2 <- read.cmpd.test('Aaron_Inhibition test and IC50_310118.xlsx', nb_replic = 3, nb_id = 0, start_sheet = 2)
df.cmpd250418 <- read.cmpd.test('Inhibition_250418_result.xlsx', nb_replic = 4, nb_id = 0)
df.cmpd1701_ <- read.cmpd.test('New compounds_170119_.xlsx', nb_replic = 4, nb_id = 0, start_sheet = 5)
df.cmpd130619 <- read.cmpd.test('AKI samples mPPase Spring 2019.xlsx', nb_replic = 4, nb_id = 0)
df.cmpd180120 <- read.cmpd.test('Testing_130220.xlsx', nb_replic = 4, nb_id = 0, skip_col = 1, skip_row = c(3, 3, 3, 6), start_sheet = 2)
df.cmpd181119 <- read.cmpd.test('Loic_testing_181119.xlsx', nb_replic = 4, nb_id = 0)
df.cmpd010716 <- read.cmpd.test('Retest_July_2016_J&K.xlsx', nb_replic = 4, nb_id = 0)
df.cmpd011020.1 <- read.cmpd.test('IC50 Alex cpds 2020.xlsx', nb_replic = 3, nb_id = 0)
df.cmpd011020.2 <- read.cmpd.test('New compounds_Oct2020.xlsx', nb_replic = 4, nb_id = 0)
df.BLZ.mppase <- read.cmpd.test('IC50_TmPPase_format.xlsx', nb_replic = 3, nb_id = 0)
df.BLZ.Pfppase <- read.cmpd.test('IC50_PfPPase_format.xlsx', nb_replic = 3, nb_id = 0)
df.UR8.mppase <- read.cmpd.test('UR-8_IC50_TmPPase_format.xlsx', nb_replic = 3, nb_id = 0)
df.UR8.Pfppase <- read.cmpd.test('UR-8_IC50_PfPPase_format.xlsx', nb_replic = 3, nb_id = 0)
df.cmpd021221 <- read.cmpd.test('Cpds_021221.xlsx', nb_replic = 4, nb_id = 0)

# Plot of the IC50 curves
IC50.new.cmpd <- plot.IC50(df.out, plot_curve = pics, save_pic = save_p)
IC50.new.cmpd.2 <- plot.IC50(df.cmpd2205, plot_curve = pics, save_pic = save_p)
IC50.new.cmpd.3 <- plot.IC50(df.cmpd1701, plot_curve = pics, save_pic = save_p)
IC50.new.cmpd.4 <- plot.IC50(df.cmpd0601.1, plot_curve = pics, save_pic = save_p)
IC50.new.cmpd.5 <- plot.IC50(df.cmpd0601.2, plot_curve = pics, save_pic = save_p)
IC50.new.cmpd.6 <- plot.IC50(df.cmpd250418, plot_curve = pics, save_pic = save_p)
IC50.new.cmpd.7 <- plot.IC50(df.cmpd1701_, plot_curve = pics, save_pic = save_p)
IC50.new.cmpd.8 <- plot.IC50(df.cmpd130619, plot_curve = pics, save_pic = save_p)
res.IC50 <- plot.IC50(df.ki, plot_curve = pics, save_pic = save_p)
IC50.new.cmpd.180120 <- plot.IC50(df.cmpd180120, plot_curve = pics, save_pic = save_p)
IC50.new.cmpd.181119 <- plot.IC50(df.cmpd181119, plot_curve = pics, save_pic = save_p)
IC50.new.cmpd.010716 <- plot.IC50(df.cmpd010716, plot_curve = pics, save_pic = save_p)
cmpd011020.1 <- plot.IC50(df.cmpd011020.1, plot_curve = pics, save_pic = save_p)
cmpd011020.2 <- plot.IC50(df.cmpd011020.2, plot_curve = pics, save_pic = save_p)
BLZ.mppase <- plot.IC50(df.BLZ.mppase, plot_curve = pics, save_pic = save_p)
BLZ.Pfppase <- plot.IC50(df.BLZ.Pfppase, plot_curve = pics, save_pic = save_p, target = 'fPPase')
UR8.mppase <- plot.IC50(df.UR8.mppase, plot_curve = pics, save_pic = save_p)
UR8.Pfppase <- plot.IC50(df.UR8.Pfppase, plot_curve = pics, save_pic = save_p, target = 'fPPase')
cmpd021221 <- plot.IC50(df.cmpd021221, plot_curve = pics, save_pic = save_p)

# writing the IC50 in a csv file
IC50.new.cmpd.csv <- update_IC50_csv(IC50.new.cmpd, F)
IC50.new.cmpd.2.csv <- update_IC50_csv(IC50.new.cmpd.2, F)
IC50.new.cmpd.3.csv <- update_IC50_csv(IC50.new.cmpd.3, F)
IC50.new.cmpd.4.csv <- update_IC50_csv(IC50.new.cmpd.4, F)
IC50.new.cmpd.5.csv <- update_IC50_csv(IC50.new.cmpd.5, F)
IC50.new.cmpd.6.csv <- update_IC50_csv(IC50.new.cmpd.6, F)
IC50.new.cmpd.7.csv <- update_IC50_csv(IC50.new.cmpd.7, F)
IC50.new.cmpd.8.csv <- update_IC50_csv(IC50.new.cmpd.8, T, 'IC50.new.cmpd.8.csv')
res.IC50.csv <- update_IC50_csv(res.IC50, F)
IC50.new.cmpd.180120.csv <- update_IC50_csv(IC50.new.cmpd.180120, T, 'IC50.new.cmpd.180120.csv')
IC50.new.cmpd.181119.csv <- update_IC50_csv(IC50.new.cmpd.181119, F)
IC50.new.cmpd.010716.csv <- update_IC50_csv(IC50.new.cmpd.010716, T, 'IC50.new.cmpd.010716.csv')
cmpd011020.1.csv <- update_IC50_csv(cmpd011020.1, T, 'cmpd011020.1.csv')
cmpd011020.2.csv <- update_IC50_csv(cmpd011020.2, T, 'cmpd011020.2.csv')

# Plots for analysis/report/meeting to comment if not needed
all <- rbind(IC50.new.cmpd, IC50.new.cmpd.2, IC50.new.cmpd.3, IC50.new.cmpd.4,
             IC50.new.cmpd.5, IC50.new.cmpd.6, IC50.new.cmpd.7, IC50.new.cmpd.8,
             res.IC50, IC50.new.cmpd.180120, IC50.new.cmpd.181119, IC50.new.cmpd.010716)

which(all$synth.name == 'AKI-A157')


tmp_ <- read_excel('/home/dreano/Desktop/mppase/Shared/mPPaseCompoundRegister_updated_25102019.xlsx', sheet = 2)
write.table(all, file = 'ALL_dataframe.txt')
all <- read.table('ALL_dataframe.txt')
all[which(substr(all$name.cmpd, 1, 3) != 'mPP'),]
azu_code <- c('mPP-0350', 'mPP-0351', 'mPP-0352', 'mPP-0078', 'mPP-0086', 'mPP-0118',
              'mPP-0267', 'mPP-0079', 'mPP-0117', 'mPP-0091', 'mPP-0080', 'mPP-0268',
              'mPP-0085', 'mPP-0096', 'mPP-0130', 'mPP-076', 'mPP-0178', 'mPP-0177',
              'mPP-0120', 'mPP-0114', 'mPP-0087', 'mPP-0181', 'mPP-0182', 'mPP-0180',
              'mPP-0134', 'mPP-0129', 'mPP-0179', 'mPP-0133', 'mPP-0225', 'mPP-0093',
              'mPP-0090', 'mPP-0135', 'mPP-0089', 'mPP-0095', 'mPP-0226', 'mPP-0269',
              'mPP-0270', 'mPP-0084', 'mPP-0092', 'mPP-0132', 'mPP-0185', 'mPP-0184',
              'mPP-0187', 'mPP-0186', 'mPP-0131', 'mPP-0183', 'mPP-0088', 'mPP-0103',
              'mPP-0119', 'mPP-0112', 'mPP-0137', 'mPP-0138', 'mPP-0115', 'mPP-0102',
              'mPP-0100', 'mPP-0104', 'mPP-0097', 'mPP-0094', 'mPP-0082', 'mPP-081',
              'mPP-083', 'mPP-0249', 'mPP-0250', 'mPP-0251', 'mPP-0279', 'mPP-0280',
              'mPP-0281', 'mPP-0282', 'mPP-0283', 'mPP-0284', 'mPP-0285', 'mPP-0253',
              'mPP-0254', 'mPP-0255', 'mPP-0078', 'mPP-0086', 'mPP-0118', 'mPP-0267',
              'mPP-0079', 'mPP-0117', 'mPP-0091', 'mPP-0080', 'mPP-0268', 'mPP-0085',
              'mPP-0096', 'mPP-0130', 'mPP-0076', 'mPP-0178', 'mPP-0177', 'mPP-0120',
              'mPP-0114', 'mPP-0087', 'mPP-0181', 'mPP-0182', 'mPP-0182', 'mPP-0180',
              'mPP-0134', 'mPP-0129', 'mPP-0179', 'mPP-0133', 'mPP-0225', 'mPP-0093',
              'mPP-0090', 'mPP-0135', 'mPP-0089', 'mPP-0095', 'mPP-0226', 'mPP-0269',
              'mPP-0270', 'mPP-0084', 'mPP-0098', 'mPP-0092', 'mPP-0132', 'mPP-0185',
              'mPP-0184', 'mPP-0187', 'mPP-0186', 'mPP-0131', 'mPP-0183', 'mPP-0088',
              'mPP-0103', 'mPP-0112', 'mPP-0137', 'mPP-0138')

azu_all <- all[match(all$name.cmpd, azu_code, nomatch = 0), 'IC50']
Azu_ic1 <- replace(azu_all, is.na(azu_all), as.numeric(5000))
Azu_ic2 <- azu_all[!is.na(azu_all)]
oth_all <- all[-match(all$name.cmpd, azu_code, nomatch = 0), 'IC50']
oth_ic1 <- replace(oth_all, is.na(oth_all), as.numeric(5000))
oth_ic2 <- oth_all[!is.na(oth_all)]
df_oth2 <- data.frame(cbind(x = as.numeric(oth_ic2), y = 'oth'))
df_azu2 <- data.frame(cbind(x = as.numeric(Azu_ic2), y = 'azu'))
df2 <- (rbind(df_oth2, df_azu2))


ggplot(df3, aes(x = lx, colour = y, fill = y)) +
  # geom_histogram(binwidth=10)+
  geom_histogram(aes(y = ..density..), binwidth = 0.05, position = 'dodge2', alpha = 0.11) +
  geom_density(alpha = 0.1) +
  theme_bw()
ylim(0, 0.032) +
  # scale_x_continuous(trans = 'log2')+
  xlim(-10, 500)

df2.1 <- df2[df2$y == 'oth',]
ggplot(df2.1, aes(x = as.numeric(as.character(x)), colour = y, fill = y)) +
  # geom_histogram(binwidth=10)+
  geom_histogram(aes(y = ..density..), alpha = 0.01) +
  geom_density(alpha = 0.3) +
  theme_bw() +
  # scale_x_continuous(trans = 'log2')+
  xlim(-10, 700)
df2.2 <- df2[df2$y != 'oth',]
ggplot(df2.2, aes(x = as.numeric(as.character(x)), colour = y, fill = y)) +
  # geom_histogram(binwidth=10)+
  geom_histogram(position = "identity", aes(y = ..density..), alpha = 0.01) +
  geom_density(alpha = 0.3) +
  theme_bw() +
  # scale_x_continuous(trans = 'log2')+
  xlim(-10, 700)


df_oth1 <- data.frame(cbind(x = as.numeric(oth_ic1), y = 'oth'))
df_azu1 <- data.frame(cbind(x = as.numeric(Azu_ic1), y = 'azu'))
df1 <- (rbind(df_oth1, df_azu1))
ggplot(df1, aes(x = as.numeric(as.character(x)), colour = y)) +
  geom_density() +
  theme_bw()
xlim(-30, 200)


# #
#
# res.IC50.clust = clust.df.out(res.IC50)
# write.xlsx(res.IC50.clust,"Summary_IC50_030719.xlsx",row.names = F)
# update.summary.IC50("Summary_IC50_030719.xlsx",'Summary_IC50_030719.xlsx',IC50.new.cmpd.7)
# update.summary.IC50("Summary_IC50_030719.xlsx",'Summary_IC50_030719.xlsx',IC50.new.cmpd.8)
# update.summary.IC50("Summary_IC50_030719.xlsx",'Summary_IC50_030719.xlsx',IC50.new.cmpd.6)
# update.summary.IC50("Summary_IC50_030719.xlsx",'Summary_IC50_030719.xlsx',IC50.new.cmpd.5)
# update.summary.IC50("Summary_IC50_030719.xlsx",'Summary_IC50_030719.xlsx',IC50.new.cmpd.4)
# update.summary.IC50("Summary_IC50_030719.xlsx",'Summary_IC50_030719.xlsx',IC50.new.cmpd.3)
# update.summary.IC50("Summary_IC50_030719.xlsx",'Summary_IC50_250520.xlsx',IC50.new.cmpd.010716)
# tmp = update.summary.IC50("Summary_IC50_030719.xlsx",'Summary_IC50_030719.xlsx',IC50.new.cmpd)

## Excel to csv for maestro
setwd('Downloads/')
cmpd.tab <- read.csv('mPPaseCompoundRegister_to_Schrodinger(1).csv')
cmpd.tab <- cmpd.tab[0:405,]
colnames(cmpd.tab) <- c("Formula", 'SMILES', 'Generic code', 'Synthetic code', 'Isoxasole ms code', 'M (g/mol)',
                        'Exact mass (g/mol)', 'clogP', 'Mol. Top. PSA', 'Tested', 'Paper 1 code', 'm (mg)', 'Notes',
                        'Screening data exist', 'Screening 50 uM', 'Screening 20 uM', 'Screening 5 uM', 'Screening 1 uM',
                        'Toxicity', 'old IC50 (uM)', 'Supplier', 'Purity (% + method)', 'Identity (method)',
                        'Aggregation 50 uM', 'Aggregation 20 uM', "Aggregation 5 uM", "Aggregation 1 uM",
                        "new IC50 [CI95%] (uM)")
new_ic50 <- str_replace(cmpd.tab$`new IC50 [CI95%] (uM)`, '#N/A', '')
new_ic50 <- str_replace(new_ic50, 'NA  \\[ NA - NA ]', '')
new_ic50 <- str_replace(new_ic50, '.\nNA  \\[ NA - NA ]', '')

for (i in seq(1, length(new_ic50))) {

  if (str_detect(new_ic50[i], '\\[')) {
    new_ic50[i] <- (str_split(new_ic50[i], '  \\[')[[1]][1])

  }
}

cmpd.tab$`new IC50 [CI95%] (uM)` <- str_replace(cmpd.tab$`new IC50 [CI95%] (uM)`, '.\n', ' / ')
cmpd.tab$`Purity (% + method)` <- str_replace(cmpd.tab$`Purity (% + method)`, '�', ' >=')
cmpd.tab$`Purity (% + method)` <- str_replace(cmpd.tab$`Purity (% + method)`, '\n', ' / ')
cmpd.tab$`Identity (method)` <- str_replace(cmpd.tab$`Identity (method)`, '\n', ' / ')
# str_replace(cmpd.tab$new.IC50..CI95....uM.[2],'NA','')
# cmpd.tab$new.IC50..CI95....uM.[2]

cmpd.out <- cbind(cmpd.tab[, 1:27], 'IC50' = new_ic50, "new IC50 [CI95%] (uM)" = cmpd.tab[, 28])
write.csv(cmpd.out, '20200220_mPPaseCompoundRegister_clean.csv', quote = F, row.names = F, na = '', sep = '')

cmpd.tab$new.IC50..CI95....uM.[23]
str_replace(cmpd.tab, '.\n', ' / ')


# t4 = update.cmpd.register('Summary_IC50_260619.xlsx','Shared/mPPaseCompoundRegister_updated_14052019.xlsx','Shared/mPPaseCompoundRegister_updated_28062019.xlsx')
#
#
#
#
#
# t1 = data.frame(x=c(1,2,3,4,5,6),y=c(1,2,3,4,5,6))
#
# t2 = data.frame(l=c(3,5,6,41,4),x=c(1,2,3,5,6))
#
#
# merge(t2,t1,all.y = T)


# df.tot = rbind(df.ki,df.out)

# IC50.new.cmpd = plot.IC50(df.tot,plot_curve = T, save_pic = F)
# IC50.new.cmpd.clust =  clust.df.out(IC50.new.cmpd)
# write.xlsx(IC50.new.cmpd.clust,"Summary_IC50_090319_clust.xlsx",row.names = F)
# df.cmpd2604 = read.cmpd.test('IC50_determination_of_compounds.xlsx',nb_replic = 3,nb_id = 240)
# df.cmpd2205= read.cmpd.test('KV_inhibition assay (Sep-Nov2017).xlsx',nb_replic = 4,nb_id = 0)
#
# IC50.new.cmpd = plot.IC50(df.cmpd2205,plot_curve = T, save_pic = T)
# print(IC50.new.cmpd)
# # write.xlsx(IC50.new.cmpd.clust,"Summary_IC50_8points.xlsx",row.names = F)
#

# update.summary.IC50("Summary_IC50_060519.xlsx",'Summary_IC50_060519.xlsx',IC50.new.cmpd)

# plot.IC50(df.out)


# df.ki.2.0[which(df.ki.2.0$synth.name == 'TL4-63'),]


# X = split(df.ki.2,df.ki.2$id)
# nb.NA = 0
# for(i in 1:length(X)){
#
#   df.tmp = X[[i]]
#   np.tmp = nplr(df.tmp$x,df.tmp$y/100)
#   cmp.name =  df.tmp$gx.name[1]
#   IC.50 = getEstimates(np.tmp,targets = 0.5,conf.level = 0.9)
#   plot(np.tmp, pcol="grey40", lcol="skyblue1",
#        showEstim=.5, showInfl = F,showPoints = T, showCI = T,showGOF = FALSE,
#        cex.main=1.5, showSDerr = F,
#        lwd = 3,
#        main = paste("mPPase response to",cmp.name),
#        conf.level = .90,
#        unit = ' nM',
#        xlab=expression(Log[10](conc)),
#        ylab='mPPase activity')
#   # if (length(df.tmp)>4){
#   #   np =
#   #     print
# }
# }

# fit = nls(y~SSlogis(log(x),Asym,xmid,scal),data = data.frame(x=tmp.6.1$x,y=tmp.6$y))
# plot(tmp.6$x,tmp.6$y,log = 'x',xlim = c(0.1,500000),ylim=c(0,110))
#
# x= seq(0.1, 500000, length.out = 1000)
# y= predict(fit, newdata = list(x = seq(0.1, 500000, length.out = 1000)),interval = 'confidence')
# matlines(x,y,lwd=2)
#
# np1= nplr(tmp.6$x,tmp.6$y/100)
# plot(np1, pcol="grey40", lcol="skyblue1", showEstim=.5, showInfl=TRUE,
#      main="Default 'nplr' plot", cex.main=1.5)


# plot(npall, pcol="grey40", lcol="skyblue1",
#      showEstim=.5, showInfl = F,showPoints = T, showCI = T,showGOF = FALSE,
#      cex.main=1.5, showSDerr = F,
#      method="sdw",lwd = 3,
#      main = "Nppase. Response to NJ1-58",
#      conf.level = .90,
#      unit = ' nM',
#      xlab=expression(Log[10](conc)),
#      ylab='')


# library(ggplot2)
#
# ggplot(data=X[[100]][2:4,],aes(x=x,y=as.numeric(y),fill=x))+
#   geom_bar(stat = 'identity')+
#   geom_errorbar(aes(ymin=as.numeric(y)-as.numeric(stdev),
#                     ymax=as.numeric(y)+as.numeric(stdev)),
#                 position = position_dodge(19),
#                 width=.05
#                 )+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5),
#         axis.text = element_text(size = 11))+
#   theme(panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5))+
#   geom_text(stat='identity', aes(label=round(as.numeric(y), 3)), vjust=-1, size=3.5)


#
# for (i in seq(1:length(X))){
#   a = (nrow(X[[i]]))
#   if(a>4) print(X[[i]])
#   }
#
# id = c()
# nb=0
# l_std_0 = which(test$stdev==0)
# for (i in seq(2,length(l_std_0 ))){
#   nb= nb+1
#   l.nb= rep(nb,l_std_0[i]-l_std_0[i-1])
#   id = c(id,l.nb)
#   if (i == length(l_std_0 )){
#     nb = nb + 1
#     l.nb= rep(nb,nrow(test)-l_std_0[i]+1)
#     id = c(id,l.nb)
#   }
# }
# test_2=cbind(test,id)


#

# df.tmp = as.data.frame(not.summary.IC50)
# df.tmp[which(!is.na.data.frame(not.summary.IC50),T)]
#
#
#
# indices = data.frame(ind=which(!is.na.data.frame(not.summary.IC50),T))
# not.summary.IC50
#
# summary.IC50 = estimate.IC50[[1]][,0:5]
# summary.IC50[which(!is.na(summary.IC50$X__2)),]
# estimate.IC50.2 =read.xls("Desktop/mppase/preliminary_IC50s-2018.xlsx",'G1andG2',
#                         na.strings = c('','<NA>','NA',"#DIV/0!"))
# # estimate.IC50.2 = read_xlsx("Desktop/mppase/preliminary_IC50s-2018.xlsx",
#                             # col_names = F)
# # estimate.IC50[which(estimate.IC50 == '',T)] = NA
# summary.IC50.G1G2 = estimate.IC50.2[,seq(0,5)]
# summary.IC50.G1G2[which(!is.na(summary.IC50.G1G2[,2])),]
#
# estimate.IC50[-is.na(estimate.IC50)]
#
# data[rowSums(is.na(data)) != ncol(data),]
#
#
# src_path_2 = "/home/dreano/Desktop/mppase/Shared/mPPaseCompoundRegister_updated 23022018.xlsx"
#
# l_cmpd_name= data.frame(read_excel(src_path_2,2))
# test = apply(l_cmpd_name['Synthetic'],2, parse_cpd)
# df.test = cbind(l_cmpd_name, test)
# colnames(df.test) = c(colnames(l_cmpd_name),'l_name')
#
# l_cmpd_name[grep('C1/',df.test$l_name),'Generic']
#
# np = plot.IC50(df.cmpd2205,plot_curve = T, save_pic = T)
#
# n=5
# png('/home/dreano/Desktop/mppase/Picture_230519/curve/Test_4.png',
#     width=678*n, height=400*n, res=100*n, family="serif")
# plot(np, pcol="grey40", lcol="skyblue1",
#      showEstim = F, showInfl = F,
#      showPoints = T, showCI = T,
#      showGOF = F, cex.main=1.5,
#      showSDerr = F, lwd = 3,
#      main = 'p.title',
#      conf.level = 0.95,
#      unit = ' nM',
#      xlab = bquote(Log[10]~.('x.title')),
#      ylab = 'mPPase activity',ylim=c(0,1.2),xlim=c(-10,-4), cex.lab=1.25)
# text(labels = lgd.ic ,x=-9,y=0.1, font=2, col='skyblue4')
# 
# dev.off()


