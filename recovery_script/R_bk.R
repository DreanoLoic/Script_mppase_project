# library("gdata", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
1561457053449:library('readxl')
1561457053577:library('gsubfn')
1561457053783:library('tidyr')
1561457053890:library("ChemmineR")
1561457054244:library('ggplot2')
1561457054253:library('nplr')
1561457054276:library('xlsx')
1561457054821:parse_cmpd = function(x){
1561457054821:out = paste(gsubfn(".", list(" " = "", "-" = "", "+" = "","("="/(" ),
1561457054821:x),'/',sep='')
1561457054821:return(out)
1561457054821:}
1561457054821:substramb <- function(x, n=4){
1561457054821:paste(substr(x,0,3),
1561457054821:substr(x, nchar(x)-n+1, nchar(x)),sep='')
1561457054821:}
1561457054822:gx.cmpd.name = function(synth.name,
1561457054822:df.list.cmpd=data.frame(
1561457054822:read_excel("/home/dreano/Desktop/mppase/Shared/mPPaseCompoundRegister_updated_13062019_LD_NJ.xlsx",2))){
1561457054822:# src_path = "/home/dreano/Desktop/mppase/Shared/mPPaseCompoundRegister_updated 23022018.xlsx"
1561457054822:# df.list.cmpd = data.frame(read_excel(src_path,2))
1561457054822:parse.synth.name = gsubfn(".", list(" " = "", "-" = "", "+" = "" ),synth.name)
1561457054822:l.gx.names = apply(df.list.cmpd['Synthetic'],2, parse_cmpd)
1561457054822:df.gx.names = cbind(df.list.cmpd, l.gx.names)
1561457054822:colnames(df.gx.names) = c(colnames(df.list.cmpd),'l.synth.name')
1561457054823:if (substr(synth.name,1,4)=='mPP-'){
1561457054823:synth.id = which(df.list.cmpd$Generic==synth.name)
1561457054823:gx.name = df.list.cmpd[synth.id,c('SMILES','Generic','Paper')]
1561457054823:}else if(substr(parse.synth.name,0,3)=='Amb'){
1561457054823:df.amb = data.frame(l.synth.name=
1561457054823:substramb(
1561457054823:as.character(
1561457054823:df.gx.names$l.synth.name)))
1561457054824:gx.name = df.list.cmpd[grep(paste(parse.synth.name,
1561457054824:'/',
1561457054824:sep=''),
1561457054824:df.amb$l.synth.name),
1561457054824:c(2,3,4)]
1561457054824:}else{
1561457054824:gx.name = df.list.cmpd[grep(paste(parse.synth.name,
1561457054825:'/',
1561457054825:sep=''),
1561457054825:df.gx.names$l.synth.name),
1561457054825:c(2,3,4)]
1561457054825:}
1561457054825:if (nrow(gx.name) == 0){
1561457054826:gx.name = data.frame(SMILES=NA,Generic=synth.name,Paper=NA)
1561457054826:}
1561457054826:return(gx.name)
1561457054826:}
1561457054827:parsing_IC50 = function(src_path = "~/Desktop/mppase/preliminary_IC50s-2018.xlsx"){
1561457054827:src_path_1 = "/home/dreano/Desktop/mppase/Shared/mPPaseCompoundRegister_updated_13062019_LD_NJ.xlsx"
1561457054827:df.list.cmpd = data.frame(read_excel(src_path_1,2))
1561457054827:l_sheets = excel_sheets(src_path)
1561457054827:estimate.IC50 = lapply(l_sheets,read_excel,
1561457054827:path=src_path,
1561457054827:col_names=FALSE)
1561457054827:df.out = data.frame()
1561457054827:nb_noname = 0
1561457054828:nb_id = 0
1561457054828:for ( i in 1:(length(l_sheets))) {
1561457054828:summary.IC50 = data.frame(estimate.IC50[[i]][,seq(1,5)])
1561457054828:not.summary.IC50 = data.frame(estimate.IC50[[i]][,-seq(1,5)])
1561457054828:id.mean = as.data.frame(which(not.summary.IC50=='mean',T) )
1561457054828:if (nrow(id.mean)==0) {
1561457054828:next
1561457054828:}
1561457054828:id.mean = id.mean[order(id.mean$row),]
1561457054829:for (j in 1:nrow(id.mean)){
1561457054829:name.cmpd = not.summary.IC50[id.mean$row[j]+2,id.mean$col[j]-8]
1561457054829:if(is.na(name.cmpd)| name.cmpd==1000 | name.cmpd == 100){
1561457054829:name.cmpd = not.summary.IC50[id.mean$row[j]+2,id.mean$col[j]-9]
1561457054829:if(is.na(name.cmpd)| name.cmpd==1000 | name.cmpd == 100){
1561457054829:name.cmpd = not.summary.IC50[id.mean$row[j]+2,id.mean$col[j]-10]
1561457054829:}
1561457054830:if(length(name.cmpd)==0) {
1561457054830:nb_noname = nb_noname + 1
1561457054830:next
1561457054830:}
1561457054830:if(is.na(name.cmpd) | name.cmpd == 0) {
1561457054831:nb_noname = nb_noname + 1
1561457054831:next
1561457054831:}
1561457054831:id.range = seq(max(id.mean$col[j]-10,0)
1561457054831:,id.mean$col[j])
1561457054831:id.col = which(not.summary.IC50[id.mean$row[j]+1,]==0 |
1561457054832:not.summary.IC50[id.mean$row[j]+1,]==0.1)
1561457054832:id.rm.1 = which(id.col > id.mean$col[j] |
1561457054832:id.col < id.mean$col[j] - 10)
1561457054832:if (length(id.rm.1) > 0){
1561457054832:id.col = id.col[-id.rm.1]
1561457054833:}
1561457054833:id.na = which(is.na(not.summary.IC50[,id.mean$col[j]]))
1561457054833:id.row = id.na[which(id.na>id.mean$row[j])]
1561457054833:id.seq.1 =  seq(id.mean$row[j]+1,id.na[which(id.na>id.mean$row[j])][1]-1)
1561457054833:x.axis = not.summary.IC50[id.seq.1,id.col]
1561457054834:x.axis = as.numeric(replace(x.axis,x.axis==0,0.1))
1561457054834:id.100= head(which(not.summary.IC50[id.mean$row[j]+1,]==100),-1)
1561457054834:id.rm  =  which(is.na(not.summary.IC50[id.mean$row[j]+2,id.100])|
1561457054834:id.100 >= id.mean$col[j] |
1561457054835:id.100 < id.mean$col[j] -10 )
1561457054835:if(length(id.rm)>0){
1561457054835:id.100 = id.100[-id.rm]
1561457054835:}
1561457054836:x.axis = rep(x.axis,length(id.100))
1561457054836:y.axis = as.numeric(gather(not.summary.IC50[id.seq.1,id.100])$value)
1561457054836:nb_id = nb_id + 1
1561457054836:df.translate = gx.cmpd.name(as.character(name.cmpd),df.list.cmpd)
1561457054837:df.tmp = data.frame(df.translate$Generic,
1561457054837:df.translate$Paper,
1561457054837:name.cmpd,
1561457054838:x.axis,
1561457054838:y.axis,
1561457054838:nb_id,
1561457054839:df.translate$SMILES)
1561457054839:colnames(df.tmp) = c('gx.name','paper.name','synth.name','x','y','id','SMILES')
1561457054839:df.out = rbind(df.out,df.tmp)
1561457054839:}
1561457054840:print(paste(nb_noname,'compound whith no name'))
1561457054841:return(df.out)
1561457054841:}
1561457054842:plot.IC50 = function(df.result = parsing_IC50(), conf.level = 0.9, showSDerr = F,
1561457054842:showInfl = F, showPoints = T, showCI = T, plot_curve=F,
1561457054842:showGOF = F, showEstim=F, save_pic = F, up.lim=5e6, n=5){
1561457054842:set.seed(123)
1561457054842:setwd('~/Desktop/mppase/')
1561457054842:df.return = data.frame()
1561457054842:l.result =  split(df.result,df.result$id)
1561457054842:for (i in 1:length(l.result)) {
1561457054842:df.cmpd = l.result[[i]]
1561457054842:df.cmpd$x = df.cmpd$x * 10^-9
1561457054843:np.cmpd = nplr(df.cmpd$x,df.cmpd$y/100) # Y axis scale from 0 to 1
1561457054843:name.cmpd = df.cmpd$gx.name[1]
1561457054843:paper.name = df.cmpd$paper.name[1]
1561457054843:synt.name = df.cmpd$synth.name[1]
1561457054843:smiles.cmpd = df.cmpd$SMILES[1]
1561457054843:id.cmpd = df.cmpd$id
1561457054843:if (is.na(paper.name)) {
1561457054844:p.title = paste("mPPase response to",name.cmpd)
1561457054844:x.title = paste('[ ',name.cmpd,' ] (M)',sep = '')
1561457054844:r.path= "./Picture_230519/"
1561457054844:}else{
1561457054844:p.title = paste("mPPase response to §",paper.name," ( ",name.cmpd, " )",
1561457054844:sep = '')
1561457054844:x.title = paste('[ §',paper.name,' ] (M)',sep='')
1561457054845:r.path= "./Picture_230519/Paper/"
1561457054845:}
1561457054845:IC.50 = getEstimates(np.cmpd,targets = 0.5,conf.level = conf.level)*10^9
1561457054845:lgd.ic = paste('IC50 =',
1561457054845:formatC(IC.50[[3]],
1561457054845:format = "e",
1561457054845:digits = 1),
1561457054846:'nM \n[',
1561457054846:formatC(IC.50[[2]],
1561457054846:format = "e",
1561457054846:digits = 1),
1561457054846:',',
1561457054847:formatC(IC.50[[4]],
1561457054847:format = "e",
1561457054847:digits = 1),']')
1561457054847:if(IC.50$y != 5e+08 | IC.50$x == Inf | IC.50$x < 1e-5 | IC.50$x > up.lim ){
1561457054847:IC.50 = data.frame(y=0,x.05=NA,x=NA,x.90=NA) # UPDATE >X ???
1561457054847:lgd.ic = paste('IC50 =',
1561457054848:'NA',
1561457054848:'\n[',
1561457054848:'NA',
1561457054848:',',
1561457054848:'NA',']')
1561457054849:}
1561457054849:df.IC50 = data.frame( name.cmpd, paper.name ,synt.name, IC.50$x, IC.50[2], IC.50[4], smiles.cmpd)
1561457054849:colnames(df.IC50) = c("name.cmpd", 'paper.name', 'synth.name', "IC50", "IC50.low", "IC50.up", 'SMILES')
1561457054849:if(plot_curve){
1561457054849:plot(np.cmpd, pcol="grey40", lcol="skyblue1",
1561457054850:showEstim = showEstim, showInfl = showInfl,
1561457054850:showPoints = showPoints, showCI = showCI,
1561457054850:showGOF = showGOF, cex.main=1.5,
1561457054850:showSDerr = showSDerr, lwd = 3,
1561457054851:main = p.title,
1561457054851:conf.level = conf.level,
1561457054851:unit = ' nM',
1561457054852:xlab = bquote(Log[10]~.(x.title)),
1561457054852:ylab = 'mPPase activity',ylim=c(0,1.2),xlim=c(-10,-4), cex.lab=1.2)
1561457054852:text(labels = lgd.ic ,x=-9,y=0.1, font=2, col='skyblue4')
1561457054852:if(save_pic){
1561457054853:png(paste(r.path,"curve/",'Curve_IC50_',name.cmpd,'_',synt.name,'_',id.cmpd,'.png',sep=''),
1561457054853:width=595.3684*n, height=400*n, res=100*n, family="serif")
1561457054853:plot(np.cmpd, pcol="grey40", lcol="skyblue1",
1561457054853:showEstim = showEstim, showInfl = showInfl,
1561457054854:showPoints = showPoints, showCI = showCI,
1561457054854:showGOF = showGOF, cex.main=1.5,
1561457054854:showSDerr = showSDerr, lwd = 3,
1561457054854:main = p.title,
1561457054855:conf.level = conf.level,
1561457054855:unit = ' nM',
1561457054855:xlab = bquote(Log[10]~.(x.title)),
1561457054856:ylab = 'mPPase activity',ylim=c(0,1.2),xlim=c(-10,-4), cex.lab=1.2)
1561457054856:text(labels = lgd.ic ,x=-9,y=0.1, font=2, col='skyblue4')
1561457054856:dev.off()
1561457054857:}
1561457054857:split.df.cmpd = split(df.cmpd[4:5],df.cmpd$x)
1561457054857:df.barplot = data.frame()
1561457054857:for (j in 1:length(split.df.cmpd)) {
1561457054858:x = round(log10(mean(split.df.cmpd[[j]]$x)), digits =3)
1561457054858:y = split.df.cmpd[[j]]$y/100
1561457054858:y = y[!is.na(y)]
1561457054859:y.mean = mean(y)
1561457054859:se = sd(y)/sqrt(length(y))
1561457054859:ci = 2*se
1561457054860:df.barplot= rbind(df.barplot,data.frame(x,y.mean,se,ci))
1561457054860:}
1561457054860:p2= ggplot(df.barplot,aes(x=as.factor(x),y=y.mean,fill=as.factor(x)))+
1561457054861:geom_bar(position=position_dodge(), stat="identity",colour="Grey70") +
1561457054861:geom_errorbar(aes(ymin=y.mean-ci, ymax=y.mean+ci), size=1.2,
1561457054862:width=.2,                    # Width of the error bar
1561457054862:position=position_dodge(0.9))+
1561457054862:theme_bw()+
1561457054863:theme(plot.title = element_text(hjust = 0.5, size=16, face="bold"),
1561457054863:axis.text.x = element_text(angle = 0, hjust = 0.5,vjust = 0.5),
1561457054863:axis.text = element_text(size = 14), legend.position = 'none',
1561457054864:axis.title=element_text(size=16),
1561457054864:panel.grid.minor = element_blank())+
1561457054865:ggtitle(p.title)+
1561457054865:xlab( bquote(Log[10]~.(x.title)))+
1561457054865:ylab('mPPase activity')+
1561457054866:scale_fill_brewer()+
1561457054866:annotate("Text",x=1,y=0.1,label=lgd.ic,
1561457054867:size=4.5, fontface=2, colour='skyblue4')
1561457054867:plot(p2)
1561457054867:if(save_pic){
1561457054868:ggsave(filename = paste(r.path,"barplot/",'Barplot_IC50_',name.cmpd,'_',synt.name,'_',id.cmpd,'.png',sep=''),
1561457054868:plot = p2,device='png', width=7.0729166667, height=4.75)}
1561457054869:}
1561457054869:df.return = rbind(df.return,df.IC50)
1561457054870:}
1561457054870:return(df.return[order(df.return$IC50),])
1561457054870:}
1561457054872:read.cmpd.test = function(src_path = "~/Desktop/mppase/New compounds_170119.xlsx",
1561457054872:nb_replic = 4,nb_id=0,start_sheet=1, end_sheet=0){
1561457054872:src_path_1 = "/home/dreano/Desktop/mppase/Shared/mPPaseCompoundRegister_updated_14052019.xlsx"
1561457054872:df.list.cmpd = data.frame(read_excel(src_path_1,2))
1561457054872:l_sheets = excel_sheets(src_path)
1561457054872:cmpd.data = lapply(l_sheets,read_excel,
1561457054872:path=src_path,
1561457054872:col_names=FALSE)
1561457054872:df.out = data.frame()
1561457054872:df.ini = data.frame(gx.name=as.character(),
1561457054873:paper.name=as.character(),
1561457054873:synth.name=as.character(),
1561457054873:x=as.numeric(),
1561457054873:y=as.numeric(),
1561457054873:SMILES=as.character())
1561457054873:nb_noname = 0
1561457054873:if (end_sheet ==0) end_sheet = length(l_sheets)
1561457054873:for ( i in start_sheet:end_sheet) {
1561457054874:df.full = data.frame(cmpd.data[[i]])
1561457054874:act.id = data.frame(which(df.full=='% Activity',T))
1561457054874:if (nrow(act.id)==0) next
1561457054874:row.id = seq(act.id$row+2,act.id$row+9)
1561457054874:plate.scheme = df.full[row.id,seq(2,12,nb_replic)]
1561457054874:col.id = seq(act.id$col+1,act.id$col+12)
1561457054875:inib.data = df.full[row.id,col.id]
1561457054875:for(r in seq(1,nrow(inib.data))){
1561457054875:for (c in seq(nb_replic+1,ncol(inib.data))) {
1561457054875:exp.name = plate.scheme[r,ceiling(c/nb_replic)]
1561457054876:if(is.na(exp.name)) next()
1561457054876:exp.info = strsplit(exp.name,' ')[[1]]
1561457054876:if(is.na(exp.info[1])) next()
1561457054876:name.cmpd = paste(exp.info[1:length(exp.info)-1], collapse = ' ')
1561457054876:y= inib.data[r,c]
1561457054877:x = as.numeric(gsubfn(".", list("u" = "", "M" = "", "n" = "" ),
1561457054877:exp.info[length(exp.info)])) * 1000 # convert to nM #
1561457054877:df.translate = gx.cmpd.name(as.character(name.cmpd),df.list.cmpd)
1561457054878:df.tmp = data.frame(df.translate$Generic,
1561457054878:df.translate$Paper,
1561457054878:name.cmpd,
1561457054878:x,
1561457054878:as.numeric(y),
1561457054879:df.translate$SMILES)
1561457054879:colnames(df.tmp) = c('gx.name','paper.name','synth.name','x','y','SMILES')
1561457054879:df.ini = rbind(df.ini,df.tmp)
1561457054879:}
1561457054880:df.split = split(df.ini,df.ini$synth.name)
1561457054881:for (l in seq(1,length(df.split))) {
1561457054881:df.split[[l]] = rbind(data.frame(df.split[[l]][1,c(1,2,3,6)],x=0.1,y=100),
1561457054881:df.split[[l]])
1561457054881:df.split[[l]] = cbind(df.split[[l]],id=nb_id+l)
1561457054882:df.out = rbind(df.out,df.split[[l]][,c(1,2,3,5,6,7,4)])
1561457054882:}
1561457054882:return(df.out)
1561457054882:}
1561457054883:clust.df.out = function(df.no.clust){
1561457054883:rm.name = c()
1561457054883:df.no.clust = df.no.clust[order(df.no.clust$IC50),]
1561457054883:df.out = data.frame()
1561457054883:for (name in df.no.clust$name.cmpd) {
1561457054883:if (name %in% rm.name) {
1561457054884:next()
1561457054884:}
1561457054884:rm.name = cbind(rm.name,name)
1561457054884:df.out = rbind(df.out,df.no.clust[which(df.no.clust$name.cmpd== name),])
1561457054884:}
1561457054884:return(df.out)
1561457054884:}
1561457054884:update.summary.IC50 = function(file_src, file_out, df.new){
1561457054884:df.old = read_excel(file_src,1)
1561457054885:df.update = rbind(df.old, df.new)
1561457054885:df.out = as.data.frame(clust.df.out(df.update))
1561457054885:write.xlsx(df.out,file_out,row.names = F)
1561457054885:return(df.out)
1561457054885:}
1561457054885:update.cmpd.register = function(ic50_sum,cpd_reg,file_out){
1561457054885:df.sum= read_excel(ic50_sum,1)
1561457054885:df.register= read_excel(cpd_reg,2)
1561457054885:instr= read_excel(cpd_reg,1)
1561457054885:df.1 = data.frame(Generic = df.sum$name.cmpd, new_Ic50 = round(df.sum$IC50*(10^-3),3))
1561457054886:df.conc1 = mutate(group_by(df.1,Generic),new_Ic50 = paste0(new_Ic50, collapse = ", "))
1561457054886:df.1.clean = df.conc1[!duplicated.data.frame(df.conc1),]
1561457054886:df.2 = merge(df.register,df.1.clean)
1561457054886:write.xlsx(instr,file_out, sheetName = 'Instructions')
1561457054886:write.xlsx(df.2,file_out,row.names = T,append = T, sheetName = 'CompoundRegister')
1561457054886:return(df.2)
1561457054886:}
1561457054887:update.cmpd.register('Summary_IC50_120619.xlsx','Shared/mPPaseCompoundRegister_updated_14052019.xlsx','Shared/mPPaseCompoundRegister_updated_13062019.xlsx')
1561457054924:df.ki = parsing_IC50()
1561457102905:res.IC50= plot.IC50(df.ki,plot_curve = F, save_pic = F)
1561457113753:# write.table(res.IC50,"Summary_IC50_240519.csv",col.names = T,row.names = F)
1561457113753:res.IC50.clust = clust.df.out(res.IC50)
1561457114058:write.xlsx(res.IC50.clust,"Summary_IC50_240519.xlsx",row.names = F)
1561457115420:# write.table(res.IC50.clust,"Summary_IC50_clust.csv",col.names = T,row.names = F)
1561457115420:pics = F
1561457115420:# df.out = read.cmpd.test(nb_id = 215)
1561457115421:df.out = read.cmpd.test('IC50_determination_of_compounds.xlsx',nb_replic = 3,nb_id = 240)
1561457160712:IC50.new.cmpd = plot.IC50(df.out,plot_curve = pics, save_pic = pics)
1561457161177:df.cmpd2205= read.cmpd.test('KV_inhibition assay (Sep-Nov2017).xlsx',nb_replic = 4,nb_id = 0)
1561457238346:IC50.new.cmpd.2 = plot.IC50(df.cmpd2205,plot_curve = pics, save_pic = pics)
1561457239755:df.cmpd1701= read.cmpd.test('New compounds_170119.xlsx',nb_replic = 4,nb_id = 0)
1561457289002:IC50.new.cmpd.3 = plot.IC50(df.cmpd1701,plot_curve = pics, save_pic = pics)
1561457289713:df.cmpd0601.1= read.cmpd.test('Aaron_Inhibition test and IC50_310118.xlsx',nb_replic = 4,nb_id = 0,end_sheet = 1)
1561457303103:IC50.new.cmpd.4 = plot.IC50(df.cmpd0601.1,plot_curve = pics, save_pic = pics)
1561457303201:df.cmpd0601.2=read.cmpd.test('Aaron_Inhibition test and IC50_310118.xlsx',nb_replic = 3,nb_id = 0,start_sheet = 2)
1561457318229:IC50.new.cmpd.5 = plot.IC50(df.cmpd0601.2,plot_curve = pics, save_pic = pics)
1561457318380:df.cmpd250418=read.cmpd.test('Inhibition_250418_result.xlsx',nb_replic = 4,nb_id = 0)
1561457328549:IC50.new.cmpd.6 = plot.IC50(df.cmpd250418,plot_curve = pics, save_pic =pics)
1561457328706:df.cmpd1701_= read.cmpd.test('New compounds_170119_.xlsx',nb_replic = 4,nb_id = 0 , start_sheet =  5)
1561457346605:IC50.new.cmpd.7 = plot.IC50(df.cmpd1701_,plot_curve = pics, save_pic = pics)
1561457346835:df.cmpd130619=read.cmpd.test('AKI samples mPPase Spring 2019.xlsx',nb_replic = 4,nb_id = 0)
1561457395980:IC50.new.cmpd.8 = plot.IC50(df.cmpd130619,plot_curve = pics, save_pic = pics)
