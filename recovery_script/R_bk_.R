ta", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
library('readxl')
library('gsubfn')
library('tidyr')
library("ChemmineR")
library('ggplot2')
library('nplr')
library('xlsx')
parse_cmpd = function(x){
out = paste(gsubfn(".", list(" " = "", "-" = "", "+" = "","("="/(" ),
x),'/',sep='')
return(out)
}
substramb <- function(x, n=4){
paste(substr(x,0,3),
substr(x, nchar(x)-n+1, nchar(x)),sep='')
}
gx.cmpd.name = function(synth.name,
df.list.cmpd=data.frame(
read_excel("/home/dreano/Desktop/mppase/Shared/mPPaseCompoundRegister_updated_13062019_LD_NJ.xlsx",2))){
# src_path = "/home/dreano/Desktop/mppase/Shared/mPPaseCompoundRegister_updated 23022018.xlsx"
# df.list.cmpd = data.frame(read_excel(src_path,2))
parse.synth.name = gsubfn(".", list(" " = "", "-" = "", "+" = "" ),synth.name)
l.gx.names = apply(df.list.cmpd['Synthetic'],2, parse_cmpd)
df.gx.names = cbind(df.list.cmpd, l.gx.names)
colnames(df.gx.names) = c(colnames(df.list.cmpd),'l.synth.name')
if (substr(synth.name,1,4)=='mPP-'){
synth.id = which(df.list.cmpd$Generic==synth.name)
gx.name = df.list.cmpd[synth.id,c('SMILES','Generic','Paper')]
}else if(substr(parse.synth.name,0,3)=='Amb'){
df.amb = data.frame(l.synth.name=
substramb(
as.character(
df.gx.names$l.synth.name)))
gx.name = df.list.cmpd[grep(paste(parse.synth.name,
'/',
sep=''),
df.amb$l.synth.name),
c(2,3,4)]
}else{
gx.name = df.list.cmpd[grep(paste(parse.synth.name,
'/',
sep=''),
df.gx.names$l.synth.name),
c(2,3,4)]
}
if (nrow(gx.name) == 0){
gx.name = data.frame(SMILES=NA,Generic=synth.name,Paper=NA)
}
return(gx.name)
}
parsing_IC50 = function(src_path = "~/Desktop/mppase/preliminary_IC50s-2018.xlsx"){
src_path_1 = "/home/dreano/Desktop/mppase/Shared/mPPaseCompoundRegister_updated_13062019_LD_NJ.xlsx"
df.list.cmpd = data.frame(read_excel(src_path_1,2))
l_sheets = excel_sheets(src_path)
estimate.IC50 = lapply(l_sheets,read_excel,
path=src_path,
col_names=FALSE)
df.out = data.frame()
nb_noname = 0
nb_id = 0
for ( i in 1:(length(l_sheets))) {
summary.IC50 = data.frame(estimate.IC50[[i]][,seq(1,5)])
not.summary.IC50 = data.frame(estimate.IC50[[i]][,-seq(1,5)])
id.mean = as.data.frame(which(not.summary.IC50=='mean',T) )
if (nrow(id.mean)==0) {
next
}
id.mean = id.mean[order(id.mean$row),]
for (j in 1:nrow(id.mean)){
name.cmpd = not.summary.IC50[id.mean$row[j]+2,id.mean$col[j]-8]
if(is.na(name.cmpd)| name.cmpd==1000 | name.cmpd == 100){
name.cmpd = not.summary.IC50[id.mean$row[j]+2,id.mean$col[j]-9]
if(is.na(name.cmpd)| name.cmpd==1000 | name.cmpd == 100){
name.cmpd = not.summary.IC50[id.mean$row[j]+2,id.mean$col[j]-10]
}
if(length(name.cmpd)==0) {
nb_noname = nb_noname + 1
next
}
if(is.na(name.cmpd) | name.cmpd == 0) {
nb_noname = nb_noname + 1
next
}
id.range = seq(max(id.mean$col[j]-10,0)
,id.mean$col[j])
id.col = which(not.summary.IC50[id.mean$row[j]+1,]==0 |
not.summary.IC50[id.mean$row[j]+1,]==0.1)
id.rm.1 = which(id.col > id.mean$col[j] |
id.col < id.mean$col[j] - 10)
if (length(id.rm.1) > 0){
id.col = id.col[-id.rm.1]
}
id.na = which(is.na(not.summary.IC50[,id.mean$col[j]]))
id.row = id.na[which(id.na>id.mean$row[j])]
id.seq.1 =  seq(id.mean$row[j]+1,id.na[which(id.na>id.mean$row[j])][1]-1)
x.axis = not.summary.IC50[id.seq.1,id.col]
x.axis = as.numeric(replace(x.axis,x.axis==0,0.1))
id.100= head(which(not.summary.IC50[id.mean$row[j]+1,]==100),-1)
id.rm  =  which(is.na(not.summary.IC50[id.mean$row[j]+2,id.100])|
id.100 >= id.mean$col[j] |
id.100 < id.mean$col[j] -10 )
if(length(id.rm)>0){
id.100 = id.100[-id.rm]
}
x.axis = rep(x.axis,length(id.100))
y.axis = as.numeric(gather(not.summary.IC50[id.seq.1,id.100])$value)
nb_id = nb_id + 1
df.translate = gx.cmpd.name(as.character(name.cmpd),df.list.cmpd)
df.tmp = data.frame(df.translate$Generic,
df.translate$Paper,
name.cmpd,
x.axis,
y.axis,
nb_id,
df.translate$SMILES)
colnames(df.tmp) = c('gx.name','paper.name','synth.name','x','y','id','SMILES')
df.out = rbind(df.out,df.tmp)
}
print(paste(nb_noname,'compound whith no name'))
return(df.out)
}
plot.IC50 = function(df.result = parsing_IC50(), conf.level = 0.9, showSDerr = F,
showInfl = F, showPoints = T, showCI = T, plot_curve=F,
showGOF = F, showEstim=F, save_pic = F, up.lim=5e6, n=5){
set.seed(123)
setwd('~/Desktop/mppase/')
df.return = data.frame()
l.result =  split(df.result,df.result$id)
for (i in 1:length(l.result)) {
df.cmpd = l.result[[i]]
df.cmpd$x = df.cmpd$x * 10^-9
np.cmpd = nplr(df.cmpd$x,df.cmpd$y/100) # Y axis scale from 0 to 1
name.cmpd = df.cmpd$gx.name[1]
paper.name = df.cmpd$paper.name[1]
synt.name = df.cmpd$synth.name[1]
smiles.cmpd = df.cmpd$SMILES[1]
id.cmpd = df.cmpd$id
if (is.na(paper.name)) {
p.title = paste("mPPase response to",name.cmpd)
x.title = paste('[ ',name.cmpd,' ] (M)',sep = '')
r.path= "./Picture_230519/"
}else{
p.title = paste("mPPase response to §",paper.name," ( ",name.cmpd, " )",
sep = '')
x.title = paste('[ §',paper.name,' ] (M)',sep='')
r.path= "./Picture_230519/Paper/"
}
IC.50 = getEstimates(np.cmpd,targets = 0.5,conf.level = conf.level)*10^9
lgd.ic = paste('IC50 =',
formatC(IC.50[[3]],
format = "e",
digits = 1),
'nM \n[',
formatC(IC.50[[2]],
format = "e",
digits = 1),
',',
formatC(IC.50[[4]],
format = "e",
digits = 1),']')
if(IC.50$y != 5e+08 | IC.50$x == Inf | IC.50$x < 1e-5 | IC.50$x > up.lim ){
IC.50 = data.frame(y=0,x.05=NA,x=NA,x.90=NA) # UPDATE >X ???
lgd.ic = paste('IC50 =',
'NA',
'\n[',
'NA',
',',
'NA',']')
}
df.IC50 = data.frame( name.cmpd, paper.name ,synt.name, IC.50$x, IC.50[2], IC.50[4], smiles.cmpd)
colnames(df.IC50) = c("name.cmpd", 'paper.name', 'synth.name', "IC50", "IC50.low", "IC50.up", 'SMILES')
if(plot_curve){
plot(np.cmpd, pcol="grey40", lcol="skyblue1",
showEstim = showEstim, showInfl = showInfl,
showPoints = showPoints, showCI = showCI,
showGOF = showGOF, cex.main=1.5,
showSDerr = showSDerr, lwd = 3,
main = p.title,
conf.level = conf.level,
unit = ' nM',
xlab = bquote(Log[10]~.(x.title)),
ylab = 'mPPase activity',ylim=c(0,1.2),xlim=c(-10,-4), cex.lab=1.2)
text(labels = lgd.ic ,x=-9,y=0.1, font=2, col='skyblue4')
if(save_pic){
png(paste(r.path,"curve/",'Curve_IC50_',name.cmpd,'_',synt.name,'_',id.cmpd,'.png',sep=''),
width=595.3684*n, height=400*n, res=100*n, family="serif")
plot(np.cmpd, pcol="grey40", lcol="skyblue1",
showEstim = showEstim, showInfl = showInfl,
showPoints = showPoints, showCI = showCI,
showGOF = showGOF, cex.main=1.5,
showSDerr = showSDerr, lwd = 3,
main = p.title,
conf.level = conf.level,
unit = ' nM',
xlab = bquote(Log[10]~.(x.title)),
ylab = 'mPPase activity',ylim=c(0,1.2),xlim=c(-10,-4), cex.lab=1.2)
text(labels = lgd.ic ,x=-9,y=0.1, font=2, col='skyblue4')
dev.off()
}
split.df.cmpd = split(df.cmpd[4:5],df.cmpd$x)
df.barplot = data.frame()
for (j in 1:length(split.df.cmpd)) {
x = round(log10(mean(split.df.cmpd[[j]]$x)), digits =3)
y = split.df.cmpd[[j]]$y/100
y = y[!is.na(y)]
y.mean = mean(y)
se = sd(y)/sqrt(length(y))
ci = 2*se
df.barplot= rbind(df.barplot,data.frame(x,y.mean,se,ci))
}
p2= ggplot(df.barplot,aes(x=as.factor(x),y=y.mean,fill=as.factor(x)))+
geom_bar(position=position_dodge(), stat="identity",colour="Grey70") +
geom_errorbar(aes(ymin=y.mean-ci, ymax=y.mean+ci), size=1.2,
width=.2,                    # Width of the error bar
position=position_dodge(0.9))+
theme_bw()+
theme(plot.title = element_text(hjust = 0.5, size=16, face="bold"),
axis.text.x = element_text(angle = 0, hjust = 0.5,vjust = 0.5),
axis.text = element_text(size = 14), legend.position = 'none',
axis.title=element_text(size=16),
panel.grid.minor = element_blank())+
ggtitle(p.title)+
xlab( bquote(Log[10]~.(x.title)))+
ylab('mPPase activity')+
scale_fill_brewer()+
annotate("Text",x=1,y=0.1,label=lgd.ic,
size=4.5, fontface=2, colour='skyblue4')
plot(p2)
if(save_pic){
ggsave(filename = paste(r.path,"barplot/",'Barplot_IC50_',name.cmpd,'_',synt.name,'_',id.cmpd,'.png',sep=''),
plot = p2,device='png', width=7.0729166667, height=4.75)}
}
df.return = rbind(df.return,df.IC50)
}
return(df.return[order(df.return$IC50),])
}
read.cmpd.test = function(src_path = "~/Desktop/mppase/New compounds_170119.xlsx",
nb_replic = 4,nb_id=0,start_sheet=1, end_sheet=0){
src_path_1 = "/home/dreano/Desktop/mppase/Shared/mPPaseCompoundRegister_updated_14052019.xlsx"
df.list.cmpd = data.frame(read_excel(src_path_1,2))
l_sheets = excel_sheets(src_path)
cmpd.data = lapply(l_sheets,read_excel,
path=src_path,
col_names=FALSE)
df.out = data.frame()
df.ini = data.frame(gx.name=as.character(),
paper.name=as.character(),
synth.name=as.character(),
x=as.numeric(),
y=as.numeric(),
SMILES=as.character())
nb_noname = 0
if (end_sheet ==0) end_sheet = length(l_sheets)
for ( i in start_sheet:end_sheet) {
df.full = data.frame(cmpd.data[[i]])
act.id = data.frame(which(df.full=='% Activity',T))
if (nrow(act.id)==0) next
row.id = seq(act.id$row+2,act.id$row+9)
plate.scheme = df.full[row.id,seq(2,12,nb_replic)]
col.id = seq(act.id$col+1,act.id$col+12)
inib.data = df.full[row.id,col.id]
for(r in seq(1,nrow(inib.data))){
for (c in seq(nb_replic+1,ncol(inib.data))) {
exp.name = plate.scheme[r,ceiling(c/nb_replic)]
if(is.na(exp.name)) next()
exp.info = strsplit(exp.name,' ')[[1]]
if(is.na(exp.info[1])) next()
name.cmpd = paste(exp.info[1:length(exp.info)-1], collapse = ' ')
y= inib.data[r,c]
x = as.numeric(gsubfn(".", list("u" = "", "M" = "", "n" = "" ),
exp.info[length(exp.info)])) * 1000 # convert to nM #
df.translate = gx.cmpd.name(as.character(name.cmpd),df.list.cmpd)
df.tmp = data.frame(df.translate$Generic,
df.translate$Paper,
name.cmpd,
x,
as.numeric(y),
df.translate$SMILES)
colnames(df.tmp) = c('gx.name','paper.name','synth.name','x','y','SMILES')
df.ini = rbind(df.ini,df.tmp)
}
df.split = split(df.ini,df.ini$synth.name)
for (l in seq(1,length(df.split))) {
df.split[[l]] = rbind(data.frame(df.split[[l]][1,c(1,2,3,6)],x=0.1,y=100),
df.split[[l]])
df.split[[l]] = cbind(df.split[[l]],id=nb_id+l)
df.out = rbind(df.out,df.split[[l]][,c(1,2,3,5,6,7,4)])
}
return(df.out)
}
clust.df.out = function(df.no.clust){
rm.name = c()
df.no.clust = df.no.clust[order(df.no.clust$IC50),]
df.out = data.frame()
for (name in df.no.clust$name.cmpd) {
if (name %in% rm.name) {
next()
}
rm.name = cbind(rm.name,name)
df.out = rbind(df.out,df.no.clust[which(df.no.clust$name.cmpd== name),])
}
return(df.out)
}
update.summary.IC50 = function(file_src, file_out, df.new){
df.old = read_excel(file_src,1)
df.update = rbind(df.old, df.new)
df.out = as.data.frame(clust.df.out(df.update))
write.xlsx(df.out,file_out,row.names = F)
return(df.out)
}
update.cmpd.register = function(ic50_sum,cpd_reg,file_out){
df.sum= read_excel(ic50_sum,1)
df.register= read_excel(cpd_reg,2)
instr= read_excel(cpd_reg,1)
df.1 = data.frame(Generic = df.sum$name.cmpd, new_Ic50 = round(df.sum$IC50*(10^-3),3))
df.conc1 = mutate(group_by(df.1,Generic),new_Ic50 = paste0(new_Ic50, collapse = ", "))
df.1.clean = df.conc1[!duplicated.data.frame(df.conc1),]
df.2 = merge(df.register,df.1.clean)
write.xlsx(instr,file_out, sheetName = 'Instructions')
write.xlsx(df.2,file_out,row.names = T,append = T, sheetName = 'CompoundRegister')
return(df.2)
}
update.cmpd.register('Summary_IC50_120619.xlsx','Shared/mPPaseCompoundRegister_updated_14052019.xlsx','Shared/mPPaseCompoundRegister_updated_13062019.xlsx')
df.ki = parsing_IC50()
res.IC50= plot.IC50(df.ki,plot_curve = F, save_pic = F)
# write.table(res.IC50,"Summary_IC50_240519.csv",col.names = T,row.names = F)
res.IC50.clust = clust.df.out(res.IC50)
write.xlsx(res.IC50.clust,"Summary_IC50_240519.xlsx",row.names = F)
# write.table(res.IC50.clust,"Summary_IC50_clust.csv",col.names = T,row.names = F)
pics = F
# df.out = read.cmpd.test(nb_id = 215)
df.out = read.cmpd.test('IC50_determination_of_compounds.xlsx',nb_replic = 3,nb_id = 240)
IC50.new.cmpd = plot.IC50(df.out,plot_curve = pics, save_pic = pics)
df.cmpd2205= read.cmpd.test('KV_inhibition assay (Sep-Nov2017).xlsx',nb_replic = 4,nb_id = 0)
IC50.new.cmpd.2 = plot.IC50(df.cmpd2205,plot_curve = pics, save_pic = pics)
df.cmpd1701= read.cmpd.test('New compounds_170119.xlsx',nb_replic = 4,nb_id = 0)
IC50.new.cmpd.3 = plot.IC50(df.cmpd1701,plot_curve = pics, save_pic = pics)
df.cmpd0601.1= read.cmpd.test('Aaron_Inhibition test and IC50_310118.xlsx',nb_replic = 4,nb_id = 0,end_sheet = 1)
IC50.new.cmpd.4 = plot.IC50(df.cmpd0601.1,plot_curve = pics, save_pic = pics)
df.cmpd0601.2=read.cmpd.test('Aaron_Inhibition test and IC50_310118.xlsx',nb_replic = 3,nb_id = 0,start_sheet = 2)
IC50.new.cmpd.5 = plot.IC50(df.cmpd0601.2,plot_curve = pics, save_pic = pics)
df.cmpd250418=read.cmpd.test('Inhibition_250418_result.xlsx',nb_replic = 4,nb_id = 0)
IC50.new.cmpd.6 = plot.IC50(df.cmpd250418,plot_curve = pics, save_pic =pics)
df.cmpd1701_= read.cmpd.test('New compounds_170119_.xlsx',nb_replic = 4,nb_id = 0 , start_sheet =  5)
IC50.new.cmpd.7 = plot.IC50(df.cmpd1701_,plot_curve = pics, save_pic = pics)
df.cmpd130619=read.cmpd.test('AKI samples mPPase Spring 2019.xlsx',nb_replic = 4,nb_id = 0)
IC50.new.cmpd.8 = plot.IC50(df.cmpd130619,plot_curve = pics, save_pic = pics)
