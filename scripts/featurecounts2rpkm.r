a=read.delim(file("stdin"),skip=1)
a[,7] = a[,7]/sum(a[,7])*10**6 /a$Length *10**3
write.table(a,stdout(),sep='\t',quote=F,row.names=F)
