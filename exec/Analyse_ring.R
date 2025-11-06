## Analyse ring data
foo = table(d.comb$unitstationindex,d.comb$unitreaderindex)
for(i in 1:8)
{
  foo2 = foo[foo[,i]>0,]
  show(table(rowSums(foo2)))
}
