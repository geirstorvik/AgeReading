fM = array(NA,c(4,Nsim,Y,A))
fM[1,,,] = fMbar[,1,,]
fM[2,,,] = fM1[,1,,]
fM[3,,,] = fM3[,1,,]
fM[4,,,] = fM6[,1,,]
Mlab = c("bar","1","3","6")
df = data.frame(N=as.vector(fM),year=NA,age=NA,sim=NA,M=as.factor(rep(Mlab,300)))
df$M = as.factor(df$M)
i=0
for(y in 1:Y)
  for(a in 1:A)
    for(m in 1:Nsim)
      for(M in 1:4)
    {
      i = i+1
      df$N[i] = abs(fM[M,m,y,a]-Truef[y,a])/Truef[y,a]
      df$M[i] = Mlab[m]
      df$sim[i] = m
      df$year[i] = y
      df$age[i] = a
      df$True[i] = Truef[y,a]
      }
df$M = as.factor(df$M)
df$age = as.factor(df$age)

df$year = as.factor(df$year)
df = within(df,M<-relevel(M,ref="bar"))
fit = lm(N~year+M,data=df)
