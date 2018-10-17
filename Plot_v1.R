par(mfrow =  c(2,2))  

plot(t_daily,type = 'l')

plot(t1_sluice,type = 'l')

plot(SM1[,2],type = 'l',ylim = c(0,200))
lines(SM1[,3],col='green')
lines(SM2[,5])

plot(ET1[,3],type = 'l',col='green',ylim = c(0,7))
lines(ET1[,2])
lines(ET2[,5])
lines(PET_max,col='red')

plot(t1_stage[12:110],type = 'l',ylim=c(0,4))
lines(t1_field_stage$T1_stage.m.[12:110],col = 'purple')

NSE.default(t1_stage[12:110],t1_field_stage$T1_stage.m.[12:110])
dev.off()
plot(AQ1,col='green',type = 'l',ylim=c(0,400))
lines(AQ2)

plot(AQ2,type = 'l',ylim=c(0,400))

plot(GW1_used[,3],type = 'l')

AQ1old = AQ1

#Water Balance terms in terms of volumes
p1_vol=sum(t_daily)*(1/1000)*(t1_const$max_catch)
#These numbers are little off because the Fallow HRU size reduces as tank size increases
ro1_sum=apply(RF1,2,sum)
ro1_vol=sum(ro1_sum*(1/1000)*HRU_Areas1)

IF1_sum=apply(IF1,2,sum)
IF1_vol=sum(IF1_sum*(1/1000)*HRU_Areas1)

ET1_sum=apply(ET1,2,sum)
ET1_vol=sum(ET1_sum*(1/1000)*HRU_Areas1)

DP1_sum=apply(DP1,2,sum)
DP1_vol=sum(DP1_sum*(1/1000)*HRU_Areas1)

GW1_vol=sum(GW1_used)*(1/1000)*t1_const$max_catch

DP1_vol+ET1_vol+IF1_vol+ro1_vol
p1_vol+GW1_vol
sum(t1_ET)+sum(t1_GW)+sum(t1_sluice)+sum(t1_spill)

t1_stagedf = as.data.frame(t1_stage[12:110])
t1_stagedf$act = t1_field_stage$T1_stage.m.[13:111]
t1_stagedf$dates = as.Date(index(t_daily[12:110]))
t1_stagedf$rain = coredata(t_daily[12:110])
colnames(t1_stagedf)[1]='model'

c5=ggplot()+
  #geom_point(data= db_merge, aes(x = db_merge$Sl_x, y = db_merge$Sl_y,colour=factor(db_merge$region),shape = factor(db_merge$region)),stroke=1.4,size = 2)+
  geom_line(data= t1_stagedf, aes(x = dates, y = model,col='a'),size = 1)+
  geom_line(data= t1_stagedf, aes(x = dates, y = act,col='b'),size = 1)+
  #geom_point(data= new_RJ, aes(x = Join_Out_4, y = out_stnm_5),size = 1.5)+
  annotate(geom="text", label="Scatter plot",
           color="red")+
  #geom_text(x = 1, y = 1, label = 'NSE = 0.8') +
  theme(axis.text=element_text(size=16,colour = 'black'),
        axis.title=element_text(size=16, face='bold'),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        plot.title=element_text(size=20,face = "bold"),
        legend.position="right",
        legend.text=element_text(size=18),legend.title=element_blank(),legend.key.size = unit(1, "cm"),
        panel.background = element_blank(),axis.line = element_line(size = 1, colour = "black")) +
  scale_color_manual(name="",labels=c('Model','Measured'),values=c('tomato4','darkslateblue'))+
  #scale_shape_manual(name="",labels=c('Decreasing','Increasing'),values=c(16,17))+
  labs(x='Time', y='Tank Stage (m)',title='Tank Stage: Model vs Measured')

#ggsave('Stage_oct1.png',plot=c5,scale=1, dpi=300,width =16.73,height = 12.7,units = "cm")

c6=ggplot()+
  #geom_point(data= db_merge, aes(x = db_merge$Sl_x, y = db_merge$Sl_y,colour=factor(db_merge$region),shape = factor(db_merge$region)),stroke=1.4,size = 2)+
  geom_bar(data= t1_stagedf, aes(x = dates, y = rain),stat="identity",size = 1)+
  theme(axis.text=element_text(size=16,colour = 'black'),
        axis.title=element_text(size=16, face='bold'),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        plot.title=element_text(size=20,face = "bold"),
        legend.position="bottom",
        legend.text=element_text(size=18),legend.title=element_blank(),legend.key.size = unit(1, "cm"),
        panel.background = element_blank(),axis.line = element_line(size = 1, colour = "black")) +
  #scale_color_manual(name="",labels=c('Model','Measured'),values=c('tomato4','darkslateblue'))+
  #scale_shape_manual(name="",labels=c('Decreasing','Increasing'),values=c(16,17))+
  labs(x='Time', y='Rainfall (mm)',title='Rainfall (mm)')
#ggsave('rain_oct1.png',plot=c6,scale=1, dpi=300,width =16.73,height = 12.7,units = "cm")

