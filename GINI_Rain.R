install.packages('ineq')
require(ineq)

########################## Long term rainfall
input_file = "Tirumungulam_rainfall_modified_v2.csv"
input_file = "Tirumungulam_rainfall.csv"

input = read.csv(input_file,header=TRUE)
in_date = as.Date(input[,1],"%Y-%m-%d")
x = zoo(input[,2],in_date)
t = window(x,start=as.Date(input[1,1]))

x2=subset(t, t>0)
years = unique(year(index(x2)))
df_x=as.data.frame(x2)
df_x$year=year(index(x2))
i=1
gini=vector()
for (i in 1:length(years)){
  cur_set=subset(df_x,df_x$year==years[i])
  gini[i] = Gini(cur_set$x2)
}

gini2=(gini-min(gini))/(max(gini)-min(gini))

x3=subset(df_x,df_x$year==years[14])
x4=subset(df_x,df_x$year==years[17])

