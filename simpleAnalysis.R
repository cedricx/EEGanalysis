require('reshape2')
require('ggplot2')
require('rasterVis')
require('cowplot')

data <- read.csv('~/Google Drive/NeuroFlow/Sample CSV.csv')
timepoints <- c(1,which(as.numeric(data$text)>1))
timepoint_meaning <- c("Calibration","Baseline 1","Meditation","A_V","Stress","Baseline 2")
num_tps <- length(timepoints)/2
num_freq <- dim(data)[2]-3

#set up processed data table

uniData <- matrix(NA,nrow = num_tps, ncol = num_freq)
colnames(uniData) <- colnames(data)[(1:num_freq)+3]
rownames(uniData) <- timepoint_meaning
biData <- matrix(NA,nrow = num_tps, ncol = num_freq/2)
colnames(biData) <- c("delta","theta","alpha","beta","gamma")
rownames(biData) <- timepoint_meaning

finalTable <- NULL
finalTable$mean <- uniData
finalTable$mad <- uniData
finalTable$std <- uniData
finalTable$cor <- biData
finalTable$diff <- biData
finalTable$bimean <- biData


#set up time points of interest
for (i in 1:num_tps) {
  tp_1 <- timepoints[i*2-1]+1
  tp_2 <- timepoints[i*2]-1
  #for each freq band on EACH side
  for (j in (1:num_freq)+3) {
    finalTable$mean[i,j-3] <- mean(data[tp_1:tp_2,j]) #calculate ave. freq. bands
    finalTable$mad[i,j-3] <- mad(data[tp_1:tp_2,j]) #calculate variance based on median absolute deviation
    finalTable$std[i,j-3] <- sd(data[tp_1:tp_2,j]) #calculate variance based on standard deviation
  }
  #for each freq band on BOTH sides
  for (k in seq(4,13,by=2)) {
    finalTable$cor[i,k/2-1] <- cor(data[tp_1:tp_2,k],data[tp_1:tp_2,k+1]) #calculate pearson correlations
    finalTable$diff[i,k/2-1] <- mean(data[tp_1:tp_2,k] - data[tp_1:tp_2,k+1]) #asymmetry Right minus Left
    finalTable$bimean[i,k/2-1] <- mean(c(data[tp_1:tp_2,k], data[tp_1:tp_2,k+1])) #total mean
  }
}

#plot the results
resultsPlot <- NULL
plotName <- c("Average","Variance(mad)","Variance(std)","Correlation","Asymmetry (R-L)","Total Average")
LR_order <- c(seq(1,10,by=2),seq(2,10,by=2))
for (i in 1:3) {
resultsPlot[[i]] <- levelplot(finalTable[[i]][,LR_order], par.settings = YlOrRdTheme(),
          xlab="Conditions",ylab = "Frequency Bands",
          scales=list(x=list(rot=45, tck = 0),y=list(tck = 0)),
          main = list(label = plotName[i]))
}
for (i in 4:5) {
  datalim <- max(abs(finalTable[[i]]))
  resultsPlot[[i]] <- levelplot(finalTable[[i]], 
                                par.settings = BuRdTheme(),
                                at = seq(-datalim,datalim,length.out = 20),
                                xlab="Conditions",ylab = "Frequency Bands",
                                scales=list(x=list(rot=45, tck = 0),y=list(tck = 0)),
                                main = list(label = plotName[i]))
}
for (i in 6) {
  datalim <- max(abs(finalTable[[i]]))
  resultsPlot[[i]] <- levelplot(finalTable[[i]], 
                                par.settings = YlOrRdTheme(),
                                xlab="Conditions",ylab = "Frequency Bands",
                                scales=list(x=list(rot=45, tck = 0),y=list(tck = 0)),
                                main = list(label = plotName[i]))
}

tpPlot <-function(calcTable, title) {
  foo <- melt(calcTable)
  colnames(foo) <- c("Conditions","Frequency","relPower")
  foo$Conditions<-factor(foo$Conditions,levels =foo$Conditions )
  plot<-ggplot(foo,aes(x = Conditions, y = relPower, color = Frequency, group = Frequency)) + 
    geom_point() + geom_line() +
    ggtitle(title)
}

for (i in seq_along(finalTable)) {
  resultsTPplots[[i]] <- tpPlot(finalTable[[i]],plotName[i])
}


#graph analysis
make_eeg_graph<-function(data,tp_1,tp_2){
eeg_graph <- matrix(NA,num_freq,num_freq)
  for (i in (1:num_freq)+3) {
    for (j in (1:num_freq)+3) {
      eeg_cor <- cor(data[tp_1:tp_2,i],data[tp_1:tp_2,j])
      eeg_graph[i-3,j-3] <- eeg_cor
    }
  }
return(eeg_graph)
}

eeg_graphs <- NULL
for (i in 1:num_tps) {
  tp_1 <- timepoints[i*2-1]+1
  tp_2 <- timepoints[i*2]-1
  eeg_graphs[[i]]<-make_eeg_graph(data,tp_1,tp_2)
}

graphLabels <- c("d_R","t_R",'a_R','b_R','g_R',"d_L","t_L",'a_L','b_L','g_L')
plot_eeg_graphs <- function(eeg_graph,title){
  levelplot(eeg_graph[LR_order,LR_order],par.settings = BuRdTheme(),
            scales=list(x=list(at= 1:10, labels=graphLabels,rot=90, tck = 0 ),
                        y=list(at= 1:10, labels=graphLabels, tck = 0)),
            xlab="freq bands",ylab= "freq bands", main = list(label=title))
}

eeg_graphs_plots<-lapply(seq_along(eeg_graphs),function(i) plot_eeg_graphs(eeg_graphs[[i]], timepoint_meaning[i]))


#save plots
#resultsplots
save_plots <- function(plot_to_save,gen_name, w,h, type){
  if (type == "ggplot") {
    plotname <- paste('~/Google Drive/NeuroFlow/figures/',gen_name,'_',plot_to_save$labels$title,".pdf",sep="")
  }  else if (type == "levelplot") {
    plotname <- paste('~/Google Drive/NeuroFlow/figures/',gen_name,'_',plot_to_save$main$label,".pdf",sep="")
  }
  pdf(file = plotname, width = w, height = h)
  print(plot_to_save)
  dev.off()
}
lapply(seq_along(resultsPlot), function(plot) save_plots(resultsPlot[[plot]],'level', 5,5, "levelplot")  )
lapply(seq_along(resultsTPplots), function(plot) save_plots(resultsTPplots[[plot]],'tp', 7,4, "ggplot")  )
lapply(seq_along(eeg_graphs_plots), function(plot) save_plots(eeg_graphs_plots[[plot]],'net', 5,5, "levelplot")  )
