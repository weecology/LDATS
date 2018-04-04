
model_list = r_LDA


plot.LDA_list <- function(model_list){

model <- model_list[[8]]

model@gamma




}

cols<-c(rgb(0,0,1), rgb(0,1,0))
cols<-(rgb(runif(10,0,1), runif(10,0,1), runif(10,0,1)) )

plot.LDA_VEM <- function(model, cols){

gamma <- model@gamma
ntopics <- ncol(gamma)

slotNames(model)



par(fig = c(0, 1, 0, 0.8))
par(mar = c(4, 4, 1, 1))
plot(gamma[ , 1], type = "n", bty = "L", xlab = "", ylab = "", las = 1,
     ylim = c(0, 1))
mtext(side = 1, "Sample", line = 2.75, cex = 1.5)
mtext(side = 2, "Proportion", line = 2.75, cex = 1.5)

for (i in 1:ntopics){
  points(gamma[ , i], col = cols[i], type = "l", lwd = 1)
}


}





