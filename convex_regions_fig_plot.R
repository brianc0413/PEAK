# visualization plots for regions of interest

line_segment_x <- data.frame(x = 1:200/200, y=c(rep(0.5,100), rep(0, 100)))
line_segment_y <- data.frame(x=c(rep(0.5,100), rep(0, 100)), x = 1:200/200)

xy_line <- data.frame(x = 1:200/200, y=1:200/200)

plot(xy_line$x, xy_line$y, type="l", xlab = "Mean of Arm 1", ylab = "Mean of Arm 2") + text(0.75,0.25, "Region 2", col = "blue") + text(0.25,0.75, "Region 1", col = "red")


plot(x=0.5, y=0.5, type="l", xlab = "Mean of Arm 1", ylab = "Mean of Arm 2", xlim = c(0,1),
     ylim = c(0,1)) + text(0.75,0.75, "Region 2", col = "blue") + text(0.25,0.75, "Region 1", col = "red") + 
  text(0.75, 0.25, "Region 3", col = "orange")+ text(0.25, 0.25, "Region 4", col = "purple") + abline(h = 0.5) + abline(v = 0.5)


plot(x=0.5, y=0.5, type="l", xlab = "Mean of Arm 1", ylab = "Mean of Arm 2", xlim = c(0,1),
     ylim = c(0,1)) + text(0.75,0.75, "Region 2", col = "blue") + text(0.25,0.25, "Region 1", col = "red") + 
  segments(x0 = 0, y0=0.5, x1 = 0.5, y1=0.5) +  segments(y0 = 0, x0=0.5, y1 = 0.5, x1=0.5)


plot(x=0.5, y=0.5, type="l", xlab = "Mean of Arm 1", ylab = "Mean of Arm 2", xlim = c(0,1),
     ylim = c(0,1)) + text(0.25,0.75, "Region 2", col = "blue") + text(0.25,0.25, "Region 1", col = "red") + 
  text(0.75, 0.25, "Region 3", col="orange") + 
  segments(x0 = 0, y0=0.5, x1 = 0.5, y1=0.5) +  segments(y0 = 0, x0=0.5, y1 = 0.5, x1=0.5) + 
    segments(x0=0.5, y0 = 0.5, x1 = 1, y1=1)
