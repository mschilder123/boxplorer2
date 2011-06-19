mencoder -nosound -ovc lavc -lavcopts vcodec=mjpeg -o boxplorer-30fps.avi -mf type=tga:fps=30 mf://frame-*.tga 
