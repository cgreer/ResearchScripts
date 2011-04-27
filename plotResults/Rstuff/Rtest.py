from rpy2.robjects import r
from rpy2.robjects.packages import importr
gDevice = importr('grDevices')

gDevice.png(file="file.png", width=512, height=512)
# plotting code here

x = [1,2,3]
y = [1,2,3]

r.plot(x,y, xlab = "X")

gDevice.dev_off()
