if 'patterns' not in globals():
   patterns=[[1,1], [2,4], [3,9]]

if 'x' not in globals():
   var('x')

if 'f' not in globals():
   f(x) = x

k(x) = e^(-x^2/2)/sqrt(2*pi)

def gp(n):
   del patterns[:]
   for i in range(n):
      ran=[gauss(0,1) for i in range(2)]
      patterns.append([ran[0], N(k(ran[0])) + 0.01*ran[1]])
   return point2d(patterns, size=12)


def fit(lam):
    l=[]
    for i in range(len(patterns)):
       xi=patterns[i][0]
       l.append(k((x-xi)/lam))
    d(x)=sum(l)
    del l[:]
    for i in range(len(patterns)):
      xi=patterns[i][0]
      l.append(k((x-xi)/lam)/d(x))
    r=0
    for i in range(len(patterns)):
      yi=patterns[i][1]
      r += l[i]*yi
    return r

def lfunctions(lam):
    l=[]
    for i in range(len(patterns)):
       xi=patterns[i][0]
       l.append(k((x-xi)/lam))
    d(x)=sum(l)
    del l[:]
    for i in range(len(patterns)):
      xi=patterns[i][0]
      l.append(k((x-xi)/lam)/d(x))
    r=0
    return l

colours = ["red", "blue", "orange", "green", "purple"]

def fit_and_plot(log_lam):
   g = point2d(patterns, size=15)
   lam = 10^log_lam
   f=fit(lam)
   var('t')
   g = g + plot(f(x=t), [t,-3,3])
   return g

def show_lfunctions():
   g = point2d(patterns, size=30)
   j = 0
   for i in range(1,3):
      l=lfunctions(1/2^i)
      for f in l:
         g = g + plot(f, [x, -5, 5], color=colours[j % 5])
         j += 1
   return g
shl=show_lfunctions

def in_sample_error(lam):
   sum = 0
   for p in patterns:
      xx = p[0]
      y = p[1]
      y_hat = fit(lam)(x=xx)
      sum += (y_hat - y)^2
   sum = sum/len(patterns)
   return sum

def gcve(lam):
   return in_sample_error(lam)/(1.0 - degrees(lam)/len(patterns))^2

def degrees(lam):
   sum = 0
   for i in range(len(patterns)):
      sum += N(lfunctions(lam)[i](patterns[i][0]))
   return sum

def show_validation_error(min, max):
   '''
   plots the validation error as a function of log(lambda), where
   lambda is the bandwidth. the plot is displayed for log(lamda)
   between min and max.
   '''
   num_points = 100
   alist = []
#   blist = []
   for i in range(num_points):
      log_lam = min + i*(max-min)/(num_points-1)
      lam = 10^log_lam
      y =  generalized_cross_validation_error(patterns, lam)
#      y1 = in_sample_error(lam)
#      y1 = gcve(lam)
      alist.append([log_lam, y])
#      blist.append([log_lam, y1+1])
   print "Error at max = %f" % (y)
   return list_plot(alist)
shve=show_validation_error

def show_degrees(min, max):
   num_points = 100
   alist = [];
   for i in range(num_points):
      log_lam = min + i*(max-min)/(num_points-1)
      lam = 10^log_lam
      y = degrees(lam)
      alist.append([log_lam, y])
   return list_plot(alist)
shd=show_degrees

   

