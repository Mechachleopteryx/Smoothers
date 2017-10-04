attach smoothers.spyx

if 'patterns' not in globals():
#   patterns=[[-3,-1], [-2, 0], [-1, 1.5], [0,2], [1,1.5],[2,0],[3,-1]]
   patterns=[[-2,1],[-1,-0.5],[0,1.5],[1,-0.5],[2,1]]
if 'raw_patterns' not in globals():
   raw_patterns=[]

if 'x' not in globals():
   var('x')

if 'means' not in globals():
    means = range(2)

if 'sds' not in globals():
    sds= range(2)

colours = ["red", "blue", "orange", "green", "purple"]


def string_to_floats(str):
    '''
    converts a string of space separated floats into a list of floats
    '''
    ret = []
    lst = str.split()
    for i in lst:
        ret.append(float(i))
    return ret

def ssv_to_list_of_lists(filename):
    '''
    A function to convert a space separated data file <filename> into
    a list (row) of lists (cols)
    '''
    infile = open(filename)
    ret = []
    for line in infile:
#        ret.append(tuple(string_to_floats(line)))
        ret.append((string_to_floats(line)))
    return ret

def save_pattern_to_file(filename):
   outfile = open(filename, 'w')
   for p in patterns:
      outfile.write("%.13f %.13f\n" % tuple(p))
   outfile.close()

def normalize_patterns():
   '''
   replace patterns with a normalized version (original is saved as
   <raw_patterns>)
   '''
#   compute the mean and standard deviation of the inputs and outputs
#   storing them in <means> and <sds> lists
   del means[:]
   del sds[:]
   for i in range(2):
      m = 0
      for p in patterns:
         m = m + p[i]
      m = m/len(patterns)
      means.append(m)
      sd = 0
      for p in patterns:
         sd = sd + (p[i]-m)^2
      sd = sqrt(sd/(len(patterns) - 1))
      sds.append(sd)
   del raw_patterns[:]
   for p in patterns:
      raw_patterns.append(p)
   for i in range(len(patterns)):
      for j in range(2):
         patterns[i][j]=n((patterns[i][j] - means[j])/sds[j])

def unnormalize():
   '''
   replace patterns with its original unnormalized version 
   '''
   del raw_patterns[:]
   for p in raw_patterns:
      patterns.append(raw_patterns)

def f(t):
   return t*sin(1/(0.1 + t^2))
#   return tanh(t)

def generate_patterns(n):
   '''
   Erases the global list patterns, writing n new elements randomly
   generated according to some formula (look and see what it is!)
   '''
   del patterns[:]
   for i in range(n):
      ran=[gauss(0,1) for j in range(2)]
      xx = -3.0+(6.0*i)/n
      patterns.append([xx, f(xx)+0.05*ran[1]])
#      patterns.append([25 + 4*ran[0], 0.5*(1000 + 100*N(5*tanh(ran[0])) + 70*ran[1])])
   return point2d(patterns, size=12)
gp = generate_patterns

def show_smoother(log_h, order = 1, xmin = -5, xmax = 5, show_data = 1):
    ''' 
    Plots the linear smoother fitted to the list of patterns called
    <patterns>, using a bandwidth of 10^log_h, and using locally
    linear smoothing for order = 1 and locally constant smoothing for
    order = 0
    ''' 
    g = point2d(patterns, size=30, color="black")
    h = 10^log_h
    def u(x): return prediction(patterns, h, order, x)
    if show_data == 1:
       return g + plot(u, (xmin, xmax), color="blue")
    else:
       return plot(u, (xmin, xmax), color="blue")
shs = show_smoother

def show_order_one_lfunctions(h, xmin = -5, xmax = 5):
   ''' 
   plots all the basis functions for the locally linear smoother
   fitting the global variable <patterns>
   '''
   g = list_plot(patterns, size=30, color='black')
   for j in range(len(patterns)):
      def u(x): return order_one_lfunction(j, patterns, h, x)
      g = g + plot(u, [xmin, xmax], color=colours[j % 5])
   return g

def show_smoothers(min, max, order = 1, xmin = -5, xmax = 5):
   ''' 
   '''
   num_smoothers = 3
   g = list_plot(patterns, size=50, color='black')
   for i in range(num_smoothers):
      log_h = min + i*(max-min)/(num_smoothers-1)
      h = 10^log_h
      def u(x): return prediction(patterns, h, order, x)
      g = g + plot(u, [xmin, xmax], color=colours[i % 5], legend_label = " $h = %.3f$" % h)
   return g
shss=show_smoothers

def show_validation_error(min, max, order = 1):
   '''
   plots the generalized cross-validation error for a linear smoother
   fitting <patterns>, as a function of log(h), where h is the
   bandwidth. the plot is displayed for log(h) between min and max.
   Order = 0 is locally constant smoothing and 1 for locally linear smoothing
   '''
   num_points = 20
   alist = []
   for i in range(num_points):
      log_h = min + i*(max-min)/(num_points-1)
      h = 10^log_h
      y =  gcv_error(patterns, h, order)
#      alist.append([h, y])
      alist.append([log_h, y])
   print "Error at max = %f" % y
   return list_plot(alist)#, scale = 'semilogx')
shve=show_validation_error

def validation_error(log_h, order = 1):
   '''
   returns the generalized cross-validation error associated with
   <patterns>; <log_h> is log_10 of the bandwidth h
   '''
   h=10^log_h
   return gcv_error(patterns, h, order)

def show_degrees(min, max, order = 1):
   ''' 
   plots the effective number of degrees of freedom as a function of
   log_10 of the bandwidth, for <patterns>
   '''
   num_points = 100
   alist = [];
   for i in range(num_points):
      log_h = min + i*(max-min)/(num_points-1)
      h = 10^log_h
      y = effective_degrees(patterns, h, order)
      alist.append([h, y])
#      alist.append([log_h, y])
   print "degrees at min = %f, at max = %f" % (alist[0][1], alist[-1][1])
   return list_plot(alist, scale='semilogx')
shd=show_degrees

   

