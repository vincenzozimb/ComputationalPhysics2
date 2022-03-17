import csv 
import matplotlib.pyplot as plt 

## functions ##

# read csv  
def read_csv(filename):
    with open(filename, 'r') as csvfile: 
        reader = csv.reader(csvfile, delimiter= '\t')
        data = list(reader)
    return data


# plot
def plot_func(data):
    x = [float(i[0]) for i in data]

    j = [float(i[1]) for i in data]
    y = [float(i[2]) for i in data]

    plt.plot(x,j,label = '$j_l$')
    plt.plot(x,y,label = '$n_l$')
    
    plt.legend()
    # plt.grid()
    plt.xlabel('$x$')
    plt.savefig('plotname.png')
    plt.show()

#############################################################


filename='data_gsl.csv'
data = read_csv(filename)
plot_func(data)

filename='data.csv'
data = read_csv(filename)
plot_func(data)