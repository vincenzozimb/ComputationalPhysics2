import csv 
import matplotlib.pyplot as plt 


# functions  
def read_csv(filename):
    with open(filename, 'r') as csvfile: 
        reader = csv.reader(csvfile, delimiter= '\t')
        data = list(reader)
    return data

# plot: 
def plot_func(data):
    x = [float(i[0]) for i in data]
    pot = [float(i[1]) for i in data]
    real = [float(i[2]) for i in data]
    imag = [float(i[3]) for i in data]
    plt.plot(x,pot,label = 'potential')
    plt.plot(x,real,label = 'real part')
    plt.plot(x,imag,label = 'imag part')
    plt.legend()
    plt.xlabel('x')
    plt.savefig('plot.png')
    plt.show()

################

filename='solution.csv'
data = read_csv(filename)
plot_func(data)



