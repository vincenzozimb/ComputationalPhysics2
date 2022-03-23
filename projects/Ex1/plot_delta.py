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
    sigma = [float(i[1]) for i in data]
    plt.plot(x,sigma)
    plt.xlabel('E')
    plt.ylabel('$DeltaI(E)')
    plt.savefig('plot.png')
    plt.show()

################

filename='delta.csv'
data = read_csv(filename)
plot_func(data)

