import matplotlib.pyplot as plt

names = []
marks = []

f = open('gnuplot.txt','r')
for row in f:
    row = row[:-1]
    row = row.split(' ')
    names.append(row[0])
    marks.append(row[1])
plt.plot(marks, '.')
plt.xlabel('G', fontsize = 12)
plt.ylabel('x', fontsize = 12)

plt.title('Gaussian')
plt.legend()
plt.show()
