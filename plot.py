import matplotlib.pyplot as plt
import numpy as np

entrada = input().split() # linha de comentario
dados = []
dados.append([])

for i in range(12):
    entrada = input().split()
    numeroNo = int(entrada[0])
    dados.append(entrada[1:])

#adotando no 14 como origem
nos = [None] * 15

nos[1] = [15, 40]
nos[2] = [5, 40]
nos[3] = [25, 30]
nos[4] = [15, 30]
nos[5] = [5, 30]
nos[6] = [-5, 30]
nos[7] = [25, 20]
nos[8] = [15, 20]
nos[9] = [5, 20]
nos[10] = [-5, 20]
nos[11] = [15,10]
nos[12] = [5,10]
nos[13] = [20,0]
nos[14] = [0,0]
print(dados)
nosorig = nos.copy()

barras = [[1, 2],
          [1, 4],
          [1, 5],
          [2, 5],
          [2, 4],
          [3, 4],
          [3, 7],
          [3, 8],
          [4, 5],
          [4, 8],
          [4, 9],
          [5, 6],
          [5, 9],
          [5, 8],
          [6, 9],
          [6, 10],
          [7, 8],
          [8, 9],
          [8, 11],
          [8, 12],
          [9, 10],
          [9, 11],
          [9, 12],
          [11, 12],
          [11, 13],
          [14, 11],
          [12, 13],
          [14, 12]]

for j in range(int(len(dados[1])/2)):
    plt.cla()
    nos = nosorig
    for k in range(1, 12):
        print(k)
        nos[k][0] += 100*float(dados[k][2*(j-1)])
        nos[k][1] += 100*float(dados[k][2*(j-1)+1])
        print(nos[k][0], nos[k][1])

    for i in range(28):
        x_values = [nos[barras[i][0]][0], nos[barras[i][1]][0]]
        y_values = [nos[barras[i][0]][1], nos[barras[i][1]][1]]
        plt.plot(x_values, y_values)

    plt.savefig(str(j)+'.png')