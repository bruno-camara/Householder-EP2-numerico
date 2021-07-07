#Algoritmo de Householder
import numpy as np
import math

def norma(a):
    somaquadrada = 0
    for i in range (len(a)):
        somaquadrada += math.pow(math.fabs(a[i]),2)
    return math.sqrt(somaquadrada)

def escalar(wnT, a0):
    #Calcula o produto escalar entre dois numpy arrays
    return sum(wnT*a0)

def getAn(A, n):
    '''Essa função retorna uma matriz com os a's minúsculos em colunas para cada iteração de wn'''
    tamanho = len(A)
    an = np.array([])
    an = A[n:len(A), n-1:]
    return an

def getWnT(an):
    an = an.T
    wnT = np.array([an[0][0] + np.sign(an[0][0])*norma(an[0])]) #Primeiro Valor de wn
    for i in range(1, len(an[0])):
        wnT = np.append(wnT, an[0][i])

    return wnT.T

def linhaMatrizHwnA(wn, A, n): #n é o numero da linha a ser preenchida
    tamanhoA = len(A)
    tamanhoWn = len(wn)
    linha = np.zeros((len(A), 1)) #Valor Inicial
    for i in range(len(linha)):
        if (i >= tamanhoA - tamanhoWn):
            linha[i] = wn[i-(tamanhoA - tamanhoWn)]
        else:
            linha[i] = A[i][n]
    return linha

def getHwnA(an, wn, A):
    an = an.T
    HwnA = np.array([])
    for i in range(an.shape[0]):
        HwnA = np.append(HwnA, an[i] - 2*(escalar(wn, an[i])/escalar(wn, wn))*wn)
    HwnA = HwnA.reshape((an.shape[0], an.shape[1]))
    HwnA = np.insert(HwnA.T, 0, A[0])

    HwnA = HwnA.reshape(A.shape)
    return HwnA


def main():
    #Primeiro Exemplo feito "na mão"
    A = np.array([[2, -1, 1, 3], [-1, 1, 4, 2], [1, 4, 2, -1], [3, 2, -1, 1]])

    an = getAn(A, 1)
    wn = getWnT(an)
    HwnA = getHwnA(an, wn, A)
    print("A matriz Hw1A é \n", HwnA)
    
    #Atenção os valores passados nessas funções estão como apresentados no exemplo

if __name__ == '__main__':
    main()



#Para concatenar matrizes e usando a função pra expandir a linha de Hw1a1
    # linha0 = linhaMatrizHwnA(wn, A, 0)
    # linha1 = linhaMatrizHwnA(wn, A, 1)
    # print(linha0)
    # print(linha1)
    # print(np.concatenate((linha0, linha1), axis=1))