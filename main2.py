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

def getWn(A, n):
    '''
    Essa função retorna o wn como vetor coluna
    '''
    an = A[n:len(A), n-1:] #an são as colunas abaixo da diagonal principal dependendo da iteração
    an = an.T
    wnT = np.array([an[0][0] + np.sign(an[0][0])*norma(an[0])]) #Primeiro Valor de wn
    for i in range(1, len(an[0])):
        wnT = np.append(wnT, an[0][i])

    return wnT.T

def getHwnA(an, wn, A):
    an = an.T
    HwnA = np.array([])
    for i in range(an.shape[0]):
        HwnA = np.append(HwnA, an[i] - 2*(escalar(wn, an[i])/escalar(wn, wn))*wn)
    HwnA = HwnA.reshape((an.shape[0], an.shape[1]))
    HwnA = np.insert(HwnA.T, 0, A[0])

    HwnA = HwnA.reshape(A.shape)
    return HwnA

def multiplicaHwnDireita(X, wn, linha, coluna):
    for i in range(linha, len(X)):
        X[i, coluna:] = X[i, coluna:] - 2*(escalar(X[i,coluna:], wn)/escalar(wn, wn))*wn.T
    return X

def multiplicaHwnEsquerda(X, wn, iteracao):
    for i in range(iteracao-1, len(X)):
        X[iteracao:, i] = X[iteracao:, i] - 2*(escalar(wn, X[iteracao:, i])/escalar(wn, wn))*wn
    return X

def rebaterColunaNaLinha(matriz, n):
    matriz[n-1, n:] = matriz[n: , n-1]
    return matriz

def householder(A):
    Ht = np.identity(len(A))
    HwnAHwn = A
    for i in range(1, len(A)-1):
        print(i)
        wn = getWn(HwnAHwn, i)
        HwnA = multiplicaHwnEsquerda(HwnAHwn, wn, i)
        HwnAHwn = multiplicaHwnDireita(HwnA, wn, i, i) #Tanto a linha quanto a coluna tem que ser o valor da iteração

        HwnAHwn = rebaterColunaNaLinha(HwnAHwn, i)
        print(HwnAHwn)
        Ht = multiplicaHwnDireita(Ht, wn, 0, i)
        print(HwnAHwn)
    return (HwnAHwn, Ht)

def main():
    entrada = np.array([[2, -1, 1, 3], [-1, 1, 4, 2], [1, 4, 2, -1], [3, 2, -1, 1]])
    A = entrada.astype(float)
    print(A)
    HwnAHwn, Ht = householder(A)
    print("HwnAHwn:")
    print(HwnAHwn)
    print("\nHt:")
    print(Ht)

    #Atenção os valores passados nessas funções estão como apresentados no exemplo

if __name__ == '__main__':
    main()

#Aula no moodle pra tirar duvida

#Para concatenar matrizes e usando a função pra expandir a linha de Hw1a1
    # linha0 = linhaMatrizHwnA(wn, A, 0)
    # linha1 = linhaMatrizHwnA(wn, A, 1)
    # print(linha0)
    # print(linha1)
    # print(np.concatenate((linha0, linha1), axis=1))