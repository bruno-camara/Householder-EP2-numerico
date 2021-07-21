#Algoritmo de Householder
import numpy as np
import math

#teste
from numpy import linalg as LA

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
        wn = getWn(HwnAHwn, i)
        HwnA = multiplicaHwnEsquerda(HwnAHwn, wn, i)
        HwnAHwn = multiplicaHwnDireita(HwnA, wn, i, i) #Tanto a linha quanto a coluna tem que ser o valor da iteração
        HwnAHwn = rebaterColunaNaLinha(HwnAHwn, i)
        Ht = multiplicaHwnDireita(Ht, wn, 0, i)
    return (HwnAHwn, Ht)

def testeA():
    A = np.array([[2, 4, 1, 1], [4, 2, 1, 1], [1, 1, 1, 2], [1, 1, 2, 1]])
    A = A.astype(float)
    print(A)
    HwnAHwn, Ht = householder(A)
    #Ep1
    w, v = LA.eig(A)
    print(w)
    print(v)
    #Verificar A*v = lambda*v
    #Verificar se a matriz formada pelos autovetores é ortogonal
        #Basta multiplicar a matriz de autovetores pela sua transposta e verificar se resulta na identidade
    
def testeB():
    A = np.zeros((20,20))
    for i in range(20):
        for j in range(20):
            A[i][j] = 20-max(i,j)
    print(A)
    w, v = LA.eig(A)
    print(w)
    #print(v)
    w_correcao = np.array([])
    for i in range(1, 21):
        w_correcao = np.append(w_correcao, (1/2)*pow((1 - math.cos((2*i - 1)*math.pi/(2*20 + 1))), -1))
    print(w_correcao)

def getKij(A, E, L, theta):
    c = math.cos(math.radians(theta))
    s = math.sin(math.radians(theta))
    matriz = np.array([[c**2, c*s, -1*c**2, -1*c*s], [c*s, s**2, -1*c*s, -1*s**2], [-1*c**2, -1*c*s, c**2, c*s], [-1*c*s, -1*s**2, c*s, s**2]])
    Kij = (A*E/L) * matriz
    return Kij

def writeToK(K, Kij, i, j):
    K[2*i-1][2*i-1] += Kij[0][0]
    K[2*i-1][2*i] += Kij[0][1]
    K[2*i-1][2*j-1] += Kij[0][2]
    K[2*i-1][2*j] += Kij[0][3]
    
    K[2*i][2*i-1] += Kij[1][0]
    K[2*i][2*i] += Kij[1][1]
    K[2*i][2*j-1] += Kij[1][2]
    K[2*i][2*j] += Kij[1][3]
    
    K[2*j-1][2*i-1] += Kij[2][0]
    K[2*j-1][2*i] += Kij[2][1]
    K[2*j-1][2*j-1] += Kij[2][2]
    K[2*j-1][2*j] += Kij[2][3]

    K[2*j][2*i-1] += Kij[3][0]
    K[2*j][2*i] += Kij[3][1]
    K[2*j][2*j-1] += Kij[3][2]
    K[2*j][2*j] += Kij[3][3]
    return K

def tarefaC():
    entrada = input().split(" ")
    totalNos = int(entrada[0])
    nosNaoFixos = int(entrada[1])
    numBarras = int(entrada[2])

    print(totalNos)
    entrada = input().split(" ")
    densidade = float(entrada[0])
    secaoA = float(entrada[1])
    elasticidadeE = float(entrada[2])
    
    K = np.zeros((totalNos*2+1, totalNos*2+1))
    M = np.zeros((totalNos*2+1, totalNos*2+1))
    for i in range(numBarras):
        entrada = input().split(" ")
        no1 = int(entrada[0])
        no2 = int(entrada[1])
        angulo = float(entrada[2])
        comprimento = float(entrada[3])
        Kij = getKij(secaoA, elasticidadeE, comprimento, angulo)
        K = writeToK(K, Kij, no1, no2)
        M[2*no1][2*no1] += 0.5*densidade*secaoA*comprimento
        M[2*no1-1][2*no1-1] += 0.5*densidade*secaoA*comprimento
        M[2*no2][2*no2] += 0.5*densidade*secaoA*comprimento
        M[2*no2-1][2*no2-1] += 0.5*densidade*secaoA*comprimento
    K = K[1:nosNaoFixos*2, 1:nosNaoFixos*2]
    M = M[1:nosNaoFixos*2, 1:nosNaoFixos*2]

    
def main():
    #Primeiro Exemplo feito "na mão"
    entrada = np.array([[2, -1, 1, 3], [-1, 1, 4, 2], [1, 4, 2, -1], [3, 2, -1, 1]])
    A = entrada.astype(float)
    HwnAHwn, Ht = householder(A)
    #testeA()
    #testeB()
    tarefaC()

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