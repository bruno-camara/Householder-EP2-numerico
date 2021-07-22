# Exercicio Programa 2 - Calculo Numerico

# Bruno Carneiro Camara - 11257230
# Leonardo Akira Shimabukuro - 9838053

#Algoritmo de Householder
import numpy as np
import math
import argparse
import sys

#teste
from ep1.ep1 import QRAlgorithm
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
    entrada = input().split()
    dimensao = int(entrada[0])
    A = np.empty((dimensao, dimensao))
    for i in range(dimensao):
        entrada = input().split(" ")
        for j in range(dimensao):
            A[i][j] = float(entrada[j])
    matrizoriginal = np.copy(A) #Para a verificação ao final
    print("Matriz de entrada: ")
    print(A)
    

    HwnAHwn, Ht = householder(A)
    print("\nMatriz tridiagonalizada: ")
    print(HwnAHwn)

    w, v = QRAlgorithm(HwnAHwn, Ht) #EP1
    #Verificar A*v = lambda*v
    erro = 0
    for i in range(len(w)):
        print("\nVerificando A*v = w*v para "+str(i+1)+"º autovetor")
        print("V = ", v[:,i])
        print("A*v = ", np.matmul(matrizoriginal, v[:, i]))
        print("w*v = ", w[i]*v[:, i])
        erro = max(max(np.matmul(matrizoriginal, v[:, i])-w[i]*v[:,i]), erro)
    print("O maior erro entre A*v e w*v foi "+str(erro))
    
    print("\nAutovalores:")
    print(w)
    print("\nMatriz de autovetores:")
    print(v)

    #Verificar se a matriz formada pelos autovetores é ortogonal
    print("\nVerificando a ortogonalidade:")
    print(np.matmul(v, v.T))
    print("Como a multiplicação da matriz de autovetores pela sua transposta é a própria identidade a matriz de autovetores é ortogonal")
    
def testeB():
    entrada = input().split()
    dimensao = int(entrada[0])
    A = np.empty((dimensao, dimensao))
    for i in range(dimensao):
        entrada = input().split()
        for j in range(dimensao):
            A[i][j] = float(entrada[j])
    matrizoriginal = np.copy(A)
    #A = np.zeros((20,20))
    #for i in range(20):
    #    for j in range(20):
    #        A[i][j] = 20-max(i,j)
    print("Matriz de entrada: ")
    print(A)

    HwnAHwn, Ht = householder(A)
    print("\nMatriz tridiagonalizada: ")
    print(HwnAHwn)

    w, v = QRAlgorithm(HwnAHwn, Ht) #EP1
    #Verificar A*v = lambda*v
    erro = 0
    for i in range(len(w)):
        print("\nVerificando A*v = w*v para "+str(i+1)+"º autovetor")
        print("V = ", v[:,i])
        print("A*v = ", np.matmul(matrizoriginal, v[:, i]))
        print("w*v = ", w[i]*v[:, i])
        erro = max(max(np.matmul(matrizoriginal, v[:, i])-w[i]*v[:,i]), erro)
    print("O maior erro entre A*v e w*v foi "+str(erro))

    print("\nAutovalores obtidos pelo algoritmo")
    print(w)

    w_correcao = np.array([])
    for i in range(1, 21):
        w_correcao = np.append(w_correcao, (1/2)*pow((1 - math.cos((2*i - 1)*math.pi/(2*20 + 1))), -1))
    print("\nAutovalores obtidos pela fórmula:")
    print(w_correcao)

    print("\nMatriz de autovetores:")
    print(v)

    #Verificar se a matriz formada pelos autovetores é ortogonal
    print("\nVerificando a ortogonalidade:")
    print(np.matmul(v, v.T))
    print("Como a multiplicação da matriz de autovetores pela sua transposta é a própria identidade a matriz de autovetores é ortogonal")

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

def tarefaC(isBonus):
    entrada = input().split()
    totalNos = int(entrada[0])
    nosNaoFixos = int(entrada[1])
    numBarras = int(entrada[2])

    entrada = input().split()
    densidade = float(entrada[0])
    secaoA = float(entrada[1])
    elasticidadeE = float(entrada[2])*10**9
    
    K = np.zeros((totalNos*2+1, totalNos*2+1))
    M = np.zeros((totalNos*2+1))
    for i in range(numBarras):
        entrada = input().split()
        no1 = int(entrada[0])
        no2 = int(entrada[1])
        angulo = float(entrada[2])
        comprimento = float(entrada[3])
        Kij = getKij(secaoA, elasticidadeE, comprimento, angulo)
        K = writeToK(K, Kij, no1, no2)
        M[2*no1] += 0.5*densidade*secaoA*comprimento
        M[2*no1-1] += 0.5*densidade*secaoA*comprimento
        M[2*no2] += 0.5*densidade*secaoA*comprimento
        M[2*no2-1] += 0.5*densidade*secaoA*comprimento
    K = K[1:nosNaoFixos*2+1, 1:nosNaoFixos*2+1]
    M = M[1:nosNaoFixos*2+1]

    Ktil = np.zeros((nosNaoFixos*2, nosNaoFixos*2))
    for i in range(len(K)):
        Ktil[i, :] = K[i, :] * M[i]**(-1/2)
    for i in range(len(K)):
        Ktil[:, i] = Ktil[:, i] * M[i]**(-1/2)

    np.set_printoptions(suppress=True)

    matrizoriginal = np.copy(Ktil)
    Ktil, Ht = householder(Ktil)

    autoval, autovet = QRAlgorithm(Ktil, Ht) # EP1
    # autoval, autovet = LA.eig(matrizoriginal)

    #for i in range(len(autoval)):
    #    print(i+1, autoval[i])
    #    print(i+1, )

    autoval = np.power(autoval, (1/2))
    indexMenores = [0] * 5
    valorMenores = [0] * 5
    zipped_list = zip(autoval, autovet)
    sorted_pairs = sorted(zipped_list)
    tuples = zip(*sorted_pairs)
    autoval, autovet = [ list(tuple) for tuple in tuples]
    for i in range(len(autovet)):
        autovet[i] = M[i]**(-1/2) * autovet[i]
    if not isBonus:
        print("As 5 menores frequências são: ")
        for i in range(5):
            print("\n["+str(i+1)+"] Freq: ", autoval[i])
            print("["+str(i+1)+"] Modo: ", autovet[i])
    else:
        print("# i", end=" ")
        for i in range(args.numquadros):
            print("Xi("+str(i*args.tempo)+") Yi("+str(i*args.tempo)+")", end=" ")
        print()
        for i in range(nosNaoFixos):
            print(i+1, end=" ")
            for j in range(args.numquadros):
                print(autovet[args.frequencia-1][2*i]*math.cos(autoval[args.frequencia-1]*j*args.tempo), end=" ")
                print(autovet[args.frequencia-1][2*i+1]*math.cos(autoval[args.frequencia-1]*j*args.tempo), end=" ")
            print()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="EP2 de MAP3121")
    parser.add_argument('tarefa', choices=['a','b','c','d'], help='qual tarefa a ser executada (a, b, c ou d)')
    parser.add_argument('-f', '--frequencia', type=int, choices=[1,2,3,4,5], help="qual das 5 menores frequências de vibração será usada para gerar as imagens", default=1)
    parser.add_argument('-n', '--numquadros', type=int, help="quantas figuras devem ser geradas", default=4)
    parser.add_argument('-t', '--tempo', type=float, help="intervalo de tempo entre figuras", default=200)
    # parser.add_argument()
    args = parser.parse_args()
    if args.tarefa == "a":
        testeA()
    elif args.tarefa == "b":
        testeB()
    elif args.tarefa == "c":
        tarefaC(False)
    elif args.tarefa == "d":
        tarefaC(True)
