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
    '''
    Essa função retorna uma matriz com os a's minúsculos em colunas para cada iteração de wn
    Retorna uma partição da matriz A, que será usada na iteração n para calcular o HwnA
    '''
    an = np.array([])
    an = A[n:len(A), n-1:]
    return an

def getWnT(an):
    '''
    Essa função retorna o wn como vetor coluna
    '''
    an = an.T
    wnT = np.array([an[0][0] + np.sign(an[0][0])*norma(an[0])]) #Primeiro Valor de wn
    for i in range(1, len(an[0])):
        wnT = np.append(wnT, an[0][i])

    return wnT.T

#A função abaixo não está sendo usada
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

def getAlfa(A):
    alfa = A[1: , 0]
    return alfa

def getAHwnAlfa(alfa, wn): #nome da função ta errado, não tá?
    Hwnalfa = alfa - 2*(escalar(wn, alfa)/escalar(wn, wn))*wn
    return Hwnalfa

def getHwnAnHwn(HwnA, wn, n):
    '''
    Essa função seria a multiplicação de Hwn pela direita.
    n é o número da iteração
    '''
    HwnAnwn = np.array([])
    for i in range(1, len(HwnA)):
        HwnAnwn = np.append(HwnAnwn, HwnA[i,n:] - 2*(escalar(HwnA[i,n:], wn)/escalar(wn, wn))*wn.T)
    HwnAnwn = HwnAnwn.reshape([len(HwnA)-n, len(HwnA)-n])
    return HwnAnwn

def juntarHwnAHwn(HwnAlfa, HwnAnHwn, n):
    HwnAHwn = np.identity(len(HwnAnHwn) + n)
    HwnAHwn[n-1, n:] = HwnAlfa
    HwnAHwn[n:, n-1] = HwnAlfa.T
    HwnAHwn[n: , n: ] = HwnAnHwn
    return HwnAHwn



def main():
    #Primeiro Exemplo feito "na mão"
    A = np.array([[2, -1, 1, 3], [-1, 1, 4, 2], [1, 4, 2, -1], [3, 2, -1, 1]])

    #Primeira iteração
    an = getAn(A, 1)
    wn = getWnT(an)
    HwnA = getHwnA(an, wn, A)
    alfa = getAlfa(A)
    Hw1alfa = getAHwnAlfa(alfa, wn)
    print("Hw1alfa é: \n", Hw1alfa)
    Hw1A1Hw1 = getHwnAnHwn(HwnA, wn, 1) #1 sendo o número da iteração
    print('Hw1A1Hw1 é: \n', Hw1A1Hw1)
    Hw1AHw1 = juntarHwnAHwn(Hw1alfa, Hw1A1Hw1, 1)
    print('Hw1AHw1 é: \n', Hw1AHw1)

    Ht = np.identity(len(A)) #Iniciando a matriz de vetores
    Ht = getHwnAnHwn(Ht, wn, 1) #Multiplicando Hw1 pela esquerda - só isso que essa função faz
    print('Ht é: \n', Ht) #Falta deixar ela 4x4. Aqui ela está 3x3
    
    #Fim da primeira iteração

    #Atenção os valores passados nessas funções estão como apresentados no exemplo

if __name__ == '__main__':
    main()



#Para concatenar matrizes e usando a função pra expandir a linha de Hw1a1
    # linha0 = linhaMatrizHwnA(wn, A, 0)
    # linha1 = linhaMatrizHwnA(wn, A, 1)
    # print(linha0)
    # print(linha1)
    # print(np.concatenate((linha0, linha1), axis=1))