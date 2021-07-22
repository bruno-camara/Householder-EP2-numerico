import numpy as np
import math
#teste
from numpy import linalg as LA

test = np.array([[2, 1, 0, 0],
                 [1, 2, 1, 0],
                 [0, 1, 2, 1],
                 [0, 0, 1, 2]])

test2 = np.array([[5/math.sqrt(5), 4/math.sqrt(5), 1/math.sqrt(5), 0],
                  [0, 3/math.sqrt(5), 2/math.sqrt(5), 0],
                  [0, 1, 2, 1],
                  [0, 0, 1, 2]])

def getGivensMatrix(matriz, linha):
    """Esta função recebe uma matriz quadrada e o número do elemento da diagonal
    principal (matriz[linha][linha]). A função retorna uma matriz de rotação de
    Givens que zera o valor do elemento abaixo do elemento da diagonal principal
    (matriz[linha][linha+1])."""
    alfa = matriz[linha][linha]
    beta = matriz[linha+1][linha]
    if abs(alfa) > abs(beta):
        tau = -beta/alfa
        c = 1/math.sqrt(1+tau**2)
        s = c*tau
    else:
        tau = -alfa/beta
        s = 1/math.sqrt(1+tau**2)
        c = s*tau
    givensRotationMatrixQ = np.identity(len(matriz))
    givensRotationMatrixQ[linha][linha] = c
    givensRotationMatrixQ[linha][linha+1] = -s
    givensRotationMatrixQ[linha+1][linha] = s
    givensRotationMatrixQ[linha+1][linha+1] = c
    return givensRotationMatrixQ

def getQRDecompositionm(matriz):
    '''Esta função recebe uma matriz tridiagonal simétrica e retorna uma tupla
    contendo as matrizes Q e R da decomposição QR.'''
    tamanho = len(matriz)
    Qt = np.identity(tamanho) # Esta variável contém a matriz Q transposta
    for i in range(tamanho-1):
        Qi = getGivensMatrix(matriz, i)
        matriz = Qi @ matriz
        Qt = Qi @ Qt
    return Qt.transpose(), matriz

def getCSParameters(alpha, beta):
    if (math.fabs(alpha) > math.fabs(beta)):
        tau = -beta/alpha
        c = 1/math.sqrt(1+tau**2)
        s = c*tau
    else:
        tau = -alpha/beta
        s = 1/math.sqrt(1+tau**2)
        c = s*tau
    return c, s

def applyGivens(matriz, tamanho, i, c ,s):
    for j in range(tamanho):
        matriz1j = matriz[i][j]
        matriz2j = matriz[i+1][j]
        matriz[i][j] = c*matriz1j-s*matriz2j
        matriz[i+1][j] = s*matriz1j+c*matriz2j
    return matriz

def applyInverseGivens(matriz, tamanho, j, c, s):
    for i in range(tamanho):
        matrizi1 = matriz[i][j]
        matrizi2 = matriz[i][j+1]
        matriz[i][j] = c*matrizi1-s*matrizi2
        matriz[i][j+1] = s*matrizi1+c*matrizi2
    return matriz

def mulRightQ(matriz, tamanho, c, s):
    for j in range(tamanho-1):
        applyInverseGivens(matriz, tamanho, j, c[j], s[j])
    return matriz

def getQRDecomposition(matriz, tamanho):
    resmat = matriz
    c = []
    s = []
    for i in range(tamanho-1):
        cr, sr = getCSParameters(resmat[i][i], resmat[i+1][i])
        c.append(cr)
        s.append(sr)
        resmat = applyGivens(resmat, tamanho, i, cr, sr)
    return resmat, c, s

def getMik(matriz, n):
    '''Calcula o múltiplo da identidade µk que será utilizado no algoritmo QR
    com deslocamento espectral.'''
    dk = 0.5*(matriz[n-2][n-2] - matriz[n-1][n-1])
    if dk >= 0:
        mik = matriz[n-1][n-1] + dk - np.hypot(dk, matriz[n-1][n-2])
    else:
        mik = matriz[n-1][n-1] + dk + np.hypot(dk, matriz[n-1][n-2])
    return mik

def addDiagonal(matriz, tamanho, MiK):
    for i in range(tamanho):
        matriz[i][i] += MiK

def QRAlgorithm(matriz, V=None):
    '''Executa o algorítmo QR'''
    A = matriz.astype(float)
    np.set_printoptions(suppress=True)
    if V is None:
        V = np.identity(len(A))
    k = 0
    MiK = 0
    for i in range(len(A)-1, 0, -1):
        while True:
            if (k > 0):
                MiK = getMik(A, i+1)
            addDiagonal(A, len(A), -MiK)
            [Q, R] = getQRDecompositionm(A)
            A = (R@Q)
            addDiagonal(A, len(A), MiK)
            V = V@Q
            k += 1
            if abs(A[i][i-1]) < 0.000001:
                A[i][i-1] = 0
                break
    autoval=[]
    for i in range(len(A)):
        autoval.append(A[i][i])
    return autoval, V

#Q, R = getQRDecomposition(test)
#print(Q)
#print(R)
#print(Q @ R)
QRAlgorithm(test)