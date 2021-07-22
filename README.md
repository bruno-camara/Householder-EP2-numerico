# Householder-EP2-numerico
Código desenvolvido para o Exercício Programa da matéria Cálculo Numérico e suas Aplicações utilizando o método de tridiagonalização de matrizes de Householder e obtenção de autovalores e autovetores pelo método QR.

### Integrantes:
Bruno Carneiro Camara - 11257230
Leonardo Akira Shimabukuro - 9838053

### Como rodar?
Para rodar o programa é necessário estar no PowerShell do Windows, ou no bash do Linux
    1. Navegue até a pasta em que o arquivo ep2.py está localizado
    2. Digite no terminal: cat input-x | python ep2.py x
        * Sendo x o item a ser verificado - a, b, c ou d
        * *Observação*: o item d corresponde ao bônus e ele utiliza o arquivo input-c como entrada
    
## Entradas:
Como entrada do programa nós temos os arquivos de input disponibilizados:
    * input-a
    * input-b
    * input-c

## Saídas:
As saídas dependerão de cada item
1. item a:
    * Matriz de entrada - A
    * Matriz tridiagonalizada após passar pelo algoritmo de Householder
    * Verificação de Av = _lamba_*v
    * Autovalores obtidos pelo algortimo QR
    * Matriz de autovetores obtida pelo algoritmo QR
    * Verificação da ortogonalidade

2. item b:
    * Matriz de entrada - A
    * Matriz tridiagonalizada após passar pelo algoritmo de Householder
    * Verificação de Av = _lamba_*v
    * Autovalores obtidos pelo algortimo QR
    * Autovalores obtidos pela fórmula do enunciado para comparação com os autovalores obtidos pelo QR
    * Matriz de autovetores obtida pelo algoritmo QR
    * Verificação da ortogonalidade

3. item c:
    * 5 menores frequências e seus respectivos modos de vibração