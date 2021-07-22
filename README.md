# Householder-EP2-numerico
Código desenvolvido para o Exercício Programa da matéria Cálculo Numérico e suas Aplicações utilizando o método de tridiagonalização de matrizes de Householder e obtenção de autovalores e autovetores pelo método QR.

## Integrantes:
Bruno Carneiro Camara - 11257230<br>
Leonardo Akira Shimabukuro - 9838053

## Como rodar?
O código foi feito em python 3.7.2 <br>
Para rodar o programa é necessário estar no PowerShell do Windows, ou no bash do Linux
1. Navegue até a pasta em que o arquivo ep2.py está localizado
2. O comando python ep2.py -h exibe informacoes sobre o uso
3. Para executar o item a digite no terminal: `cat input-a | python ep2.py a`
4. Para executar o item b digite no terminal: `cat input-b | python ep2.py b`
5. Para executar o item c digite no terminal: `cat input-c | python ep2.py c`
6. Para executar o item c digite no terminal: `cat input-c | python ep2.py d | python plot.py`
    * **Observação:** para executar a tarefa d, é possível especificar os seguintes parâmetros.
        - A opção -f <opcao> especifica qual das menores frequências de vibração será utilizada para a geração das imagens (1=menor freq, 5=maior) (padrão=1)
        - A opção -n <int> especifica quantos quadros serão gerados (padrão=4)
        - A opção -t <float> especifica o intervalo de tempo entre os quadros (padrao=200)
    * **Observação 2:** Para gerar as imagens da tarefa d, é necessário ter a biblioteca matplotlib instalada.
        - É possível instalá-la utilizando o comando python -m pip install -U matplotlib

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

4. item bônus:
    * O vídeo saida.avi contém uma animação gerada utilizando o comando
        - `cat input-c | python ep2.py d -f 1 -n 120 -t 0.5 | py plot.py`
        - As imagens geradas por este comando se encontram na pasta imagens

