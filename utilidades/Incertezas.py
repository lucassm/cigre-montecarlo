# -*- coding: utf-8 -*-
"""."""
import random
from math import log, exp


def calcula_parametros(net, numero_da_linha):
    """."""
    # Tipo de linha:
    # 1 - Tronco primário.
    # 2 - Derivação lateral.
    # 3 - Linha secundária.
    tipo_de_linha = 2

    # Corrente que flui no trecho
    I = net.res_line.loading_percent[numero_da_linha] * 1e-2

    # Comprimento da linha em milhas
    L = net.line.length_km[numero_da_linha] / 1.609344

    # Sorteio do Numero de Faltas

    # Calculo da taxa de falha
    if tipo_de_linha is 1:
        A = 0.01976
        B = 3.4295969
        C = -0.009756098
    elif tipo_de_linha is 2:
        A = 0.07759
        B = 2.1522789
        C = -0.067586207
    else:
        A = 0.01402
        B = 3.7632316
        C = -0.004018433
    a = 10.0
    c = 0.5085
    r = 1 / (1 + exp(-a * (I - c)))
    yp = (A * exp(B * r) + C) * L
    y = 5.7143 * yp # Faltas momentaneas correspondem a 82,5% de todas as faltas

    # Distribuicao de Poisson
    aux = 0
    n = 1
    a = 1.0
    while aux == 0:
        Un = random.random()
        a = a * Un
        if a >= exp(-y):
            n = n + 1
        else:
            X = n - 1
            aux = 1

    parametros = list()

    for i in range(X):

        # Sorteio do tipo de falta

        # Percentuais de tipo de falta encontrados em Bordalo,
        # Rodrigues e Silva (2006)
        r = 1e2 * random.random()
        if 0.0 <= r < 3.0:
            tipo_de_falta = '3f'
        elif 3.0 <= r < 13.0:
            tipo_de_falta = '2f'
        elif 13.0 <= r < 19.0:
            tipo_de_falta = '2ft'
        else:
            tipo_de_falta = 'ft'

        # Sorteio da Posicao da Falta
        
        # Distribuicao Uniforme
        posicao = random.random()

        # Sorteio da Resistencia de Falta

        # Distribuição de Weibull
        U = random.random()
        Rmin = 1.0
        Rmax = 3.0
        intervalo_confianca = 99.0
        Finf = ((100.0 - intervalo_confianca) / 2.0) / 100.0
        Fsup = 1.0 - ((100.0 - intervalo_confianca) / 2.0) / 100.0
        Finf_barrado = 1.0 - Finf
        Fsup_barrado = 1.0 - Fsup
        b = log(log(1.0 / Finf_barrado) / log(1.0 / Fsup_barrado)) / log(Rmin / Rmax)
        a = Rmax * (-log(Fsup_barrado))**(-1.0 / b)
        R_falta = a * (-log(1.0 - U))**(1.0 / b)

        parametros.append((posicao, tipo_de_falta, R_falta))

    return parametros


if __name__ == '__main__':
    calc_falta()