#! /usr/bin/env python
# -*- coding: utf-8 -*-
from math import log, exp
import random

def potencia_pv():
    r = random.random()
    if 0 <= r < 0.491741:
        G = 0.0
        Ta = 26.35
    elif 0.491741 <= r < 0.618638:
        G = 264.683
        Ta = 30.23
    elif 0.618638 <= r < 0.713229:
        G = 468.267
        Ta = 31.87
    elif 0.713229 <= r < 0.780674:
        G = 971.098
        Ta = 34.06
    elif 0.780674 <= r < 0.854104:
        G = 712.626
        Ta = 33.19
    else:
        G = 66.022
        Ta = 27.39

    #Valores de Referência
    Gref = 1000.0
    Vocref = 32.9
    Iscref = 8.21
    Tcref = 25.0+273.0
    Vmpref = 26.3
    Impref = 7.61

    #Constantes
    Rso = 0.221
    Rsho = 415.405
    uIsc = 3.18*10**-3
    uVoc = -1.23*10**-1
    Eq = 1.124
    k = 1.38064852*10**-23
    q = 1.602176565*10**-19
    Npaineis = 15000.0
    nm = 0.97 #Perdas devido à incompatibilidade entre os módulos
    nd = 0.96 #Perdas devido à sujeira
    ninv = 0.98 #Eficiência do inversor

    #Entrada das Variáveis
    Tc = Ta + (((Tcref-273.0)-20.0)/800.0)*G #Ross, 1980
    Tc = Tc + 273.0;

    #Cálculo dos Parâmetros
    Rsh = Rsho
    Vt = k*Tc/q
    Isc = Iscref*(G/Gref) + uIsc*(Tc - Tcref)
    m1 = 54.0 #Chute inicial
    aux = 0.0 #Variável Auxiliar
    erro = 0.01

    if G != 0.0:
        while(aux == 0.0):
            Voc = Vocref + m1*Vt*log(G/Gref) + uVoc*(Tc - Tcref)
            Vmp = 0.79939*Voc #Método da Tensão Constante
            Imp = 0.92692*Isc #Método da Corrente Constante
            m = (Vmp + Imp*Rso - Voc)/(Vt*(log(abs(Isc - Vmp/Rsh - Imp)) - log(Isc - Voc/Rsh) + Imp/(Isc - Voc/Rsh)))
            if abs(m-m1) <= erro:
                aux = 1.0
            m1 = m
        Io = (Isc - Voc / Rsh) * exp(-Voc / (m * Vt))
        Rs = Rso - (m * Vt / Io) * exp(-Voc / (m * Vt))
        Il = Isc * (1.0 + Rs / Rsh) + Io * (exp(Isc * Rs / (m * Vt)) - 1.0)

        #Cálculo da Corrente de Saída
        I1 = 6.0 #Chute Inicial para a Corrente
        aux = 0.0 #Variável Auxiliar
        erro = 0.01
        V = Vmp #Considerando que o sistema possui MPPT
        while(aux == 0.0):
            I = Il - Io * (exp((V + I1 * Rs)/(m * Vt)) - 1.0) - (V + I1 * Rs)/Rsh
            if abs(I-I1) <= erro:
                aux = 1.0
            I1 = I
        P = nd * nm * ninv * (Npaineis * I * V)
    else:
        P = 0.0

    return P/(10.0**3.0)