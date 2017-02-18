#! /usr/bin/env python
# -*- coding: utf-8 -*-
from math import log, exp
import random


def potencia_eolica():
    # Distribuição de Weibull
    U = random.random()
    Vmin = 0.14
    Vmax = 14.27
    intervalo_confianca = 99.0
    Finf = ((100.0 - intervalo_confianca) / 2.0) / 100.0
    Fsup = 1.0 - ((100.0 - intervalo_confianca) / 2.0) / 100.0
    Finf_barrado = 1.0 - Finf
    Fsup_barrado = 1.0 - Fsup
    b = log(log(1.0 / Finf_barrado) / log(1.0 / Fsup_barrado)) / log(Vmin / Vmax)
    a = Vmax * (-log(Fsup_barrado))**(-1.0 / b)
    Vwind = (a * (-log(1.0 - U))**(1.0 / b)) #Velocidade do Vento

    #Cálculo da Potência de Saída
    Pnom = 3.45 * 10.0**6
    Vnom = 11.5
    d = 1.225 #Pressão de 1 atm e temperatura de 15 ºC
    A = 12469.0
    Cp = 2 * Pnom/(d * A * Vnom**3)


    if Vwind < 3.0 or Vwind > 22.5:
        Pwind = 0.0
    elif Vwind >= Vnom:
        Pwind = Pnom
    else:
        Pwind = Cp * (0.5) * d * A * Vwind**3

    return Pwind/10.0**3
