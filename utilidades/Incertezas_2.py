#! /usr/bin/env python
# -*- coding: utf-8 -*-
import random
from math import log, exp

linha = input("Entre com o tipo de linha: 1 - Tronco primário; 2 - Derivação lateral; 3 - Linha secundária\n")
I = float(input("Entre com a corrente que flui no trecho: "))
ampacidade = float(input("Entre com a ampacidade do condutor: "))
I = I/ampacidade

#Percentuais encontrados em Bordalo, Rodrigues e Silva (2006)
r = 100*random.random()
if 0 <= r < 3:
    tipo = 3
elif 3 <= r < 19:
    tipo = 2
else:
    tipo = 1

#Posição da Falta
#Distribuição de Weibull
U = random.random()
Rmin = 1.0
Rmax = 3.0
intervalo_confianca = 99.0
Finf = ((100-intervalo_confianca)/2.0)/100.0
Fsup = 1.0 - ((100-intervalo_confianca)/2.0)/100.0
Finf_barrado = 1.0 - Finf
Fsup_barrado = 1.0 - Fsup
b = log(log(1.0/Finf_barrado)/log(1.0/Fsup_barrado))/log(Rmin/Rmax)
a = Rmax*(-log(Fsup_barrado))**(-1.0/b)
R_falta = a*(-log(1.0-U))**(1/b)

#Posição da Falta
#Distribuição Uniforme
posicao = random.random()

#Número de Faltas
#Cálculo da taxa de falha
if linha is 1:
   A = 0.01976
   B = 3.4295969
   C = -0.009756098
elif linha is 2:
   A = 0.07759
   B = 2.1522789
   C = -0.067586207
else:
   A = 0.01402;
   B = 3.7632316
   C = -0.004018433
a = 10
c = 0.5085
r = 1/(1+exp(-a*(I-c)))
yp = A*exp(B*r) + C
y = 5.7143*yp #Faltas momentâneas correspondem a 82,5% de todas as faltas

#Distribuição de Poisson
aux = 0
n = 1
a = 1
while aux == 0:
    Un = random.random()
    a = a*Un
    if a >= exp(-y):
        n = n + 1
    else:
        X = n - 1
        aux = 1

#Resultados
print("Quantidade de Faltas: " + str(X))
print("Posição da Falta: " + str(posicao))
print("Tipo de Falta: " + str(tipo))
print("Resistência de Falta: " + str(R_falta))

