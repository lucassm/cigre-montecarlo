import numpy as np
import matplotlib.pyplot as plt

# Funcao de Protecao de Sobrecorrente Temporizada

# para essa funcao entra-se com o tipo de curva, o DIAL de tempo
# e a corrente medida no secundario do TC e retorna-se o tempo 
# de atuacao

# t_op = DIAL * beta / (I**alpha - 1.0)

class Rele(object):
    def __init__(self, dial, i_pickup, curva):
        
        self.dial = dial
        self.i_pickup = i_pickup

        if curva == 'ni':
            self.alpha = 0.02
            self.beta = 0.14
        elif curva == 'mi':
            self.alpha = 1.0
            self.beta = 13.5
        elif curva == 'ei':
            self.alpha = 2.0
            self.beta = 80.0

    def t_op_atuacao(self, i):
        
        if i > 5.0 * self.i_pickup:
            return 0.1
        else:
            M = i / self.i_pickup
            return self.dial * self.beta / (M**self.alpha - 1.0)

    def plot_coordenograma(self):
        i = np.arange(self.i_pickup+0.1, 4.0 * self.i_pickup, 0.01)
        M = i / self.i_pickup
        t = self.dial * self.beta / (M**self.alpha - 1.0)

        plt.loglog(i, t, linewidth=3)
        plt.grid(True, which='both', linestyle='-')
        plt.axis([0.0, 10.0 * self.i_pickup, 0.0, 1e3])
        plt.show()
