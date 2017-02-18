"""Cigre-MonteCarlo.

    Este software foi desenvolvido para o trabalho final da
    disciplina de Analise de Disturbios de Tensao de Curta Duracao
    ministrada no Programa de Pos Graduacao da Universidade
    Federal do Ceara.

    A base do software e a biblioteca pandapower que disponibiliza
    funcionalidades como representacao da rede e de seus elementos
    principais e calculo de fluxo de carga.
    Os parametros de niveis de tensao mediante curto-circuito foram
    desenvolvimentos proprios da equipe:

    Os envolvidos com o desenvolvimento deste software sao:

    - Fellipe Souto
    - Lucas Melo
    - Nestor Fontinele
    - Francisco Neto

    Esse software e livre e pode ser utilizado assim como esta
    ou com alteracoes de qualquer forma.
    Caso voce venha a utilizar esse software, de alguma forma,
    pedimos que considere nosso esforco de desenvolvimento e
    cite os desenvolvedors.
"""

import pandapower as pp
import pandapower.networks as pn
import networkx as nx
import numpy as np
import pandas as pd
from tqdm import tqdm

from utilidades.Incertezas import calcula_parametros
from utilidades.Potencia_Eolica import potencia_eolica
from utilidades.Potencia_PV import potencia_pv

from utilidades.protecao import Rele

# configuracoes do pandas

pd.set_option('display.max_columns', 15)
pd.set_option('display.max_rows', 20)
pd.set_option('precision', 3)


def simula_curto(net, n_linha, p_line, tipo, z_cc, barras, barras_gd, i_cc_gd):
    """Realiza a simulacao da falta na rede."""
    # -------------------------------- #
    # calcula a impedancia equivalente #
    # -------------------------------- #
    linha = net.line.query(expr=str(n_linha))
    
    # carrega a topologia da rede em forma de grafos

    net_graph = pp.topology.create_nxgraph(net)

    # caminho da barra em falta ate a barra da subestacao
    c1 = nx.shortest_path(net_graph, linha.from_bus, 0) 
    c2 = nx.shortest_path(net_graph, linha.to_bus, 0)

    if len(c1) < len(c2):
        c = c1
    else:
        c = c2

    b = c[0]  # barra imediatamente a montante da linha sob falta
    c.pop()
    c.reverse()
    
    d = list(c)
    d.pop(0)

    # caso a linha selecionada esteja conectada a barra da subestacao
    # nao sera necesario calcular as impedancias equivalestes das
    # outras barras
    if d != []:
        # percorre as linhas entre a barra em falta ate a barra da subestacao
        for i, j in zip(c, d):
            l = net.line.query('to_bus==' + str(i) + ' and ' + 'from_bus==' + str(j))
            if len(l) == 0:
                l = net.line.query('to_bus==' + str(j) + ' and ' + 'from_bus==' + str(i))

            z_linha = float(l.length_km) * complex(l.r_ohm_per_km, l.x_ohm_per_km)
            z0_linha = float(l.length_km) * complex(l.r0_ohm_per_km, l.x0_ohm_per_km)

            if net.bus.zeq_ohms[j] == 0.0j:
                net.bus.zeq_ohms[j] = net.bus.zeq_ohms[i] + z_linha

            if net.bus.z0_eq_ohms[j] == 0.0j:
                net.bus.z0_eq_ohms[j] = net.bus.z0_eq_ohms[i] + z0_linha

        # impedancia equivalente ate a barra da linha onde ocorreu a falta
        zeq_ohms = net.bus.zeq_ohms[j]
        z0_eq_ohms = net.bus.z0_eq_ohms[j]
    else:
        zeq_ohms = net.bus.zeq_ohms[b]
        z0_eq_ohms = net.bus.z0_eq_ohms[b]

    # impedancias da barra da linha ate o ponto onde ocorreu a falta
    zeq_ohms_ = p_line * float(linha.length_km) * complex(linha.r_ohm_per_km,
                                                          linha.x_ohm_per_km)
    z0_eq_ohms_ = p_line * float(linha.length_km) * complex(linha.r0_ohm_per_km,
                                                            linha.x0_ohm_per_km)
    # impedancia de seq. posit e neg totais da falta
    zeq_ohms_total = zeq_ohms + zeq_ohms_
    z0_eq_ohms_total = z0_eq_ohms + z0_eq_ohms_

    # ----------------------------- #
    # calculo das tensoes pre-falta #
    # ----------------------------- #
    bus = net.bus.query(str(b))
    vn_bus = bus.vn_kv * 1e3 / np.sqrt(3)

    bus_ = net.res_bus.query(str(b))
    vr_bus_pu_m = bus_.vm_pu
    vr_bus_dg_a = bus_.va_degree
    # ---------------------------- #
    # def. das impedancias de seq. #
    # ---------------------------- #

    z1 = zeq_ohms_total
    z2 = z1
    z0 = z0_eq_ohms_total

    # ------------------------------------------------------- #
    # calculo da contribuicao da corrente de falta do sistema #
    # de acordo com o tipo de falta selecionada.              #
    # ------------------------------------------------------- #
    if tipo == '3f':
        i_cc = abs((vn_bus * vr_bus_pu_m) / (z1 + z_cc))
    elif tipo == 'ft':
        i_cc = 3.0 * abs((vn_bus * vr_bus_pu_m) / (z0 + z1 + z2 + 3.0 * z_cc))
    elif tipo == '2f':
        i_cc = np.sqrt(3.0) * abs((vn_bus * vr_bus_pu_m) / (z1 + z2 + z_cc))
    elif tipo == '2ft':
        i1 = (vn_bus * vr_bus_pu_m) / (z1 + (z2 * (z0 + 3.0 * z_cc)) / (z2 + z0 + 3.0 * z_cc))
        i2 = - i1 * (z0 + 3.0 * z_cc) / (z0 + 3.0 * z_cc + z2)
        i0 = - i1 * (z2) / (z0 + 3.0 * z_cc + z2)

        a = -0.5 + (np.sqrt(3) / 2.0) * 1.0j
        a2 = a**2

        A = np.array([[1.0, 1.0, 1.0], [1.0, a2, a], [1.0, a, a2]])
        Is = np.array([i0, i1, i2])

        i_cc = np.max(np.abs(np.dot(A, Is)))

    c.reverse()
    d.reverse()

    # --------------------------------------------------------- #
    # calcula as correntes de falta circulantes na rede         #
    # --------------------------------------------------------- #
    
    net.bus['v_cc_volts'] = np.zeros(len(net.bus))

    if barras_gd != []:
        # ------------------------------------------------- #
        # calcula o caminho entre a gd e a linha sob falta  #
        # ------------------------------------------------- #
        caminhos_gd = list()
        for i in barras_gd:
            # caminho da barra da gd ate a barra da linha em falta
            # mais proxima da gd
            c1 = nx.shortest_path(net_graph, linha.from_bus, i) 
            c2 = nx.shortest_path(net_graph, linha.to_bus, i)

            if len(c1) < len(c2):
                c_gd = c1
            else:
                c_gd = c2
            
            caminhos_gd.append(c_gd)

        # --------------------------------------------------------- #
        # Calcula as correntes de falta que passam pelos trafos     #
        # --------------------------------------------------------- #
        
        # os vizinhos do no zero sao os alimentadores 
        aliment = nx.neighbors(net_graph, 0)
        aliment = {i:[] for i in aliment}

        # percorre os alimentadores para associar as suas respectivas GDs
        for i in aliment.keys():
            for j in barras_gd:
                c1 = nx.shortest_path(net_graph, j, 0)
                if i in c1:
                    aliment[i].append(j)

        for i in net.trafo.index:
            # associa a corrente de falta do sistema
            # ao transformador que alimenta o alimentador
            # que contem a linha sob falta
            if net.trafo.lv_bus[i] in c:
                net.trafo.i_cc[i] += i_cc
            # percorre as GDs para associar suas respectivas correntes 
            # de falta
            for k, j in enumerate(barras_gd):
                c1 = nx.shortest_path(net_graph, j, b)
                if net.trafo.lv_bus[i] in c1:
                    net.trafo.i_cc[i] += i_cc_gd[k]

        # --------------------------------------------------------- #
        # Insere as correntes de falta causadas pelas GD            #
        # --------------------------------------------------------- #

        # insere a corrente de falta devido a geracao distribuida
        # nas linhas entre a GD e a falta
        for k, cam in enumerate(caminhos_gd):
            for i, j in zip(cam, cam[1:]):
                l = net.line.query('to_bus==' + str(i) + ' and ' + 'from_bus==' + str(j))
                if len(l) == 0:
                    l = net.line.query('to_bus==' + str(j) + ' and ' + 'from_bus==' + str(i))
                    if len(l) == 0:
                        continue
                net.line.i_cc[int(l.index[0])] += i_cc_gd[k]
    
    else: # if barras_gd != []:
        for i in net.trafo.index:
            # associa a corrente de falta do sistema
            # ao transformador que alimenta o alimentador
            # que contem a linha sob falta
            if net.trafo.lv_bus[i] in c:
                net.trafo.i_cc[i] += i_cc

    # ---------------------------------------------------------------- #
    # Insere a corrente de falta causada pela geracao principal        #                  #
    # ---------------------------------------------------------------- #
    for i, j in zip(d, c[1:]):

        l = net.line.query('to_bus==' + str(i) + ' and ' + 'from_bus==' + str(j))
        if len(l) == 0:
            l = net.line.query('to_bus==' + str(j) + ' and ' + 'from_bus==' + str(i))
        
        net.line.i_cc[int(l.index[0])] += i_cc

    # --------------------------------------------------------------------- #
    # Calcula os afundamentos de tensao causados pelas correntes de falta   #                  #
    # --------------------------------------------------------------------- #

    # ---------------------------------------------------------- #
    # Calculo da tensao nas duas extremidades
    # da linha em falta ou seja:
    # vb1 = (pos) * (len_linha) * (z_linha) * icc_1
    # vb2 = (1.0 - pos) * (len_linha) * (z_linha) * icc_2
    # 
    # Onde vb1 e a tensao na barra mais proxima da SE e vb2
    # e a tensap na barra mais distante da SE, logo
    # vb2 != 0.0 somente se houver GD a jusante, de tal modo que
    # icc_1 != icc_2
    # ---------------------------------------------------------- #

    # calcula as quedas de tensao nas barras de extremidade
    # da linha em falta
    b1 = net.bus.query(str(linha.from_bus))
    b2 = net.bus.query(str(linha.to_bus))

    lines_1 = net.line.query('from_bus==' + str(linha.from_bus) + ' or ' + 'to_bus==' + str(linha.from_bus))

    icc_1 = 0.0
    if d != []:
        for l in zip(lines_1.from_bus, lines_1.to_bus, lines_1.i_cc):
            if l[0] != linha.from_bus and l[1] != linha.to_bus:
                icc_1 += l[2]
    else:
        icc_1 = i_cc
    # print icc_1

    lines_2 = net.line.query('from_bus==' + str(linha.to_bus) + ' or ' + 'to_bus==' + str(linha.to_bus))
    icc_2 = 0.0
    for l in zip(lines_2.from_bus, lines_2.to_bus, lines_2.i_cc):
        if l[0] != linha.from_bus and l[1] != linha.to_bus:
            icc_2 += l[2]
    # print icc_2

    z = p_line * float(linha.length_km) * complex(linha.r_ohm_per_km, linha.x_ohm_per_km)
    z_ = (1.0 - p_line) * float(linha.length_km) * complex(linha.r_ohm_per_km, linha.x_ohm_per_km)
    if icc_1 != 0.0:
        net.bus.v_cc_volts[b1.name] = abs((icc_1 + icc_2) * z_cc + icc_1 * z)
    if icc_2 != 0.0:
        net.bus.v_cc_volts[b2.name] = abs((icc_1 + icc_2) * z_cc + icc_2 * z_)

    # calcula as quedas de tensao nas barras observadas

    for i in barras:
        # caminho da barra em falta ate a barra observada
        c1 = nx.shortest_path(net_graph, linha.from_bus, i) 
        c2 = nx.shortest_path(net_graph, linha.to_bus, i)

        if len(c1) > len(c2):
            c1 = c2

        for i, j in zip(c1, c1[1:]):
            l = net.line.query('to_bus==' + str(i) + ' and ' + 'from_bus==' + str(j))
            if len(l) == 0:
                l = net.line.query('to_bus==' + str(j) + ' and ' + 'from_bus==' + str(i))
                if len(l) == 0:
                    # ---------------------------------
                    # Falta configurar o caso de curto
                    # no alimentador 2 com GD no ali-
                    # mentador 1, ou seja, a corrente
                    # de falta passa pelos trafos 
                    # ---------------------------------
                    if barras_gd != [] and j != 0:
                        if aliment[j] != []:
                            t = net.trafo.query('hv_bus==' + str(i) + ' and ' + 'lv_bus==' + str(j))
                            if len(t) == 0:
                                t = net.trafo.query('hv_bus==' + str(j) + ' and ' + 'lv_bus==' + str(i))
                            net.bus.v_cc_volts[j] = abs(net.bus.v_cc_volts[i] + t.z_ohms.values[0] * t.i_cc.values[0])
                    continue

            if l.i_cc.values[0] != 0.0:
                zeq_ohms_ = float(l.length_km) * complex(l.r_ohm_per_km,
                                                         l.x_ohm_per_km)
                net.bus.v_cc_volts[j] = abs(net.bus.v_cc_volts[i] + zeq_ohms_ * net.line.i_cc[int(l.index[0])])

    dados = list()
    for i in barras:
        dados.append([i, net.bus.v_cc_volts[i] / (net.bus.vn_kv[i] * 1e3 / np.sqrt(3))])
        # print 'barra: ', i, ' ', net.bus.v_cc_volts[i] / (net.bus.vn_kv[i] * 1e3 / np.sqrt(3))

    return dados

def inicializa_rede(barras_gd, p_pv_kw, p_eolica_kw):
    """Metodo de inicializacao da rede."""
    # seleciona a rede cigre media tensao
    net = pn.create_cigre_network_mv(with_der=False)

    if barras_gd != []:
        load = net.load.query('bus==' + str(barras_gd[0]))
        net.load.p_kw[load.index[0]] -= p_pv_kw
        net.load.q_kvar[load.index[0]] -= p_pv_kw/0.9 * np.sin(np.arccos(0.9))

        load = net.load.query('bus==' + str(barras_gd[1]))

        net.load.p_kw[load.index[0]] -= p_eolica_kw
        net.load.q_kvar[load.index[0]] -= p_eolica_kw / 0.9 * np.sin(np.arccos(0.9))

    # roda o fluxo de carga
    pp.runpp(net)

    # ------------------------------- #
    # insere a impedancia de sq. zero #
    # ------------------------------- #

    r0_ohm_per_km = 3.0 * 0.501
    x0_ohm_per_km = 3.0 * 0.716

    net.line['r0_ohm_per_km'] = r0_ohm_per_km * np.ones(len(net.line))
    net.line['x0_ohm_per_km'] = x0_ohm_per_km * np.ones(len(net.line))

    # ------------------------------- #
    # insere a impedancia equivalente #
    # ------------------------------- #

    fat = 0.25
    net.bus['zeq_ohms'] = np.zeros(len(net.bus))
    net.bus.zeq_ohms[1] = fat * (0.488 + 5.832j)
    net.bus.zeq_ohms[12] = fat * (0.488 + 5.832j)

    net.bus['z0_eq_ohms'] = np.zeros(len(net.bus))
    net.bus.z0_eq_ohms[1] = fat * (0.444 + 5.34j)
    net.bus.z0_eq_ohms[12] = fat * (0.444 + 5.34j)

    # ------------------------------------- #
    # insere a impedancia do trafo          #
    # ------------------------------------- #
    net.trafo['z_ohms'] = np.zeros(len(net.trafo))
    net.trafo.z_ohms[0] = 0.5 + 3.2j
    net.trafo.z_ohms[1] = 0.5 + 3.2j


    # ------------------------------------- #
    # insere a corrente defalta dos trafos  #
    # ------------------------------------- #
    net.trafo['i_cc'] = np.zeros(len(net.trafo))

    # ------------------------------------- #
    # insere a corrente de falta nas linhas #
    # ------------------------------------- #
    net.line['i_cc'] = np.zeros(len(net.line))

    return net

# ------------------------------------- #
# executa a simulacao de curto-circuito #
# ------------------------------------- #
if __name__ == '__main__':
    
    # simulacao de Monte Carlo
    n = 100 # numero de anos de simulacao
    
    barras_gd = [5, 10]
    # barras_gd = []

    rele = Rele(0.1, 1.3 * 150.0, 'ni')

    parametros = list()
    afundamentos = list()
    # loop sobre o numero de anos de simulacao
    for i in tqdm(range(1, n+1), ascii=True, desc='Simu. de Monte Carlo', ncols=100):
        # loop sobre cada uma das linhas da rede
        for j in tqdm(range(12), ascii=True, desc='Falta em linha', ncols=100):
            # loop sobre cada uma das faltas na linha em um ano

            net = inicializa_rede(barras_gd, potencia_pv(), potencia_eolica())
            for k in calcula_parametros(net, j):
                p_pv = potencia_pv()
                p_eolica = potencia_eolica()
                net = inicializa_rede(barras_gd, p_pv, p_eolica)
                l = [i, j]
                l.extend(k)
                parametros.append(l)

                # contribuicao de correntes de falta das GD
                alpha = 3.0
                i_cc_pv = alpha * p_pv / (np.sqrt(3) * 20)
                i_cc_eolica = alpha * p_eolica / (np.sqrt(3) * 20)
                # simula a falta para os parametros fornecidos em l e de geracao distribuida
                dados = simula_curto(net, l[1], l[2], l[3], l[4], [1, 3, 5, 10, 11], barras_gd, [i_cc_pv, i_cc_eolica])
                for m in dados:
                    m.append(i)
                    m.append(j)
                    m.append(len(parametros)-1)
                    m.append(rele.t_op_atuacao(net.trafo.i_cc.max()))
                    afundamentos.append(m)

    parametros = pd.DataFrame(parametros, columns=['ano', 'linha', 'posicao', 'tipo', 'resistencia'])
    afundamentos = pd.DataFrame(afundamentos, columns=['barra_obs', 'valor_afund', 'ano', 'linha_em_falta', 'n_aft', 'tempo'])

    writer = pd.ExcelWriter('dados_com_gd.xlsx', engine='xlsxwriter')
    afundamentos.to_excel(writer, sheet_name='afundamentos')
    parametros.to_excel(writer, sheet_name='parametros')
    writer.save()

    # Escrita dos resultados em uma planilha do excel


    # teste da funcao simula_curto(rede,
    #                              linha_em_falta,
    #                              posicao_da_falta,
    #                              tipo_de_falta,
    #                              resistencia_de_falta,
    #                              barras_monitoradas,
    #                              barras_com_gd)
    #

    # net = inicializa_rede([5, 6], potencia_pv(), potencia_eolica())
    # dados = simula_curto(net, 1, 0.8, '3f', 2.0, [1, 2, 3, 4, 5, 8, 12], [5, 6])
    # dados = simula_curto(net, 1, 0.5, '3f', 5.0, [1, 2, 3, 4, 5], [5, 4])
