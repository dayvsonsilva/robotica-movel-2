# Dayvson Silva
# 06/10/2021
from math import pi, cos, sin, sqrt, atan, tanh
import matplotlib.pyplot as plt
import numpy as np

########################################################################################################################
# Dados simulação
tsim = 50  # Tempo da simulação em segundos
Tc = 0.1  # Tempo de amostragem
ciclos = int(tsim / Tc)  # Ciclos nescessarios para completar simulação em x segundos
print('ciclos nescessarios para a simulação: ', ciclos)

# Definição de variaveis
x = np.zeros(ciclos, dtype=float)
y = np.zeros(ciclos, dtype=float)
z = np.zeros(ciclos, dtype=float)
psi = np.zeros(ciclos, dtype=float)

xf = np.zeros(ciclos, dtype=float)
yf = np.zeros(ciclos, dtype=float)
zf = np.zeros(ciclos, dtype=float)
psif = np.zeros(ciclos, dtype=float)

x1 = np.zeros(ciclos, dtype=float)
y1 = np.zeros(ciclos, dtype=float)
z1 = np.zeros(ciclos, dtype=float)
psi1 = np.zeros(ciclos, dtype=float)

x2 = np.zeros(ciclos, dtype=float)
y2 = np.zeros(ciclos, dtype=float)
z2 = np.zeros(ciclos, dtype=float)
psi2 = np.zeros(ciclos, dtype=float)

rof = np.zeros(ciclos, dtype=float)
alfaf = np.zeros(ciclos, dtype=float)
betaf = np.zeros(ciclos, dtype=float)

xtil = np.zeros(ciclos, dtype=float)
xtilponto = np.zeros(ciclos, dtype=float)

ytil = np.zeros(ciclos, dtype=float)
ytilponto = np.zeros(ciclos, dtype=float)

ztil = np.zeros(ciclos, dtype=float)
ztilponto = np.zeros(ciclos, dtype=float)

rotil = np.zeros(ciclos, dtype=float)
alfatil = np.zeros(ciclos, dtype=float)
betatil = np.zeros(ciclos, dtype=float)

Uvx = np.zeros(ciclos, dtype=float)
Uvy = np.zeros(ciclos, dtype=float)
Uvz = np.zeros(ciclos, dtype=float)
omega = np.zeros(ciclos, dtype=float)

xponto1 = np.zeros(ciclos, dtype=float)
yponto1 = np.zeros(ciclos, dtype=float)
zponto1 = np.zeros(ciclos, dtype=float)
psiponto1 = np.zeros(ciclos, dtype=float)

xponto2 = np.zeros(ciclos, dtype=float)
yponto2 = np.zeros(ciclos, dtype=float)
zponto2 = np.zeros(ciclos, dtype=float)
psiponto2 = np.zeros(ciclos, dtype=float)

xdes = np.zeros(ciclos, dtype=float)
ydes = np.zeros(ciclos, dtype=float)
zdes = np.zeros(ciclos, dtype=float)

rofdes = np.zeros(ciclos, dtype=float)
alfafdes = np.zeros(ciclos, dtype=float)
betafdes = np.zeros(ciclos, dtype=float)

rodes = np.zeros(ciclos, dtype=float)
alfades = np.zeros(ciclos, dtype=float)
betades = np.zeros(ciclos, dtype=float)

xdesponto = np.zeros(ciclos, dtype=float)
ydesponto = np.zeros(ciclos, dtype=float)
zdesponto = np.zeros(ciclos, dtype=float)

rodesponto = np.zeros(ciclos, dtype=float)
alfadesponto = np.zeros(ciclos, dtype=float)
betadesponto = np.zeros(ciclos, dtype=float)

psid = np.zeros(ciclos, dtype=float)
psitil = np.zeros(ciclos, dtype=float)
Ucref = np.zeros(ciclos, dtype=float)
Vcref = np.zeros(ciclos, dtype=float)
Wcref = np.zeros(ciclos, dtype=float)
erro = np.zeros(ciclos, dtype=float)
der_psi = np.zeros(ciclos, dtype=float)

vxref1 = np.zeros(ciclos, dtype=float)
vyref1 = np.zeros(ciclos, dtype=float)
vzref1 = np.zeros(ciclos, dtype=float)
vpsi1 = np.zeros(ciclos, dtype=float)

vxref2 = np.zeros(ciclos, dtype=float)
vyref2 = np.zeros(ciclos, dtype=float)
vzref2 = np.zeros(ciclos, dtype=float)
vpsi2 = np.zeros(ciclos, dtype=float)

Ucrefd = np.zeros(ciclos, dtype=float)
Wcrefd = np.zeros(ciclos, dtype=float)
Vcrefd = np.zeros(ciclos, dtype=float)

qdes = np.zeros((ciclos, 6))
qpontodes = np.zeros((ciclos, 6))
qtil = np.zeros((ciclos, 6))
q = np.zeros((ciclos, 6))
qpontoref = np.zeros((ciclos, 6))

########################################################################################################################
# Camada de planejamento
# Inputs
# Outputs qdes, qpontodes
# Configuração da simulação
# 1 - circulo
# 2 - oito deitado
simulacao = 1
if simulacao == 1:
    kx = ky = kz = 7  # constante proporcional para x e y
    kpsi = 7  # constante proporcional para psi
    lx = ly = lz = 0.1
    lpsi = 1
    for j in range(0, ciclos, 1):
        # xdes[j] = 3 * cos(Tc * j)
        # yd[j] = 3 * sin(Tc * j)
        xdes[j] = 3 + 5 * cos(3 * j * Tc)
        ydes[j] = 5 + 5 * sin(3 * j * Tc)
        zdes[j] = 20
        rodes[j] = 5
        alfades[j] = 0
        betades[j] = 0
        qdes[j] = np.array([xdes[j], ydes[j], zdes[j], rodes[j], alfades[j], betades[j]])

    # Derivada das trajetórias de interesse
    for j in range(0, ciclos, 1):
        # xdponto[j] = -3 * sin(Tc * j)
        # ydponto[j] = 3 * cos(Tc * j)
        xdesponto[j] = np.diff([xdes[j - 1], xdes[j]], n=1) / Tc
        ydesponto[j] = np.diff([ydes[j - 1], ydes[j]], n=1) / Tc
        zdesponto[j] = 0
        rodesponto[j] = 0
        alfadesponto[j] = 0
        betadesponto[j] = 0

        qpontodes[j] = np.array(
            [xdesponto[j], ydesponto[j], zdesponto[j], rodesponto[j], alfadesponto[j], betadesponto[j]])

if simulacao == 2:
    kx = ky = kz = 7  # constante proporcional para x e y
    kpsi = 7  # constante proporcional para psi
    lx = ly = lz = 0.1
    lpsi = 1
    for j in range(0, ciclos, 1):
        # xdes[j] = 3 * cos(Tc * j)
        # yd[j] = 3 * sin(Tc * j)
        # xdes[j] = 3 * cos(Tc * j)
        # ydes[j] = 3 * sin(Tc * j)
        xdes[j] = 3 + 5 * cos(3 * j * Tc)
        ydes[j] = 5 + 5 * sin(3 * j * Tc)
        zdes[j] = 20
        rodes[j] = 5
        alfades[j] = 0
        betades[j] = 0
        qdes[j] = np.array([xdes[j], ydes[j], zdes[j], rof[j], alfaf[j], betaf[j]])
    # Derivada das trajetórias de interesse
    for j in range(0, ciclos, 1):
        # xdesponto[j] = -3 * sin(Tc * j)
        # ydesponto[j] = 3 * cos(Tc * j)
        xdesponto[j] = - 15 * sin(3 * j * Tc)
        ydesponto[j] = 15 * cos(3 * j * Tc)
        zdesponto[j] = 0
        rodesponto[j] = 0
        alfadesponto[j] = 0
        betadesponto[j] = 0
        qpontodes[j] = np.array(
            [xdesponto[j], ydesponto[j], zdesponto[j], rodesponto[j], alfadesponto[j], betadesponto[j]])

########################################################################################################################

# dados do Ponto de controle
a = 0
alfa = pi / 2  #

k1 = 0.8417
k2 = 0.18227
k3 = 0.8354
k4 = 0.17095
k5 = 3.966
k6 = 4.001
k7 = 9.8524
k8 = 4.7295

# Dados do robô no ponto inicial
# Posição inicial do VANT 1
x1[1] = 0
y1[1] = 0
z1[1] = 0.75
psi1[1] = 0
# Posição inicial do VANT 2
x2[1] = -2
y2[1] = 1
z2[1] = 0.75
psi2[1] = 0

# loop de controle
fig = plt.figure(figsize=(60, 60))
ax = fig.add_subplot(111, projection='3d')
ax.view_init(elev=30, azim=45)

for k in range(1, ciclos, 1):
    print('k', k)
    ####################################################################################################################
    # f(x)
    # Inputs X
    # Outputs q(xf,yf,zf,rof,alfaf,betaf)
    xf[k] = x1[k]
    yf[k] = y1[k]
    zf[k] = z1[k]
    rof[k] = sqrt((x2[k] - x1[k]) ** 2 + (y2[k] - y1[k]) ** 2 + (z2[k] - z1[k]) ** 2)
    alfaf[k] = atan((z2[k] - z1[k]) / (sqrt((x2[k] - x1[k]) ** 2 + (y2[k] - y1[k]) ** 2)))
    betaf[k] = atan((y2[k] - y1[k]) / (x2[k] - x1[k]))

    q = np.array([xf[k], yf[k], zf[k], rof[k], alfaf[k], betaf[k]])

    ####################################################################################################################
    # Camada de controle
    # Inputs qdes(xdes,ydes,zdes), qpontodes(xpontodes,ypontodes,zpontodes)
    # Outputs qpontoref(xref, yref, zref)
    # Calculando qtil  qtil[k] = qdes[k] - q[k]
    xtil[k] = (xdes[k] - xf[k])
    ytil[k] = (ydes[k] - yf[k])
    ztil[k] = (zdes[k] - zf[k])
    rotil[k] = (rodes[k] - rof[k])
    alfatil[k] = (alfades[k] - alfaf[k])
    betatil[k] = (betades[k] - betaf[k])

    qtil = np.array([xtil[k], ytil[k], ztil[k], rotil[k], alfatil[k], betatil[k]]).T

    L1 = np.diag([3, 3, 1.9, 3, 3, 1.9])  # Kp
    L2 = np.diag([1, 1, 1, 1, 1, 1])  # kd

    parcial1 = np.dot(L2, qtil)
    # vector = np.vectorize(float)
    # parcial1 = vector(parcial1)
    parcial1[0] = tanh(parcial1[0])
    parcial1[1] = tanh(parcial1[1])
    parcial1[2] = tanh(parcial1[2])
    parcial1[3] = tanh(parcial1[3])
    parcial1[4] = tanh(parcial1[4])
    parcial1[5] = tanh(parcial1[5])

    qpontoref[k] = qpontodes[k] + np.dot(L1, parcial1)

    ####################################################################################################################
    # J^-1(q)
    # Inputs qpontoref(xref, yref, zref), q(xf,yf,zf,rof,alfaf,betaf)
    # Outputs xpontoref

    Jinvq = np.array([[1, 0, 0, 0, 0, 0],
                      [0, 1, 0, 0, 0, 0],
                      [0, 0, 1, 0, 0, 0],
                      [1, 0, 0, np.dot(-cos(alfaf[k]), cos(betaf[k])),
                       np.dot(rof[k], np.dot(sin(alfaf[k]), cos(betaf[k]))),
                       np.dot(rof[k], np.dot(cos(alfaf[k]), sin(betaf[k])))],
                      [0, 1, 0, np.dot(-cos(alfaf[k]), sin(betaf[k])),
                       np.dot(rof[k], np.dot(sin(alfaf[k]), sin(betaf[k]))),
                       np.dot(-rof[k], np.dot(cos(alfaf[k]), cos(betaf[k])))],
                      [0, 0, 1, sin(alfaf[k]), np.dot(rof[k], cos(alfaf[k])), 0]])

    jxqpontoref = np.dot(Jinvq, qpontoref[k])

    xpontoref = np.array([[jxqpontoref[0]],
                          [jxqpontoref[1]],
                          [jxqpontoref[2]],
                          [0],
                          [jxqpontoref[3]],
                          [jxqpontoref[4]],
                          [jxqpontoref[5]],
                          [0]])
    ####################################################################################################################
    # K^-1
    # Inputs xpontoref
    # Outputs vref

    kinv = np.array([[cos(psi1[k]), sin(psi1[k]), 0, 0, 0, 0, 0, 0],
                     [-sin(psi1[k]), cos(psi1[k]), 0, 0, 0, 0, 0, 0],
                     [0, 0, 1, 0, 0, 0, 0, 0],
                     [0, 0, 0, 1, 0, 0, 0, 0],
                     [0, 0, 0, 0, cos(psi2[k]), sin(psi2[k]), 0, 0],
                     [0, 0, 0, 0, -sin(psi2[k]), cos(psi2[k]), 0, 0],
                     [0, 0, 0, 0, 0, 0, 1, 0],
                     [0, 0, 0, 0, 0, 0, 0, 1]])

    Vref = np.dot(kinv, xpontoref)

    ####################################################################################################################

    vxref1[k] = Vref[0]
    vyref1[k] = Vref[1]
    vzref1[k] = Vref[2]
    vpsi1[k] = Vref[3]

    vxref2[k] = Vref[4]
    vyref2[k] = Vref[5]
    vzref2[k] = Vref[6]
    vpsi2[k] = Vref[7]

    ####################################################################################################################
    # Cinematica do Drone 1
    # inputs vxref, vyref1, vzref1, vpsi1
    # outputs x, y, z, psi
    # Calcular Xr, Yr e Psir (Ponto de controle)
    mie1 = np.array([vxref1[k], vyref1[k], vzref1[k], vpsi1[k]])

    Ae1 = np.array([[k1 * cos(psi1[k]), -k3 * sin(psi1[k]), 0, 0],
                    [k1 * sin(psi1[k]), +k3 * cos(psi1[k]), 0, 0],
                    [0, 0, k5, 0],
                    [0, 0, 0, k7]])
    # Ae1 = np.array([[cos(psi1[k]), sin(psi1[k]), 0, 0],
    #                 [sin(psi1[k]), cos(psi1[k]), 0, 0],
    #                 [0, 0, 1, 0],
    #                 [0, 0, 0, 1]])

    px1 = np.dot(Ae1, mie1)

    xponto1[k] = px1[0]  # Saida xponto
    yponto1[k] = px1[1]  # Saida yponto
    zponto1[k] = px1[2]  # Saida yponto
    psiponto1[k] = px1[3]  # Saida psip

    if k < ciclos - 1:
        x1[k + 1] = Tc * xponto1[k] + x1[k]  # x deve ser x
        y1[k + 1] = Tc * yponto1[k] + y1[k]
        z1[k + 1] = Tc * zponto1[k] + z1[k]
        psi1[k + 1] = Tc * psiponto1[k] + psi1[k]  # Saida p

    ####################################################################################################################
    # Cinematica do Drone 2
    # inputs vxref, vyref, vzref, vpsi
    # outputs x, y, z, psi
    # Calcular Xr, Yr e Psir (Ponto de controle)
    mie2 = np.array([vxref2[k], vyref2[k], vzref2[k], vpsi2[k]])

    # Ae2 = np.array([[cos(psi1[k]), sin(psi1[k]), 0, 0],
    #                 [sin(psi1[k]), cos(psi1[k]), 0, 0],
    #                 [0, 0, 1, 0],
    #                 [0, 0, 0, 1]])

    Ae2 = np.array([[k1 * cos(psi1[k]), -k3 * sin(psi1[k]), 0, 0],
                    [k1 * sin(psi1[k]), +k3 * cos(psi1[k]), 0, 0],
                    [0, 0, k5, 0],
                    [0, 0, 0, k7]])

    px2 = np.dot(Ae2, mie2)

    xponto2[k] = px2[0]  # Saida xponto
    yponto2[k] = px2[1]  # Saida yponto
    zponto2[k] = px2[2]  # Saida yponto
    psiponto2[k] = px2[3]  # Saida psip

    if k < ciclos - 1:
        x2[k + 1] = Tc * xponto2[k] + x2[k]  # x deve ser x
        y2[k + 1] = Tc * yponto2[k] + y2[k]
        z2[k + 1] = Tc * zponto2[k] + z2[k]
        psi2[k + 1] = Tc * psiponto2[k] + psi2[k]  # Saida p

    ####################################################################################################################
    # Preparação dos dados para plotar
    erro[k] = sqrt((xtil[k]) ** 2 + (ytil[k]) ** 2 + (ztil[k]) ** 2)
    plt.axis('auto')
    ####################################################################################################################

    # Plot saida
    plt.cla()
    # ax.scatter(xf[0:k], yf[0:k], zf[0:k])

    ax.plot(xdes[1:k], ydes[1:k], zdes[1:k], "-r", label="Trajetoria desejada")
    plt.plot(xdes[k], ydes[k], zdes[k], "xr")

    ax.plot(x1[1:k], y1[1:k], z1[1:k], "-g", label="Trajetoria Vant 1")
    plt.plot(x1[k], y1[k], z1[k], "xg")

    ax.plot(x2[1:k], y2[1:k], z2[1:k], "-b", label="Trajetoria Vant 2")
    plt.plot(x2[k], y2[k], z2[k], "xb")

    # ax.set_aspect('equal')
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    ax.legend(loc='best')

    # ax.view_init(elev=-75,azim=-90)
    plt.pause(0.01)
###################################


# plt.show()

# Plota dadas
plt.show(block=False)
# Trajetoria desejada - x, y , z
# Trajetoria executada - x, y, z
fig2, axs = plt.subplots(3, 1, num='Trajetoria desejada do robô no tempo')
axs[0].plot(xdes, label="Trajetoria desejada - Eixo x")
axs[0].plot(x1, label="Trajetoria executada - Eixo x")
axs[0].set_title('Trajetoria eixo x')
axs[0].set(ylabel='[m]')
axs[0].legend(loc='best')
axs[0].grid()

axs[1].plot(ydes, label="Trajetoria desejada - Eixo y")
axs[1].set_title('Trajetoria eixo y')
axs[1].plot(y1, label="Trajetoria executada - Eixo y")
axs[1].set_title('Trajetoria eixo y')
axs[1].set(xlabel='[Ciclos de 100 ms]', ylabel='[m]')
axs[1].legend(loc='best')
axs[1].grid()

axs[2].plot(zdes, label="Trajetoria desejada - Eixo z")
axs[2].set_title('Trajetoria eixo z')
axs[2].plot(z1, label="Trajetoria executada - Eixo z")
axs[2].set_title('Trajetoria eixo z')
axs[2].set(xlabel='[Ciclos de 100 ms]', ylabel='[m]')
axs[2].legend(loc='best')
axs[2].grid()

plt.tight_layout()

# Velocidade comandada - Ucref, Wcref
# Velocidade real  - U, W
fig3, axs = plt.subplots(4, 1, num='Velocidade do robô no tempo')
axs[0].plot(vxref1, label="Veloc. linear comandada")
axs[0].plot(vxref1, label="Veloc. linear x real - Vant1")
axs[0].plot(vxref2, label="Veloc. linear x real - Vant2")
axs[0].set_title('Velocidade linear comandada x real ')
axs[0].set(ylabel='[m/s]')
axs[0].legend(loc='best')
axs[0].grid()

axs[1].plot(vyref1, label="Veloc. linear comandada")
axs[1].plot(vyref1, label="Veloc. linear y real - Vant1")
axs[1].plot(vyref2, label="Veloc. linear y real - Vant2")
axs[1].set_title('Velocidade linear comandada x real')
axs[1].set(xlabel='[Ciclos de 100 ms]', ylabel='[rad/s]')
axs[1].legend(loc='best')
axs[1].grid()

axs[2].plot(vzref1, label="Veloc. linear comandada")
axs[2].plot(vzref1, label="Veloc. linear z real - Vant1")
axs[2].plot(vzref2, label="Veloc. linear z real - Vant2")
axs[2].set_title('Velocidade linear comandada x real')
axs[2].set(xlabel='[Ciclos de 100 ms]', ylabel='[rad/s]')
axs[2].legend(loc='best')
axs[2].grid()

axs[3].plot(vpsi1, label="Veloc. angular comandada")
axs[3].plot(vpsi1, label="Veloc. angular real")
axs[3].plot(vpsi2, label="Veloc. angular real")
axs[3].set_title('Velocidade angular comandada x real')
axs[3].set(xlabel='[Ciclos de 100 ms]', ylabel='[rad/s]')
axs[3].legend(loc='best')
axs[3].grid()

plt.tight_layout()

# Pose do robô - x, y, psi
fig4, axs = plt.subplots(4, 1, num='Posição do robô no tempo')
axs[0].plot(x1, label="posição x do Vant1")
axs[0].plot(x2, label="posição x do Vant2")
axs[0].set_title('posição x do robô no tempo')
axs[0].set(ylabel='y[m]')
axs[0].legend(loc='best')
axs[0].grid()

axs[1].plot(y1, label="posição y do Vant1")
axs[1].plot(y2, label="posição y do Vant2")
axs[1].set_title('posição y do robô no tempo')
axs[1].set(ylabel='y[m]')
axs[1].legend(loc='best')
axs[1].grid()

axs[2].plot(z1, label="posição z do Vant1")
axs[2].plot(z2, label="posição z do Vant2")
axs[2].set_title('posição z do robô no tempo')
axs[2].set(ylabel='y[m]')
axs[2].legend(loc='best')
axs[2].grid()

axs[3].plot(psi1, label="Orientação do Vant1")
axs[3].plot(psi2, label="Orientação do Vant2")
axs[3].set_title('Angulo psi do robô no tempo')
axs[3].set(xlabel='[Ciclos de 100 ms]', ylabel='radianos')
axs[3].grid()
axs[3].legend(loc='best')
plt.tight_layout()

# Erro de trajetoria - erro
plt.figure('Erro de trajetória')
plt.plot(erro, label="Erro de trajetória")
plt.title('Erro de trajetória')
plt.xlabel('[Ciclos de 100 ms]')
plt.ylabel('[m]')
plt.legend(loc='best')
plt.grid(True)

# q = np.array([xf[k], yf[k], zf[k], rof[k], alfaf[k], betaf[k]])
#
# # Dados da formação - xf, yf, zf, rof, alfaf, betaf
# fig4, axs = plt.subplots(4, 1, num='Posição do robô no tempo')
# axs[0].plot(xf, label="posição x do Vant1")
# axs[0].set_title('posição x do robô no tempo')
# axs[0].set(ylabel='y[m]')
# axs[0].legend(loc='best')
# axs[0].grid()
#
# axs[1].plot(y1, label="posição y do Vant1")
# axs[1].plot(y2, label="posição y do Vant2")
# axs[1].set_title('posição y do robô no tempo')
# axs[1].set(ylabel='y[m]')
# axs[1].legend(loc='best')
# axs[1].grid()
#
# axs[2].plot(z1, label="posição z do Vant1")
# axs[2].plot(z2, label="posição z do Vant2")
# axs[2].set_title('posição z do robô no tempo')
# axs[2].set(ylabel='y[m]')
# axs[2].legend(loc='best')
# axs[2].grid()
#
# axs[3].plot(psi1, label="Orientação do Vant1")
# axs[3].plot(psi2, label="Orientação do Vant2")
# axs[3].set_title('Angulo psi do robô no tempo')
# axs[3].set(xlabel='[Ciclos de 100 ms]', ylabel='radianos')
# axs[3].grid()
# axs[3].legend(loc='best')
# plt.tight_layout()

# # Rascunho
# fig6, axs = plt.subplots(4, 1, num='RASCUNHO ')
# plt.tight_layout()
# axs[0].plot(psi)
# axs[0].set_title('psi')
# axs[1].plot(psitil)
# axs[1].set_title('psitil')
# axs[2].plot(psid)
# axs[2].set_title('psid')
# axs[3].plot(der_psi)
# axs[3].set_title('der_psi')
#
# # rascunho 2
# fig7, axs = plt.subplots(3, 1, num='RASCUNHO2 ')
# plt.tight_layout()
# axs[0].plot(psiponto)
# axs[0].set_title('psiponto')
# axs[1].plot(xponto)
# axs[1].set_title('Xponto')
# axs[2].plot(yponto)
# axs[2].set_title('Yponto')
plt.show()
