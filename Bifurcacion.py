import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from functools import partial
from multiprocessing import Pool
import time

# Definir constantes globales
g = 981  # aceleración gravitacional en cm/s^2
a = 0.25  # radio
gamma = 0.001  # coeficiente de amortiguamiento
m0 = 2.5  # masa inicial
z0 = 0  # altura inicial
z_crit = 5e1  # altura crítica

# Definir k(m) como una función
def k_m(m):
    if m < 4.61:
        k_val = -11.4 * m + 52.5
    else:
        k_val = 0
    return k_val

# Sistema de ODEs
def system(t, y, v0):
    z, dzdt, m = y
    k = k_m(m)
    flow_rate = np.pi * a**2 * v0
    dmdt = flow_rate
    dz2dt2 = (m * g - k * z - gamma * dzdt - flow_rate * (dzdt - v0)) / m
    return [dzdt, dz2dt2, dmdt]

# Función de evento
def make_event_breakoff():
    def event_breakoff(t, y):
        z = y[0]
        return z - z_crit  # El evento se activa cuando z = z_crit
    event_breakoff.terminal = True  # Detener la integración
    event_breakoff.direction = 0    # Activar en ambas direcciones
    return event_breakoff

# Función para calcular los períodos para un V0 dado
def compute_periods_for_V0(V0):
    v0 = V0
    initial_conditions = [z0, 0, m0]
    periods = []
    N_per = 50  # Número total de períodos
    N_start = 30  # Comenzar a recopilar datos después de este período para evitar transitorios

    for i in range(N_per):
        t_span = [0, 1000]
        sol = solve_ivp(
            fun=partial(system, v0=v0),
            t_span=t_span,
            y0=initial_conditions,
            events=make_event_breakoff(),
            rtol=1e-8,
            atol=1e-10
        )
        periods.append(round(sol.t_events[0][0], 8))
        # Actualizar las condiciones iniciales para el siguiente ciclo
        m_new = sol.y[2][-1] * 0.7
        initial_conditions = [0, 0, m_new]

    # Recopilar T_n después de los transitorios
    V0_list = [V0] * (N_per - N_start)
    Tn_list = periods[N_start:]
    return V0_list, Tn_list

# Lista de valores de V0
V0_values = np.arange(0.02, 0.1, 0.0001)  # Ajusta el rango y el paso según sea necesario

# Ejecutar en paralelo usando multiprocessing

time_start = time.time()
if __name__ == "__main__":
    with Pool() as pool:
        results = pool.map(compute_periods_for_V0, V0_values)

    # Unir los resultados
    V0_list = []
    Tn_list = []
    for V0_sublist, Tn_sublist in results:
        V0_list.extend(V0_sublist)
        Tn_list.extend(Tn_sublist)

    # Graficar el diagrama de bifurcación
    plt.figure()
    plt.scatter(V0_list, Tn_list, s=10, marker='.')
    #xlog
    plt.yscale('log')
    plt.xlabel(r'$V_0$')
    plt.ylabel(r'Período $T_n$')
    plt.title(r'Diagrama de Bifurcación: $T_n$ vs $V_0$')
    plt.grid(True)
    plt.savefig('bifurcation_diagram9.pdf')
    plt.show()
    print('Tiempo de ejecución:', time.time() - time_start, 'segundos')