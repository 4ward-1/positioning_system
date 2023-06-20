import random
import matplotlib.pyplot as plt
import numpy as np

def constants():
    const = {'Fmax': 1000, 'fn': 50, 'L_B': 10, 'L_D': 20, 'm': 10000, 'W_D': 5, 'tao_s': 1,
            'k_tr': 10000*(1 / 5), 'DT': 0.01, 'Kp': 100, 'Ki': 0, 'Kd': 60, 'err_poz': 0.1}
    # good podbor: 1) 200, 10, 1500; 2) 500, 0, 2500
    # Ziver - Nikols method: 1) 1740, 246*0, 3078
    # good from Ziver - Nikols: 1) 540, 0, 3500
    # with Td: 1) 100, 0, 60
    # method LACH при использовании Td в дифференциировании
    # Пояснение: Fmax - макс допустимое усилие в лебедке, Н; fn - частота помехи датчика положения, Гц; L_B - длина яхты, м
    # L_D - расстояние между ЭМП лебедок 1-2(1-3) в координате x, м; m - масса яхты, кг; W_D - расстояние между ЭМП лебедок 2-3 в, м
    # tao_s - постоянная времени измерения положения яхты, с; k_tr - # коэффициент силы трения; DT - дискретность
    # Kp,Ki,Kd - коэффициенты ПИД - регулятора, err_poz - допустимая ошибка позиционроования
    const['Td'] = 0.15*const['Kd']  # постоянная времени дифференциатора для фильтрации ВЧ помехи
    
    # случайное начальное положение яхты (точки установки датчика положения)
    const['x_0']    = -11 #-11 #round(random.uniform((const['L_B'] / 2) + 0.1, (const['L_D'] - (const['L_B'] / 2))), 2)
    # случайная координата x, в которую нужно привести яхту (точку установки датчика положения)
    const['x_treb'] = -6 #-6 #round(random.uniform((const['L_B'] / 2) + 0.1, (const['L_D'] - (const['L_B'] / 2))), 2)
    return const

def DP_izm(x, x_iner_p, tao_s, DT):
    x_iner = (x + (tao_s / DT) * x_iner_p) / (1 + (tao_s / DT))  # апериодическое звено первого порядка
    return x_iner

def pid_regulator(dx_i, dx_i_p, Kp, Ki, Kd, Fd_i_p, Td, DT, Fmax):
    Fp_i = Kp * dx_i  # Сигнал, пропорциональный ошибке позиционирования
    Fd_i = (Td*Fd_i_p + Kd * (dx_i - dx_i_p))/(Td + DT)# Сигнал, пропорциональный производной ошибки позиционирования
    Fi_i = Ki * (dx_i_p + DT * dx_i)  # Сигнал, пропорциональный интегралу ошибки позиционирования
    Fi   = Fp_i + Fd_i + Fi_i

    if (abs(Fi) >= Fmax):
        if (Fi > 0):
            Flim = Fmax
        else:
            Flim = -Fmax
    else: Flim = Fi

    return Flim, Fd_i

def Fx_proec(F1, F2, F3, x, L_D, L_B, W_D):
    y = L_D - (abs(x) + (L_B / 2))  # расстояние между ЭМП2(3) и местом крепления лебедок 2,3 на яхте
    al = np.arctan((W_D / 2) / y)  # угол между осью x и направлением действия сил F2, F3
    F_summ_x = -(F1 + F2 * np.cos(al) + F3 * np.cos(al))
    return F_summ_x  # суммарная сила тяги, действующая на яхту вдоль оси x

def yaht_moving(F_summ_x, m, k_tr, DT, x_p, Vx_p):
    Vx = (m*Vx_p + F_summ_x*DT)/(m + k_tr*DT)
    x = x_p + DT*Vx
    return x, Vx

if __name__ == "__main__":
    const = constants()  # Характеристики дока, яхты и ПИД-регулятора системы позиционирования яхты
    dx_values   = [] # Ошибка позиционирования яхты
    time_values = [] # Время
    F1_val = [] # усилие 1
    F2_val = [] # усилие 2
    F3_val = [] # усилие 3
    x_val = []
    i = time = 0
    x = const['x_0']
    dx_izm = const['x_treb'] - const['x_0']
    while (time <= 200): #(abs(dx) >= const['err_poz']) and (time <= 300):
        time_values.append(time)
        # Начальные условия:
        if i == 0:
            x_iner_p = const['x_0'] # начальное измерение датчика положения
            x_p  = const['x_0'] # начальное положение яхты
            Vx_p = 0    # начальная скорость яхты
            Fd_1_p = 0
            Fd_2_p = 0
            Fd_3_p = 0
            if (dx_izm <= 0):
                dx_1_p = 0
                dx_2_p = -dx_izm / 2
                dx_3_p = -dx_izm / 2
            else:
                dx_1_p = -dx_izm
                dx_2_p = 0
                dx_3_p = 0
        time = i * const['DT']  # время
        # измерения датчика положения яхты
        x_iner = DP_izm(x, x_iner_p, const['tao_s'], const['DT'])
        x_iner_p = x_iner
        pom = 0.2 * np.sin(2 * np.pi * const['fn'] * time)
        x_izm = x_iner + pom
        # разность между заданной и измеренной координатой x яхты:
        dx_izm = const['x_treb'] - x_izm
        # Выбор управляющих ЭМП(1,2,3) в зависимости от знака dx
        if (dx_izm <= 0):
            dx_1 = 0
            dx_2 = -dx_izm / 2
            dx_3 = -dx_izm / 2
        else:
            dx_1 = -dx_izm
            dx_2 = 0
            dx_3 = 0

        # Формирование управляющих воздействий в виде силы тяги ЭМП(1,2,3) с использованием ПИД-регулятора и ограничением усилия величиной Fmax
        [F1, Fd_1] = pid_regulator(dx_1, dx_1_p, const['Kp'], const['Ki'], const['Kd'], Fd_1_p, const['Td'], const['DT'], const['Fmax'])
        [F2, Fd_2] = pid_regulator(dx_2, dx_2_p, const['Kp'], const['Ki'], const['Kd'], Fd_2_p, const['Td'], const['DT'], const['Fmax'])
        [F3, Fd_3] = pid_regulator(dx_3, dx_3_p, const['Kp'], const['Ki'], const['Kd'], Fd_3_p, const['Td'], const['DT'], const['Fmax'])
        Fd_1_p = Fd_1
        Fd_2_p = Fd_2
        Fd_3_p = Fd_3

        dx_1_p = dx_1
        dx_2_p = dx_2
        dx_3_p = dx_3

        # суммарная сила тяги, действующая на яхту вдоль оси x
        x_val.append(x)
        F_summ_x = Fx_proec(F1, F2, F3, x, const['L_D'], const['L_B'], const['W_D'])
        # передаточная функция яхты на воде от силы тяги к перемещению x и скорости Vx
        [x, Vx] = yaht_moving(F_summ_x, const['m'], const['k_tr'], const['DT'], x_p, Vx_p)
        x_p = x
        Vx_p = Vx
        dx_values.append(const['x_treb'] - x)
        F1_val.append(F1)
        F2_val.append(F2)
        F3_val.append(F3)
        i += 1

    plt.figure(1)
    plt.plot(time_values, dx_values)
    plt.grid(color = 'red', linewidth = 1, linestyle = '--')
    plt.xlabel('t, с')
    plt.ylabel('dx, м')
    plt.title('Ошибка позиционирования яхты')
   # plt.show()

    plt.figure(2)
    plt.plot(time_values, F1_val)
    plt.grid(color = 'red', linewidth = 1, linestyle = '--')
    plt.xlabel('t, с')
    plt.ylabel('F1, H')
    plt.title('Усилие на 1й лебедке')
    
    plt.figure(3)
    plt.plot(time_values, F2_val)
    plt.grid(color = 'red', linewidth = 1, linestyle = '--')
    plt.xlabel('t, с')
    plt.ylabel('F2, H')
    plt.title('Усилие на 2й лебедке')
    
    plt.figure(4)
    plt.plot(time_values, F3_val)
    plt.grid(color = 'red', linewidth = 1, linestyle = '--')
    plt.xlabel('t, с')
    plt.ylabel('F3, H')
    plt.title('Усилие на 3й лебедке')

    plt.figure(5)
    plt.plot(time_values, x_val)
    plt.grid(color='red', linewidth=1, linestyle='--')
    plt.xlabel('t, с')
    plt.ylabel('x, м')
    plt.title('Положение яхты')

    plt.show()

    # Пояснение
    # 1) Нашел ошибку в исходном коде - в строке 120 положение яхты записывалось как dx_values.append(x - const['x_treb']),
    #    а должно было быть dx_values.append(const['x_treb'] - x);
    # 2) В текущем коде принял направление оси x от стенки дока влево с уменьшеньшением значения x(как на рисунке),
    #    в связи с этим изменил знаки условий в строках 76 и 93 и знаки в 78, 79, 81, 95, 96, 98;
    # 3) Увеличил коэффициент усиления дифференциатора Kd до 60, при этом всплески в усилиях ЭМП исчезли;
    # 4) Fd_1, Fd_2 и Fd_3 - это значения выхода дифференциатора с постоянной времени на предыдущем шаге, а усилия на лебедках обозначены
    #    как F1,F2,F3. В связи с этим изменил значения в строках 122,123,124;
    # 5) Добавил график положения яхты