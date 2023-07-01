import numpy as np
from collections import defaultdict
import re
import pandas as pd

# Функция бинома Ньютона
def binomial_coefficients(n):
    coefficients = [1]
    for i in range(1, n + 1):
        coefficient = coefficients[-1] * (n - i + 1) // i
        coefficients.append(coefficient)
    return coefficients

# Набор для упорядочивания спектра
class Spectrum:
    def __init__(self):
        self.data = []

    def get_array(self):
        x = []
        y = []
        for elem in self.data:
            x.append(elem.x)
            y.append(elem.y)
        return np.array(x), np.array(y)

    def sort_data(self):
        self.data.sort(key=lambda x: x.x)


# Класс для одной записи спектра
class Data_spect:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def set(self, data):
        self.x, self.y = data

# Учитываем мультиплетность сигнала
class Signal:
    def __init__(self, center, freq, coupling, numbers_atoms, base_intensity):
        self.center = center
        self.freq = freq
        self.base_intensity = base_intensity
        self.coupling = [i for i in zip(coupling, numbers_atoms)]
        self.binom = {1: [1, 1], 2: [1, 2, 1], 3: [1, 3, 3, 1], 4: [1, 4, 6, 4, 1],
                      5: [1, 5, 10, 10, 5, 1], 6: [1, 6, 15, 20, 15, 6, 1],
                      7: [1, 7, 21, 35, 35, 21, 7, 1]}
        self.signals = self.refine()

    def cleaner(self):
        buffer = []
        for i in self.coupling:
            if i[0] > 2:
                buffer.append(i)
        self.coupling = buffer.copy()
    def refine(self):
        self.cleaner()
        self.coupling.sort(key=lambda x: x[0])
        centers_peaks = []
        if len(self.coupling) == 0:
            return [(self.center, 1 * self.base_intensity),]
        for idx, elem in enumerate(self.coupling):
            if idx == 0:
                centers_peaks = self.multi(self.center, self.base_intensity, elem[0], elem[1] + 1)
            else:
                buffer = []
                while len(centers_peaks) != 0:
                    peak = centers_peaks.pop()
                    out = self.multi(peak[0], peak[1], elem[0], elem[1] + 1)
                    for peaks in out:
                        buffer.append(peaks)
                centers_peaks = buffer.copy()
        return centers_peaks.copy()


    def multi(self, center: float, intensity: float, J: float, count_peaks: int):
        indent = J / self.freq
        coord_peaks = []

        # Добавление координат пиков по ppm
        if count_peaks % 2 == 1 and count_peaks != 1:
            for i in range(int(count_peaks / 2)):
                coord_peaks.append(center + (i + 1) * indent)
                coord_peaks.append(center - (i + 1) * indent)
            coord_peaks.append(center)
        elif count_peaks % 2 == 0:
            coord_peaks.append(center - 0.5 * indent)
            coord_peaks.append(center + 0.5 * indent)
            if count_peaks > 2:
                for i in range(int((count_peaks - 1) / 2)):
                    coord_peaks.append(center + (i + 1.5) * indent)
                    coord_peaks.append(center + (i - 1.5) * indent)
        else:
            coord_peaks.append(center)

        coord_peaks.sort()

        # Добавление интенсивности пиков
        if count_peaks == 1:
            return [(coord_peaks[-1], 1 * intensity), ]
        else:
            result = []
            intensity_arr = binomial_coefficients(count_peaks - 1)
            max_intensity = intensity / max(intensity_arr)
            for idx, idy in zip(coord_peaks, intensity_arr):
                result.append((idx, idy * max_intensity))
            return result

class JShiftsHandler:
    def __init__(self, filename, groups):
        self.ref_h = 31.660
        self.ref_c = 182.219
        self.ref_f = 160.636
        self.ref_o = 0
        self.ref_n = 0

        self.filename = filename
        self.groups = groups
        self.full_coupling = self.get_full_table()
        self.j_groups = self.groups_data()
        self.full_shifts = self.get_shifts()
        self.groups_shits = self.get_group_shift()
        self.numbers_atoms = self.get_numbers_atoms()

    def get_full_table(self):
        text_table = ''
        with open(self.filename, 'r') as file:
            write_table = False
            for line in file.readlines():
                if 'SUMMARY OF ISOTROPIC COUPLING CONSTANTS' in line:
                    write_table = True
                elif 'Maximum memory used throughout ' in line and write_table:
                    write_table = False

                if write_table:
                    text_table += line

        clear_text = text_table.split('\n')[2:]
        arr = []
        for i in clear_text:
            if '                ' in i:
                pattern = r'\d+\s*[a-zA-Z]+'
                result = re.findall(pattern, i)
                arr.append(result)
            else:
                arr.append(i.split())

        arr = list(filter(lambda x: len(x) != 0, arr))

        data = defaultdict(list)
        for i in arr:
            if i[0].isdigit():
                data[str(int(i[0]) + 1)] += [float(x) for x in i[2:]]
        df = pd.DataFrame(data=data, index=data.keys())
        return df

    def groups_data(self):
        READY_TABLE = []
        for idx1, i in enumerate(self.groups):
            READY_TABLE.append([])
            for idx2, j in enumerate(self.groups):
                filter_row = self.full_coupling.loc[j]
                filter_column = filter_row.filter(items=[str(x) for x in i])
                READY_TABLE[idx1].append(round(np.array(filter_column).mean(), 3))

        return pd.DataFrame(READY_TABLE)

    def get_shifts(self):
        data = []
        with open(self.filename, 'r') as file:
            reading = False
            for line in file.readlines():
                if '  Nucleus  Element    Isotropic     Anisotropy' in line:
                    reading = True
                if line == '\n':
                    reading = False
                if reading:
                    data.append(line.split())
            data.pop(1)
        df = pd.DataFrame(columns=data[0], data=data[1:])
        df['Nucleus'] = df['Nucleus'].astype(int)
        df['Isotropic'] = df['Isotropic'].astype(float)
        df['Anisotropy'] = df['Anisotropy'].astype(float)

        df['Nucleus'] += 1

        # Вносим разницу между референсом и полученными значениями
        df.loc[df['Element'] == 'H', 'Isotropic'] -= self.ref_h
        df.loc[df['Element'] == 'H', 'Isotropic'] *= -1

        df.loc[df['Element'] == 'O', 'Isotropic'] -= self.ref_o
        df.loc[df['Element'] == 'O', 'Isotropic'] *= -1

        df.loc[df['Element'] == 'C', 'Isotropic'] -= self.ref_c
        df.loc[df['Element'] == 'C', 'Isotropic'] *= -1

        df.loc[df['Element'] == 'F', 'Isotropic'] -= self.ref_f
        df.loc[df['Element'] == 'F', 'Isotropic'] *= -1

        return df

    def get_numbers_atoms(self):
        GROUPS = {}
        for idx, i in enumerate(self.groups):
            GROUPS[str(idx + 1)] = len(i)
        return GROUPS


    def get_group_shift(self):
        READY_GROUP_SHIFT = {}
        for idx1, i in enumerate(self.groups):
            filter_row = self.full_shifts[self.full_shifts['Nucleus'].isin([int(x) for x in i])]
            READY_GROUP_SHIFT[str(idx1 + 1)] = round(np.array(filter_row['Isotropic']).mean(), 3)

        return pd.DataFrame(data=READY_GROUP_SHIFT, index=[0])