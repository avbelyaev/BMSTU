import math
import random
import matplotlib.pyplot as plt

FUNC = lambda x: math.sin(x) + 0.5  # исходая функция
K_STEPS = 100  # количество отсчетов
R_WINDOW_SIZE = 3  # размер скользящего окна
SIGMA_NOISE_AMPLITUDE = 0.25  # разброс значений зашемленного сигнала
X_MIN = 0  # интервал от
X_MAX = math.pi  # интервал до


# шум с нормальным распределением
def uniform_noise(signal: float) -> float:
    noise = random.uniform(-SIGMA_NOISE_AMPLITUDE, SIGMA_NOISE_AMPLITUDE)
    return signal + noise


# Среднее гармоническое
# def harmonic_mean(k: int, m: int, alphas: list) -> float:

# Чебышева


class Signal:
    def __init__(self, func: callable, xmin: float, xmax: float, steps: int,
                 window_size: int):
        self.signal_function = func
        self.x_min = xmin
        self.x_max = xmax
        self.steps = steps
        self.window_size = window_size
        self.m = (window_size - 1) // 2

    @property
    def xs(self) -> list:
        xs = []
        x_curr = self.x_min
        step = (self.x_max - self.x_min) / self.steps
        while x_curr <= self.x_max:
            xs.append(x_curr)
            x_curr += step
        return xs

    @property
    def clear_signal(self) -> list:
        return list(map(lambda x: self.signal_function(x), self.xs))

    @property
    def noisy_signal(self) -> list:
        return list(map(lambda signal: uniform_noise(signal), self.clear_signal))

    @property
    def filtered_signal(self) -> list:
        fs = self.noisy_signal
        alphas = [0.4, 0.2, 0.4]  # TODO wtf is this

        filtered = []
        filtered.extend(self.clear_signal[:self.m])

        k = self.m
        while k < self.steps - self.m:
            window_sum = sum([alphas[0] / fs[k - 1],
                              alphas[1] / fs[k],
                              alphas[2] / fs[k + 1]])
            harmonic_mean_value = 1 / window_sum
            filtered.append(harmonic_mean_value)
            k += 1
        filtered.extend((self.clear_signal[k:]))
        return filtered

    def plot(self):
        plt.subplot(311)
        plt.title('Clear signal', fontdict={'size': 10})
        plt.plot(self.xs, self.clear_signal, c='green')

        plt.subplot(312)
        plt.title('Noisy signal', fontdict={'size': 10})
        plt.plot(self.xs, self.noisy_signal, c='red')

        plt.subplot(313)
        plt.title('Filtered signal', fontdict={'size': 10})
        plt.plot(self.xs, self.filtered_signal, c='blue')

        plt.tight_layout()
        plt.savefig('plot.png')


if __name__ == '__main__':
    Signal(
        func=FUNC,
        xmin=X_MIN,
        xmax=X_MAX,
        steps=K_STEPS,
        window_size=R_WINDOW_SIZE
    ).plot()
