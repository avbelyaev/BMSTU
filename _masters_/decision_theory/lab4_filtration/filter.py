import math
import random
import matplotlib.pyplot as plt
import numpy as np

# @formatter:off
SIGNAL_FUNC = lambda x: math.sin(x) + 0.5       # исходая функция
K_STEPS = 100                                   # количество отсчетов
R_WINDOW_SIZE = 5                               # размер скользящего окна
SIGMA_NOISE_AMPLITUDE = 0.25                    # разброс значений зашемленного сигнала
X_MIN = 0                                       # интервал от
X_MAX = math.pi                                 # интервал до
NAIVE_CENTRAL_ALPHA = 0.5                       # наивное значение alpha
NAIVE_LAMBDA = 0.5                              # наивное значение лямбды
# @formatter:on


# шум с нормальным распределением
def uniform_noise(signal: float) -> float:
    noise = random.uniform(-SIGMA_NOISE_AMPLITUDE, SIGMA_NOISE_AMPLITUDE)
    return signal + noise


# критерий зашумленности сигнала, метрика чебышева
# aka максимальная разность двух рядомстоящих элементов
def noisy_criteria(signal: list) -> float:
    max_diff = -999999
    for k in range(1, K_STEPS):
        diff = abs(signal[k] - signal[k - 1])
        if diff > max_diff:
            max_diff = diff
    return max_diff


# критейрий близости сигналов, метрика чебышева
def closeness_criteria(signal1: list, signal2: list) -> float:
    max_diff = -999999
    for k in range(0, K_STEPS):
        diff = abs(signal1[k] - signal2[k])
        if diff > max_diff:
            max_diff = diff
    return max_diff


# метрика расстояния до идеальной точки, метрика чебышева
def ideal_point_distance_metric(omega: float, delta: float) -> float:
    return max(omega, delta)


class Signal:
    def __init__(self, func: callable):
        self.signal_function = func
        self.m = (R_WINDOW_SIZE - 1) // 2
        self.central_alpha = NAIVE_CENTRAL_ALPHA
        self.lambda_weight = NAIVE_LAMBDA

    @property
    def xs(self):
        xs = []
        x_curr = X_MIN
        step = (X_MAX - X_MIN) / K_STEPS
        while x_curr <= X_MAX:
            xs.append(x_curr)
            x_curr += step
        return xs

    def clear_signal(self) -> list:
        return list(map(lambda x: self.signal_function(x), self.xs))

    def noisy_signal(self) -> list:
        return list(map(lambda signal: uniform_noise(signal), self.clear_signal()))

    def filtered_signal(self, signal: list) -> list:
        filtered = []
        filtered.extend(signal[:self.m])  # TODO по краям - изначальный или шумный?

        # @formatter:off
        if 3 == R_WINDOW_SIZE:
            alphas = [
                (1 - self.central_alpha) / 2,       # i - 1
                self.central_alpha,                 # i
                (1 - self.central_alpha) / 2        # i + 1
            ]
        elif 5 == R_WINDOW_SIZE:
            alphas = [
                (1 - self.central_alpha) / 4,       # i - 2
                (1 - self.central_alpha) / 4,       # i - 1
                self.central_alpha,                 # i
                (1 - self.central_alpha) / 4,       # i + 1
                (1 - self.central_alpha) / 4        # i + 2
            ]
        else:
            raise ValueError('invalid window size')
        # @formatter:on

        k = self.m
        while k < K_STEPS - self.m:
            # Среднее гармоническое
            if 3 == R_WINDOW_SIZE:
                window_sum = sum([alphas[0] / signal[k - 1],
                                  alphas[1] / signal[k],
                                  alphas[2] / signal[k + 1]])
            else:
                window_sum = sum([alphas[0] / signal[k - 2],
                                  alphas[1] / signal[k - 1],
                                  alphas[2] / signal[k],
                                  alphas[3] / signal[k + 1],
                                  alphas[4] / signal[k + 2]])
            harmonic_mean_value = 1 / window_sum
            filtered.append(harmonic_mean_value)
            k += 1
        filtered.extend(signal[k:])
        return filtered

    def linear_convolution_of_criteria(self):
        omega = noisy_criteria(self.filtered_signal(self.noisy_signal()))
        delta = closeness_criteria(self.filtered_signal(self.noisy_signal()), self.clear_signal())
        return self.lambda_weight * omega + (1 - self.lambda_weight) * delta

    def optimize_alphas(self):
        criterias = dict()
        # разбиваем отрезок на N частей
        alpha_range = np.linspace(start=0, stop=1, num=20)[1:-1]
        for alpha in alpha_range:
            self.central_alpha = alpha
            # пересчитываем значение j
            j = self.linear_convolution_of_criteria()
            criterias[j] = alpha

        # находим альфу,которой соответствует минимум функции J и делаем ее текущей
        minimal_by_j_criteria = sorted(criterias.keys())[0]
        self.central_alpha = criterias[minimal_by_j_criteria]

    def optimize_lambdas(self):
        def dist():
            omega = noisy_criteria(self.filtered_signal(self.noisy_signal()))
            delta = closeness_criteria(self.filtered_signal(self.noisy_signal()), self.clear_signal())
            return ideal_point_distance_metric(omega, delta)

        criterias = dict()
        lambdas_range = [1 / l for l in range(1, 20)]
        for lambda_value in lambdas_range:
            self.lambda_weight = lambda_value
            # пересчитываем значение дистанции до опт. точки
            d = dist()
            criterias[d] = lambda_value

        # находим лямбду
        minimal_by_distance = sorted(criterias.keys())[0]
        self.lambda_weight = criterias[minimal_by_distance]


def main():
    s = Signal(SIGNAL_FUNC)

    # plot original
    plt.subplot(411)
    plt.title('Clear signal', fontdict={'size': 10})
    plt.plot(s.xs, s.clear_signal(), c='orange')

    # plot noisy
    plt.subplot(412)
    plt.title('Noisy signal', fontdict={'size': 10})
    plt.plot(s.xs, s.noisy_signal(), c='red')

    # plot naive filtered
    plt.subplot(413)
    plt.title('Naive filtered signal', fontdict={'size': 10})
    plt.plot(s.xs, s.filtered_signal(signal=s.noisy_signal()), c='green')

    # find optimal alphas, lambdas
    s.optimize_alphas()
    s.optimize_lambdas()

    # plot with optimal alphas
    plt.subplot(414)
    plt.title('Optimally filtered signal', fontdict={'size': 10})
    plt.plot(s.xs, s.filtered_signal(signal=s.noisy_signal()), c='blue')

    plt.tight_layout()
    plt.savefig('plot.png')


if __name__ == '__main__':
    main()
