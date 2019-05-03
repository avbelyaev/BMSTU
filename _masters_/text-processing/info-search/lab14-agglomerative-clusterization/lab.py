import matplotlib.pyplot as plot
import numpy as np
import scipy.cluster.hierarchy as sch


POINTS = [
    (0.6, 1.9),
    (1.8, 1.6),
    (2.7, 2.0),
    (3.0, 2.1),
    (3.0, 2.6),
    (3.1, 4.5),
    (3.8, 0.6),
    (4.2, 2.7)
]


class Cluster:
    def __init__(self):
        self.points = []

    def add_point(self, point: tuple):
        self.points.append(point)

    def coordinates(self):
        return list(point[0] for point in self.points), list(point[1] for point in self.points)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f'Cluster {self.points}'


# fck data science
def single_link(distance_function):
    linkage = sch.linkage(distance_function, method='single')
    return sch.fcluster(linkage, 0.05 * distance_function.max(), 'distance')


# fck data science
def complete_link(distance_function):
    linkage = sch.linkage(distance_function, method='complete')
    return sch.fcluster(linkage, 5000 * distance_function.min(), 'distance')


def main():
    distance_function = sch.distance.pdist(np.array(POINTS), metric='cosine')
    # clusterized = single_link(distance_function)
    clusterized = complete_link(distance_function)

    # (-1) for each index since numeration after `fclutser` starts from 1
    clusterized = [cluster_indx - 1 for cluster_indx in clusterized]

    clusters_num = max(clusterized) + 1
    clusters = [Cluster() for _ in range(clusters_num)]

    i = 0
    for cluster_index in clusterized:
        clusters[cluster_index].add_point(POINTS[i])
        i += 1

    for c in clusters:
        xs, ys = c.coordinates()
        plot.scatter(xs, ys)
    plot.show()


if __name__ == '__main__':
    main()
