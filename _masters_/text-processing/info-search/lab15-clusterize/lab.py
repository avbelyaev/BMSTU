import matplotlib.pyplot as plot
import math

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

INIT_CLUSTER_CENTERS = [
    (1.8, 1.6),
    (3.0, 2.6)
]


class Point:
    def __init__(self, x: float, y: float):
        self.x = x
        self.y = y
        self.cluster_index = None

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f'[{self.x} {self.y}]'


class Cluster:
    def __init__(self, center: Point):
        self.points = []
        self.center = center

    def add_point(self, p: Point):
        self.points.append(p)

    def remove_point(self, p: Point):
        if p in self.points:
            self.points.remove(p)

    def update_center(self):
        if 0 != len(self.points):
            x_avg = sum(list(map(lambda p: p.x, self.points))) / len(self.points)
            y_avg = sum(list(map(lambda p: p.y, self.points))) / len(self.points)
            self.center = Point(x_avg, y_avg)

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return f'Cluster {len(self.points)}: {self.center}'


def reshuffle(points: list, clusters: list) -> bool:
    def dist_between(a: Point, b: Point):
        x = math.fabs(a.x - b.x)
        y = math.fabs(a.y - b.y)
        return math.sqrt(x**2 + y**2)

    smth_has_changed = False
    for p in points:

        curr_point = p
        min_dist = float('inf')
        new_cluster_index = 0

        for j in range(len(clusters)):
            curr_center = clusters[j].center
            distance = dist_between(curr_point, curr_center)
            if distance <= min_dist:
                new_cluster_index = j
                min_dist = distance

        if curr_point.cluster_index != new_cluster_index:

            # remove point from old cluster
            if p.cluster_index is not None:
                clusters[curr_point.cluster_index].remove_point(curr_point)

            # append point to new cluster
            clusters[new_cluster_index].add_point(curr_point)
            curr_point.cluster_index = new_cluster_index
            smth_has_changed = True

    return smth_has_changed


def k_means_clusterize(points: list, initial_centers: list):
    clusters = [Cluster(initial_center) for initial_center in initial_centers]

    not_converged = True
    while not_converged:
        # assign points to nearest clusters
        not_converged = reshuffle(points, clusters)

        # update cluster centers
        [c.update_center() for c in clusters]
        print(clusters)

    return clusters


def main():
    points = list(map(lambda p: Point(p[0], p[1]), POINTS))
    centers = list(map(lambda p: Point(p[0], p[1]), INIT_CLUSTER_CENTERS))

    clusters = k_means_clusterize(points, centers)

    xs1, ys1 = list(point.x for point in clusters[0].points), list(point.y for point in clusters[0].points)
    plot.scatter(xs1, ys1)

    xs2, ys2 = list(point.x for point in clusters[1].points), list(point.y for point in clusters[1].points)
    plot.scatter(xs2, ys2)

    cxs, cys = list(c.center.x for c in clusters), list(c.center.y for c in clusters)
    plot.scatter(cxs, cys, color='black')

    plot.gca().legend(('C1', 'C2', 'centers'))
    plot.show()


if __name__ == '__main__':
    main()
