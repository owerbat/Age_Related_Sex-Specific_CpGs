import numpy as np
from shapely.geometry.polygon import Polygon


def get_polygon(betas, ages, slope, intercept):
    sigma = np.std(betas)
    min_age = np.min(ages)
    max_age = np.max(ages)

    point1 = (min_age, slope*min_age+intercept-3*sigma)
    point2 = (min_age, slope*min_age+intercept+3*sigma)
    point3 = (max_age, slope*max_age+intercept+3*sigma)
    point4 = (max_age, slope*max_age+intercept-3*sigma)

    return (point1, point2, point3, point4)


def get_polygons_areas(coordinates1, coordinates2):
    polygon1 = Polygon(coordinates1)
    polygon2 = Polygon(coordinates2)

    intersection = polygon1.intersection(polygon2)
    union = polygon1.union(polygon2)

    try:
        coordinates = intersection.exterior.xy
    except AttributeError:
        coordinates = [[0, 0, 0, 0], [0, 0, 0, 0]]

    return intersection.area, union.area, coordinates


def get_polygon_xy(coordinates):
    polygon = Polygon(coordinates)
    return polygon.exterior.xy
