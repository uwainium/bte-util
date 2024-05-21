import json

import matplotlib.path as Path
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

from projections.GeographicProjection import GeographicProjection

map_scale = 7318261.522857145
radius = 16 * 32 * 4
region_conversion = 16 * 32


def is_ccw(p):
    v = p.vertices - p.vertices[0, :]
    a = np.arctan2(v[1:, 1], v[1:, 0])
    return (a[1:] >= a[:-1]).astype(int).mean() >= 0.5


class Map:
    def __init__(self, geo_file: str, projection: GeographicProjection):
        self.features = []
        data = {}
        with open(geo_file, "r") as f:
            data = json.load(f)

        if data["features"]:
            for feature in data["features"]:
                self.features.append(Feature(feature, projection))

    def get_geo(self, filter_list):
        return [poly for geo in [feature.get_geo(filter_list) for feature in self.features] for poly in geo]

    def is_point_inside(self, filter_list, points):
        contains = [False] * len(points)
        for feature in self.features:
            contains = [a or b for a, b in zip(contains, feature.are_points_inside(filter_list, points))]

        return contains


class Feature:
    def __init__(self, data, projection: GeographicProjection):
        self.name = ""
        self.name = data["properties"].get("name", self.name)
        self.name = data["properties"].get("NAME", self.name)
        self.name = data["properties"].get("ADMIN", self.name)

        self.polygons = []

        if data["geometry"]["type"] == "Polygon":
            self.polygons.append(MyPolygon(data["geometry"]["coordinates"][0], projection))

        if data["geometry"]["type"] == "MultiPolygon":
            for poly in data["geometry"]["coordinates"]:
                self.polygons.append(MyPolygon(poly[0], projection))

    def get_geo(self, filter_list):
        if self.name in filter_list or not filter_list:
            return [polygon.path.vertices for polygon in self.polygons]
        else:
            return []

    def are_points_inside(self, filter_list, points):
        contains = [False] * len(points)
        if self.name in filter_list:
            for poly in self.polygons:
                contains = [a or b for a, b in
                            zip(contains, poly.path.contains_points(points, radius=radius * poly.radius_multiplier))]

        return contains


class MyPolygon:
    def __init__(self, data, projection: GeographicProjection):
        self.radius_multiplier = 1

        geometry = data

        geometry = [list(projection.from_geo(coord[0], coord[1])) for coord in geometry]
        geometry = [[coord[0] * map_scale, coord[1] * map_scale] for coord in geometry]

        self.path = Path.Path(geometry)

        if not is_ccw(self.path):
            self.radius_multiplier = -1


class ProjectionToMap:

    def __init__(self, regions, projection: GeographicProjection, geo_file, filtered_region_list=[]):
        self.regions = regions
        self.projection = projection
        self.border_map = Map(geo_file, self.projection)

        filtered_regions = list(
            filter(lambda x: x["region_name"] in filtered_region_list or not filtered_region_list, self.regions))

        self.filtered_states = [filtered_state for filtered_region in filtered_regions for filtered_state in
                                filtered_region["states"]]

        self.points_list = self.get_points(filtered_regions)

    def get_points(self, filtered_regions):
        to_ret = []
        for i, region in enumerate(filtered_regions):
            points = []
            with open("data/" + region["region_name"] + "_chunks", "r") as infile:
                for line in infile:
                    nums = line.split(' ')
                    x = int(nums[0])
                    y = -int(nums[1])

                    points.append([x * region_conversion, y * region_conversion])

            points_inside = self.border_map.is_point_inside(region["states"], points)

            to_ret.append({'region': i, 'points': points, 'points_inside': points_inside})

        return to_ret

    def get_points_inside(self):
        # to_ret = [(i, inside_point) for (i,inside_point) in enumerate(inside_points) for i[(_region["points"], _region["points_inside"]) for _region in self.points_list]]
        to_ret = []
        for _region in self.points_list:
            points = _region["points"]
            inside_points = _region["points_inside"]

            points_inside = []

            for i, inside_point in enumerate(inside_points):
                if inside_point:
                    points_inside.append(points[i])

            to_ret.append({'region': _region["region"], 'points': points_inside})

        return to_ret

    def print(self):
        fig, axs = plt.subplots()
        axs.set_aspect('equal', 'datalim')

        border_map_poly_collection = PolyCollection(self.border_map.get_geo(self.filtered_states), facecolors='None')
        axs.add_collection(border_map_poly_collection)

        polys = []
        colors = []

        for _region in self.points_list:
            region = self.regions[_region["region"]]

            for point in _region["points"]:
                polys.append([(point[0], point[1]),
                              (point[0], point[1] + region_conversion),
                              (point[0] + region_conversion, point[1] + region_conversion),
                              (point[0] + region_conversion, point[1])])

                colors.append(region["color"])

        poly_collection = PolyCollection(polys, facecolors=colors)
        axs.add_collection(poly_collection)
        axs.autoscale_view()
        plt.show()
