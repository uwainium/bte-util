import math
import sys

from projections.GeographicProjection import GeographicProjection
from MathUtils import MathUtils
from exceptions.OutOfProjectionBoundsException import OutOfProjectionBoundsException


class DymaxionProjection(GeographicProjection):
    ARC = 2 * math.asin(math.sqrt(5 - math.sqrt(5)) / math.sqrt(10))
    Z = math.sqrt(5 + 2 * math.sqrt(5)) / math.sqrt(15)
    EL = math.sqrt(8) / math.sqrt(5 + math.sqrt(5))
    EL6 = EL / 6
    DVE = math.sqrt(3 + math.sqrt(5)) / math.sqrt(5 + math.sqrt(5))
    R = -3 * EL6 / DVE

    NEWTON = 5

    VERTICES = [
        (10.536199, 64.700000),
        (-5.245390, 2.300882),
        (58.157706, 10.447378),
        (122.300000, 39.100000),
        (-143.478490, 50.103201),
        (-67.132330, 23.717925),
        (36.521510, -50.103200),
        (112.867673, -23.717930),
        (174.754610, -2.300882),
        (-121.842290, -10.447350),
        (-57.700000, -39.100000),
        (-169.463800, -64.700000),
    ]

    ISO = [
        (2, 1, 6),
        (1, 0, 2),
        (0, 1, 5),
        (1, 5, 10),
        (1, 6, 10),
        (7, 2, 6),
        (2, 3, 7),
        (3, 0, 2),
        (0, 3, 4),
        (4, 0, 5),
        (5, 4, 9),
        (9, 5, 10),
        (10, 9, 11),
        (11, 6, 10),
        (6, 7, 11),
        (8, 3, 7),
        (8, 3, 4),
        (8, 4, 9),
        (9, 8, 11),
        (7, 8, 11),
        (11, 6, 7),
        (3, 7, 8)
    ]

    CENTER_MAP = [
        [-3, 7],
        [-2, 5],
        [-1, 7],
        [2, 5],
        [4, 5],
        [-4, 1],
        [-3, -1],
        [-2, 1],
        [-1, -1],
        [0, 1],
        [1, -1],
        [2, 1],
        [3, -1],
        [4, 1],
        [5, -1],
        [-3, -5],
        [-1, -5],
        [1, -5],
        [2, -7],
        [-4, -7],
        [-5, -5],
        [-2, -7]
    ]

    FLIP_TRIANGLE = [
        True, False, True, False, False,
        True, False, True, False, True, False, True, False, True, False,
        True, True, True, False, False,
        True, False
    ]

    CENTROIDS: list[tuple[float, float, float]] = [None] * 22

    ROTATION_MATRICES: list[
        tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]] = [None] * 22

    INVERSE_ROTATION_MATRICES: list[
        tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]] = [None] * 22

    FACE_ON_GRID = [
        -1, -1, 0, 1, 2, -1, -1, 3, -1, 4, -1,
        -1, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
        20, 19, 15, 21, 16, -1, 17, 18, -1, -1, -1,
    ]

    # begin static math shit

    @staticmethod
    def find_triangle_grid(x: float, y: float) -> int:
        xp = x / DymaxionProjection.ARC
        yp = y / (DymaxionProjection.ARC * MathUtils.ROOT3)

        row = None
        if yp > -0.25:
            if yp < 0.25:
                row = 1
            elif yp <= 0.75:
                row = 0
                yp = 0.5 - yp
            else:
                return -1
        elif yp >= -0.75:
            row = 2
            yp = -yp - 0.5
        else:
            return -1

        yp += 0.25

        xr = xp - yp
        yr = xp + yp

        gx = math.floor(xr)
        gy = math.floor(yr)

        col = 2 * gx + (1 if gy != gx else 0) + 6

        if col < 0 or col >= 11:
            return -1

        return DymaxionProjection.FACE_ON_GRID[row * 11 + col]

    @staticmethod
    def y_rot(spherical, rot: float) -> tuple[float, float]:
        c = MathUtils.spherical_to_cartesian(spherical)

        x = c[0]
        c[0] = c[2] * math.sin(rot) + x * math.cos(rot)
        c[2] = c[2] * math.cos(rot) - x * math.sin(rot)

        mag = math.sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2])
        c[0] /= mag
        c[1] /= mag
        c[2] /= mag

        return math.atan2(c[1], c[0]), math.atan2(math.sqrt(c[0] * c[0] + c[1] * c[1]), c[2])

    vertices_cartesian = [[]] * len(VERTICES)

    for i in range(22):
        CENTER_MAP[i][0] *= 0.5 * ARC
        CENTER_MAP[i][1] *= ARC * MathUtils.ROOT3 / 12

    for i in range(len(VERTICES)):
        vertex_spherical = MathUtils.geo_to_spherical(VERTICES[i])
        vertex = MathUtils.spherical_to_cartesian(vertex_spherical)
        vertices_cartesian[i] = vertex
        VERTICES[i] = vertex_spherical

    for i in range(22):
        vec1 = vertices_cartesian[ISO[i][0]]
        vec2 = vertices_cartesian[ISO[i][1]]
        vec3 = vertices_cartesian[ISO[i][2]]

        xsum = vec1[0] + vec2[0] + vec3[0]
        ysum = vec1[1] + vec2[1] + vec3[1]
        zsum = vec1[2] + vec2[2] + vec3[2]
        mag = math.sqrt(xsum * xsum + ysum * ysum + zsum * zsum)
        CENTROIDS[i] = (xsum / mag, ysum / mag, zsum / mag)

        centroid_spherical = MathUtils.cartesian_to_spherical(CENTROIDS[i])
        centroid_lambda = centroid_spherical[0]
        centroid_phi = centroid_spherical[1]

        vertex = VERTICES[ISO[i][0]]
        v = (vertex[0] - centroid_lambda, vertex[1])
        v = y_rot(v, -centroid_phi)

        ROTATION_MATRICES[i] = MathUtils.produce_zyz_rotation_matrix(-centroid_lambda, -centroid_phi,
                                                                     (math.pi / 2) - v[0])
        INVERSE_ROTATION_MATRICES[i] = MathUtils.produce_zyz_rotation_matrix(v[0] - (math.pi / 2), centroid_phi,
                                                                             centroid_lambda)

    def find_triangle(self, vector: list[float]) -> int:
        _min = sys.float_info.max
        face = 0

        for i in range(20):
            xd = DymaxionProjection.CENTROIDS[i][0] - vector[0]
            yd = DymaxionProjection.CENTROIDS[i][1] - vector[1]
            zd = DymaxionProjection.CENTROIDS[i][2] - vector[2]

            dissq = xd * xd + yd * yd + zd * zd

            if dissq < _min:
                if dissq < 0.1:
                    return i

                face = i
                _min = dissq

        return face

    def triangle_transform(self, vec: tuple[float, float, float]) -> tuple[float, float]:
        s = DymaxionProjection.Z / vec[2]

        xp = s * vec[0]
        yp = s * vec[1]

        a = math.atan((2 * yp / MathUtils.ROOT3 - DymaxionProjection.EL6) / DymaxionProjection.DVE)
        b = math.atan((xp - yp / MathUtils.ROOT3 - DymaxionProjection.EL6) / DymaxionProjection.DVE)
        c = math.atan((-xp - yp / MathUtils.ROOT3 - DymaxionProjection.EL6) / DymaxionProjection.DVE)

        return 0.5 * (b - c), (2 * a - b - c) / (2 * MathUtils.ROOT3)

    def inverse_triangle_transform_newton(self, xpp: float, ypp: float) -> tuple[float, float, float]:
        tanaoff = math.tan(MathUtils.ROOT3 * ypp + xpp)
        tanboff = math.tan(2 * xpp)

        anumer = tanaoff * tanaoff + 1
        bnumer = tanboff * tanboff + 1

        tana = tanaoff
        tanb = tanboff
        tanc = 0

        adenom = 1
        bdenom = 1

        for i in range(DymaxionProjection.NEWTON):
            f = tana + tanb + tanc - DymaxionProjection.R
            fp = anumer * adenom * adenom + bnumer * bdenom * bdenom + 1

            tanc -= f / fp

            adenom = 1 / (1 - tanc * tanaoff)
            bdenom = 1 / (1 - tanc * tanboff)

            tana = (tanc + tanaoff) * adenom
            tanb = (tanc + tanboff) * bdenom

        yp = MathUtils.ROOT3 * (DymaxionProjection.DVE * tana + DymaxionProjection.EL6) / 2
        xp = DymaxionProjection.DVE * tanb + yp / MathUtils.ROOT3 + DymaxionProjection.EL6

        xpoz = xp / DymaxionProjection.Z
        ypoz = yp / DymaxionProjection.Z

        z = 1 / math.sqrt(1 + xpoz * xpoz + ypoz * ypoz)

        return z * xpoz, z * ypoz, z

    def inverse_triangle_transform(self, x: float, y: float) -> tuple[float, float, float]:
        return self.inverse_triangle_transform_newton(x, y)

    def from_geo(self, longitude: float, latitude: float) -> tuple[float, float]:
        OutOfProjectionBoundsException.check_longitude_latitude_in_range(longitude, latitude)

        vector = MathUtils.spherical_to_cartesian((MathUtils.geo_to_spherical((longitude, latitude))))

        face = self.find_triangle(vector)

        pvec = MathUtils.mat_vec_prod_d(DymaxionProjection.ROTATION_MATRICES[face], vector)
        projected_vec = list(self.triangle_transform(pvec))

        if DymaxionProjection.FLIP_TRIANGLE[face]:
            projected_vec[0] = -projected_vec[0]
            projected_vec[1] = -projected_vec[1]

        vector[0] = projected_vec[0]

        if ((face == 15 and vector[0] > projected_vec[1] * MathUtils.ROOT3) or face == 14) and vector[0] > 0:
            projected_vec[0] = 0.5 * vector[0] - 0.5 * MathUtils.ROOT3 * projected_vec[1]
            projected_vec[1] = 0.5 * MathUtils.ROOT3 * vector[0] + 0.5 * projected_vec[1]
            face += 6

        projected_vec[0] += DymaxionProjection.CENTER_MAP[face][0]
        projected_vec[1] += DymaxionProjection.CENTER_MAP[face][1]

        return projected_vec[0], projected_vec[1]

    def to_geo(self, x: float, y: float) -> tuple[float, float]:
        face = self.find_triangle_grid(x, y)

        if face == -1:
            raise OutOfProjectionBoundsException.get()

        x -= DymaxionProjection.CENTER_MAP[face][0]
        y -= DymaxionProjection.CENTER_MAP[face][1]

        match face:
            case 14:
                if x > 0:
                    raise OutOfProjectionBoundsException.get()
            case 20:
                if -y * MathUtils.ROOT3 > x:
                    raise OutOfProjectionBoundsException.get()
            case 15:
                if x > 0 and x > y * MathUtils.ROOT3:
                    raise OutOfProjectionBoundsException.get()
            case 21:
                if x < 0 or -y * MathUtils.ROOT3 > x:
                    raise OutOfProjectionBoundsException.get()

        if DymaxionProjection.FLIP_TRIANGLE[face]:
            x = -x
            y = -y

        c = self.inverse_triangle_transform(x, y)
        x = c[0]
        y = c[1]
        z = c[2]

        vecp = MathUtils.mat_vec_prod_d(DymaxionProjection.ROTATION_MATRICES[face], [x, y, z])

        return MathUtils.spherical_to_geo(MathUtils.cartesian_to_spherical(vecp))

    def bounds(self) -> tuple[float, float, float, float]:
        return (
            -3 * DymaxionProjection.ARC, -0.75 * DymaxionProjection.ARC * MathUtils.ROOT3, 2.5 * DymaxionProjection.ARC,
            0.75 * DymaxionProjection.ARC * MathUtils.ROOT3)

    def upright(self) -> bool:
        return False

    def meters_per_unit(self) -> float:
        return math.sqrt(
            510100000000000.0 / (20 * MathUtils.ROOT3 * DymaxionProjection.ARC * DymaxionProjection.ARC / 4))

    def __str__(self) -> str:
        return "Dymaxion"
