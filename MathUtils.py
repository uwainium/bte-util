import math


class MathUtils:
    ROOT3 = math.sqrt(3)

    @staticmethod
    def geo_to_spherical(geo: tuple[float, float]) -> tuple[float, float]:
        l = math.radians(geo[0])
        p = math.radians(90 - geo[1])
        return l, p

    @staticmethod
    def spherical_to_geo(spherical: tuple[float, float]) -> tuple[float, float]:
        lon = math.degrees(spherical[0])
        lat = 90 - math.degrees(spherical[1])
        return lon, lat

    @staticmethod
    def spherical_to_cartesian(spherical: tuple[float, float]) -> list[float, float, float]:
        sin_phi = math.sin(spherical[1])
        x = sin_phi * math.cos(spherical[0])
        y = sin_phi * math.sin(spherical[0])
        z = math.cos(spherical[1])
        return [x, y, z]

    @staticmethod
    def cartesian_to_spherical(cartesian: tuple[float, float, float]) -> tuple[float, float]:
        l = math.atan2(cartesian[1], cartesian[0])
        p = math.atan2(math.sqrt(cartesian[0] * cartesian[0] + cartesian[1] * cartesian[1]), cartesian[2])
        return l, p

    @staticmethod
    def produce_zyz_rotation_matrix(a: float, b: float, c: float) -> tuple[
        tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]:
        sina = math.sin(a)
        cosa = math.cos(a)
        sinb = math.sin(b)
        cosb = math.cos(b)
        sinc = math.sin(c)
        cosc = math.cos(c)

        mat = ((cosa * cosb * cosc - sinc * sina, -sina * cosb * cosc - sinc * cosa, cosc * sinb),
               (sinc * cosb * cosa + cosc * sina, cosc * cosa - sinc * cosb * sina, sinc * sinb),
               (-sinb * cosa, sinb * sina, cosb))

        return mat

    @staticmethod
    def mat_vec_prod_d(
            matrix: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]],
            vector: list[float]) -> tuple[float, float, float]:
        result = [0.0] * len(vector)

        for i in range(len(result)):
            for j in range(len(matrix[i])):
                result[i] += matrix[i][j] * vector[j]

        return result[0], result[1], result[2]
