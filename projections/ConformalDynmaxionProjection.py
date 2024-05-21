import math
import struct

from projections.DymaxionProjection import DymaxionProjection
from MathUtils import MathUtils


class ConformalDynmaxionProjection(DymaxionProjection):
    VECTOR_SCALE_FACTOR = 1.0 / 1.1473979730192934
    SIDE_LENGTH = 256

    vx = []
    vy = []

    for i in range(SIDE_LENGTH + 1):
        vx.append([0] * (SIDE_LENGTH + 1 - i))
        vy.append([0] * (SIDE_LENGTH + 1 - i))

    f = open("./data/conformal", "rb")

    for v in range(SIDE_LENGTH + 1):
        for u in range(SIDE_LENGTH + 1 - v):
            vx[u][v] = struct.unpack('>d', f.read(8))[0] * VECTOR_SCALE_FACTOR
            vy[u][v] = struct.unpack('>d', f.read(8))[0] * VECTOR_SCALE_FACTOR

    f.close()

    def triangle_transform(self, vec: tuple[float, float, float]) -> tuple[float, float]:
        c = list(super().triangle_transform(vec))

        x = c[0]
        y = c[1]

        c[0] /= ConformalDynmaxionProjection.ARC
        c[1] /= ConformalDynmaxionProjection.ARC

        c[0] += 0.5
        c[1] += MathUtils.ROOT3 / 6

        c = list(self.apply_newtons_method(x, y, c[0], c[1], 5))

        c[0] -= 0.5
        c[1] -= MathUtils.ROOT3 / 6

        c[0] *= ConformalDynmaxionProjection.ARC
        c[1] *= ConformalDynmaxionProjection.ARC

        return c[0], c[1]

    def inverse_triangle_transform(self, x: float, y: float) -> tuple[float, float, float]:

        x /= ConformalDynmaxionProjection.ARC
        y /= ConformalDynmaxionProjection.ARC

        x += 0.5
        y += MathUtils.ROOT3 / 6

        # some shit
        c = self.get_interpolated_vector(x, y)

        return super().inverse_triangle_transform(c[0], c[1])

    def meters_per_unit(self) -> float:
        return (40075017.0 / (2.0 * math.pi)) / ConformalDynmaxionProjection.VECTOR_SCALE_FACTOR

    def __str__(self) -> str:
        return "Conformal Dymaxion"

    def get_interpolated_vector(self, x: float, y: float) -> tuple[float, float, float, float, float, float]:
        x *= ConformalDynmaxionProjection.SIDE_LENGTH
        y *= ConformalDynmaxionProjection.SIDE_LENGTH

        v = 2 * y / MathUtils.ROOT3
        u = x - v * 0.5

        u1 = int(u)
        v1 = int(v)

        if u1 < 0:
            u1 = 0
        elif u1 >= ConformalDynmaxionProjection.SIDE_LENGTH:
            u1 = ConformalDynmaxionProjection.SIDE_LENGTH - 1

        if v1 < 0:
            v1 = 0
        elif v1 >= ConformalDynmaxionProjection.SIDE_LENGTH - u1:
            v1 = ConformalDynmaxionProjection.SIDE_LENGTH - u1 - 1

        valx1 = 0
        valy1 = 0
        valx2 = 0
        valy2 = 0
        valx3 = 0
        valy3 = 0
        y3 = 0
        x3 = 0

        flip = 1

        if (y < -MathUtils.ROOT3 * (x - u1 - v1 - 1)) or (v1 == ConformalDynmaxionProjection.SIDE_LENGTH - u1 - 1):
            valx1 = self.vx[u1][v1]
            valy1 = self.vy[u1][v1]
            valx2 = self.vx[u1][v1 + 1]
            valy2 = self.vy[u1][v1 + 1]
            valx3 = self.vx[u1 + 1][v1]
            valy3 = self.vy[u1 + 1][v1]

            y3 = 0.5 * MathUtils.ROOT3 * v1
            x3 = (u1 + 1) + 0.5 * v1
        else:
            valx1 = self.vx[u1][v1 + 1]
            valy1 = self.vy[u1][v1 + 1]
            valx2 = self.vx[u1 + 1][v1]
            valy2 = self.vy[u1 + 1][v1]
            valx3 = self.vx[u1 + 1][v1 + 1]
            valy3 = self.vy[u1 + 1][v1 + 1]

            flip = -1
            y = -y

            y3 = -(0.5 * MathUtils.ROOT3 * (v1 + 1))
            x3 = (u1 + 1) + 0.5 * (v1 + 1)

        w1 = -(y - y3) / MathUtils.ROOT3 - (x - x3)
        w2 = 2 * (y - y3) / MathUtils.ROOT3
        w3 = 1 - w1 - w2

        return (valx1 * w1 + valx2 * w2 + valx3 * w3,
                valy1 * w1 + valy2 * w2 + valy3 * w3,
                (valx3 - valx1) * ConformalDynmaxionProjection.SIDE_LENGTH,
                ConformalDynmaxionProjection.SIDE_LENGTH * flip * (2 * valx2 - valx1 - valx3) / MathUtils.ROOT3,
                (valy3 - valy1) * ConformalDynmaxionProjection.SIDE_LENGTH,
                ConformalDynmaxionProjection.SIDE_LENGTH * flip * (2 * valy2 - valy1 - valy3) / MathUtils.ROOT3)

    def apply_newtons_method(self, expected_f: float, expected_g: float, x_est: float, y_est: float, _iter: int):
        for i in range(_iter):
            c = self.get_interpolated_vector(x_est, y_est)

            f = c[0] - expected_f
            g = c[1] - expected_g

            dfdx = c[2]
            dfdy = c[3]
            dgdx = c[4]
            dgdy = c[5]

            determinant = 1 / (dfdx * dgdy - dfdy * dgdx)

            x_est -= determinant * (dgdy * f - dfdy * g)
            y_est -= determinant * (-dgdx * f + dfdx * g)

        return x_est, y_est
