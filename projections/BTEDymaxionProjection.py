import math

from MathUtils import MathUtils
from exceptions.OutOfProjectionBoundsException import OutOfProjectionBoundsException
from projections.ConformalDynmaxionProjection import ConformalDynmaxionProjection


class BTEDymaxionProjection(ConformalDynmaxionProjection):
    THETA = math.radians(-150)
    SIN_THETA = math.sin(THETA)
    COS_THETA = math.cos(THETA)
    BERING_X = -0.3420420960118339
    BERING_Y = -0.322211064085279
    ARCTIC_Y = -0.2
    ARCTIC_M = ((ARCTIC_Y - MathUtils.ROOT3 * ConformalDynmaxionProjection.ARC / 4) /
                (BERING_X - -0.5 * ConformalDynmaxionProjection.ARC))
    ARCTIC_B = ARCTIC_Y - ARCTIC_M * BERING_X
    ALEUTIAN_Y = -0.5000446805492526
    ALEUTIAN_XL = -0.5149231279757507
    ALEUTIAN_XR = -0.45
    ALEUTIAN_M = (BERING_Y - ALEUTIAN_Y) / (BERING_X - ALEUTIAN_XR)
    ALEUTIAN_B = BERING_Y - ALEUTIAN_M * BERING_X

    def from_geo(self, longitude: float, latitude: float) -> tuple[float, float]:
        c = list(super().from_geo(longitude, latitude))
        x = c[0]
        y = c[1]

        easia = self.is_eurasian_part(x, y)

        y -= 0.75 * BTEDymaxionProjection.ARC * MathUtils.ROOT3

        if easia:
            x += BTEDymaxionProjection.ARC

            t = x
            x = BTEDymaxionProjection.COS_THETA * x - BTEDymaxionProjection.SIN_THETA * y
            y = BTEDymaxionProjection.SIN_THETA * t + BTEDymaxionProjection.COS_THETA * y

        else:
            x -= BTEDymaxionProjection.ARC

        c[0] = y
        c[1] = -x

        return c[0], c[1]

    def to_geo(self, x: float, y: float) -> tuple[float, float]:
        easia = False
        if y < 0:
            easia = x > 0
        elif y > BTEDymaxionProjection.ARC / 2:
            easia = x > -MathUtils.ROOT3 * BTEDymaxionProjection.ARC / 2
        else:
            easia = y * -MathUtils.ROOT3 < x

        t = x
        x = -y
        y = t

        if easia:
            t = x
            x = BTEDymaxionProjection.COS_THETA * x + BTEDymaxionProjection.SIN_THETA * y
            y = BTEDymaxionProjection.COS_THETA * y - BTEDymaxionProjection.SIN_THETA * t
            x -= BTEDymaxionProjection.ARC
        else:
            x += BTEDymaxionProjection.ARC

        y += 0.75 * BTEDymaxionProjection.ARC * MathUtils.ROOT3

        if easia != self.is_eurasian_part(x, y):
            raise OutOfProjectionBoundsException.get()
        
        return super().to_geo(x, y)

    def is_eurasian_part(self, x: float, y: float) -> bool:
        if x > 0:
            return False

        if x < -0.5 * BTEDymaxionProjection.ARC:
            return True

        if y > MathUtils.ROOT3 * BTEDymaxionProjection.ARC / 4:
            return x < 0

        if y < BTEDymaxionProjection.ALEUTIAN_Y:
            return y < (BTEDymaxionProjection.ALEUTIAN_Y + BTEDymaxionProjection.ALEUTIAN_XL) - x

        if y > BTEDymaxionProjection.BERING_Y:

            if y < BTEDymaxionProjection.ARCTIC_Y:
                return x < BTEDymaxionProjection.BERING_X

            return y < BTEDymaxionProjection.ARCTIC_M * x + BTEDymaxionProjection.ARCTIC_B

        return y > BTEDymaxionProjection.ALEUTIAN_M * x + BTEDymaxionProjection.ALEUTIAN_B

    def bounds(self) -> tuple[float, float, float, float]:
        return (-1.5 * BTEDymaxionProjection.ARC * MathUtils.ROOT3,
                -1.5 * BTEDymaxionProjection.ARC,
                3 * BTEDymaxionProjection.ARC,
                MathUtils.ROOT3 * BTEDymaxionProjection.ARC)

    def __str__(self) -> str:
        return "BuildTheEarth Conformal Dymaxion"
