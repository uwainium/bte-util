class GeographicProjection:
    def to_geo(self, x: float, y: float) -> tuple[float, float]:
        raise NotImplementedError

    def from_geo(self, longitude: float, latitude: float) -> tuple[float, float]:
        raise NotImplementedError

    def meters_per_unit(self) -> float:
        raise NotImplementedError

    def bounds(self) -> tuple[float, float, float, float]:
        raise NotImplementedError

    def upright(self) -> bool:
        raise NotImplementedError

    def vector(self, x: float, y: float, north: float, east: float) -> tuple[float, float]:
        raise NotImplementedError

    def tissot(self, longitude: float, latitude: float) -> tuple[float, float, float, float]:
        raise NotImplementedError

    def azimuth(self, x: float, y: float, angle: float) -> float:
        raise NotImplementedError

    def properties(self) -> dict[str, object]:
        raise NotImplementedError
