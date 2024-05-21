class OutOfProjectionBoundsException(Exception):

    @staticmethod
    def get():
        return OutOfProjectionBoundsException()

    @staticmethod
    def check_in_range(x: float, y: float, max_x: float, max_y: float) -> None:
        if abs(x) > (max_x + 0.1) or abs(y) > (max_y + 0.1):
            print(x, y, max_x, max_y)
            raise OutOfProjectionBoundsException.get()

    @staticmethod
    def check_longitude_latitude_in_range(longitude: float, latitude: float) -> None:
        OutOfProjectionBoundsException.check_in_range(longitude, latitude, 180, 90)

    def __init__(self):
        super()
