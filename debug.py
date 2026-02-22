import math

def calculate_intersection(height: float, width: float, angle: float):
    # Middle of the structure
    x_middle = width / 2

    # Slope of the line
    slope = math.tan(math.radians(angle))

    y_intersection = slope * x_middle

    return x_middle, y_intersection

print(calculate_intersection(100, 200, 10))