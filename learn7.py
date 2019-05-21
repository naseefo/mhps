


# Import the os module, for the os.walk function
import os


def find(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)

print(find('r_se21', 'data\earthquakes'))