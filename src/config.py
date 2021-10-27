# config.py
# stores variables to be imported across the project

import os

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data = root + '/data/'
raw = data + '/raw/'
processed = data + '/processed/'