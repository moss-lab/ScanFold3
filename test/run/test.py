import sys
# add lib to import path
sys.path.append("../../lib")

# from (name of folder) import (name of binary)
from fold import fold

example = fold.BasePair(0, 1, 'A', 'U', 0.1, 0.2, 0.3, 0.4, 2, 2)
example.print()
