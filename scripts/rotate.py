import os, sys
from pathlib import Path
import numpy as np

# import DoZen
# add parent directory to path
try:
    # works when called externally via panel serve
    script_path = os.path.abspath(os.path.join(__file__,'..'))
    module_path = os.path.abspath(os.path.join(__file__,'..','..'))
except NameError:
    # works in notebook
    script_path = str(Path().resolve())
    module_path = str(Path().resolve().parent)
if module_path not in sys.path:
    sys.path.append(module_path)
from dozen import rotate

rotate.main(sys.argv)

