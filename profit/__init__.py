try:
  import _profit
except:
  raise Exception("Need to build the package first!")

__all__ = ['profit']

# Re-export everything
from .lib import *
from .pylib import *
