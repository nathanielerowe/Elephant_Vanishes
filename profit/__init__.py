import _profit

try:
  import _profit
except:
  raise Exception("Need to build the package first!")

__all__ = ['profit']

from .lib import *
