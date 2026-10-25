"""Design-based inference for a fixed graph under independent node sampling."""
from ._version import __version__
from .core import estimate, ResourceLimitError
__all__ = ["estimate", "ResourceLimitError"]
