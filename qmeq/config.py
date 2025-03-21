"""Module containing global configuration settings."""

# Global debug flag to control debug prints
DEBUG     = False
verbosity = 0

def debug_print(*args, **kwargs):
    """Wrapper for print that only prints if DEBUG is True."""
    if DEBUG:
        print(*args, **kwargs)

def verb_print(*args, **kwargs):
    if verbosity > 0:
        print(*args, **kwargs)

def verb_print_(level, *args, **kwargs):
    #print("Verbosity level:", verbosity, "Level:", level)
    if verbosity >= level:
        print(*args, **kwargs)
