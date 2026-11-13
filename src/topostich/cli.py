"""Compatibility import; canonical implementation is topotrace.cli."""
from topotrace.cli import *  # noqa: F401,F403

if __name__ == "__main__":
    raise SystemExit(main())
