"""
Command-line interface for the MCTS materials framework.

`app` is resolved lazily. Importing it eagerly would mean that
`from mcts_framework.cli.builders import build_mcts` - assembling a search in
a Python script, with no CLI involved - pulls in typer through this module,
and fails where typer is absent. The builder and the results writer need only
the core dependencies, so they stay importable on their own.

© 2026. Triad National Security, LLC. All rights reserved.
"""

from typing import Any

from .builders import build_mcts
from .results import save_results

__all__ = ["app", "build_mcts", "save_results"]


def __getattr__(name: str) -> Any:
    """Resolve `app` on first access, so typer is imported only if it is used."""
    if name == "app":
        from .main import app

        return app
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
