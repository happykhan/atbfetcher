"""User configuration for atbfetcher.

Reads and writes ``~/.atbfetcher/config.toml``.  Settings in the config file
serve as defaults that can always be overridden by explicit CLI flags.

Supported keys (all optional):

    [defaults]
    cache_dir = "~/.atbfetcher"   # default cache directory
    output_format = "auto"        # table|tsv|csv|json|auto
    threads = 4                   # download/decompress threads

Example usage::

    from atbfetcher.config import load_config, get_config_value

    cfg = load_config()
    cache_dir = get_config_value(cfg, "defaults", "cache_dir", fallback="~/.atbfetcher")
"""

from __future__ import annotations

import tomllib
from pathlib import Path
from typing import Any

CONFIG_DIR = Path.home() / ".atbfetcher"
CONFIG_PATH = CONFIG_DIR / "config.toml"

# Keys that are valid in each section, with their defaults
_DEFAULTS: dict[str, dict[str, Any]] = {
    "defaults": {
        "cache_dir": str(CONFIG_DIR),
        "output_format": "auto",
        "threads": None,
    }
}


def load_config(path: Path | None = None) -> dict[str, Any]:
    """Load the atbfetcher config file.

    Parameters
    ----------
    path : Path, optional
        Path to the TOML config file. Defaults to ``~/.atbfetcher/config.toml``.

    Returns
    -------
    dict
        Parsed config, or empty dict if file doesn't exist.
    """
    config_path = path or CONFIG_PATH
    if not config_path.exists():
        return {}
    try:
        with open(config_path, "rb") as f:
            return tomllib.load(f)
    except Exception:
        return {}


def get_config_value(
    config: dict[str, Any],
    section: str,
    key: str,
    fallback: Any = None,
) -> Any:
    """Get a value from the config, returning *fallback* if absent.

    Parameters
    ----------
    config : dict
        Loaded config dict (from :func:`load_config`).
    section : str
        Top-level section name (e.g. ``"defaults"``).
    key : str
        Key within the section.
    fallback : any
        Value to return when the key is absent.

    Returns
    -------
    any
        The config value, or *fallback*.
    """
    return config.get(section, {}).get(key, fallback)


def save_config(config: dict[str, Any], path: Path | None = None) -> None:
    """Write a config dict to disk as TOML.

    Parameters
    ----------
    config : dict
        Config to write.
    path : Path, optional
        Destination path. Defaults to ``~/.atbfetcher/config.toml``.
    """
    import tomli_w  # soft dependency, only needed for writing

    config_path = path or CONFIG_PATH
    config_path.parent.mkdir(parents=True, exist_ok=True)
    with open(config_path, "wb") as f:
        tomli_w.dump(config, f)


def _save_config_manual(config: dict[str, Any], path: Path | None = None) -> None:
    """Write a config dict to disk without tomli_w (manual TOML serialiser).

    Supports only string/int/float/bool/None values — sufficient for our
    simple config structure.
    """
    config_path = path or CONFIG_PATH
    config_path.parent.mkdir(parents=True, exist_ok=True)

    lines: list[str] = []
    for section, values in config.items():
        lines.append(f"[{section}]")
        for key, val in values.items():
            if val is None:
                continue
            if isinstance(val, bool):
                lines.append(f"{key} = {str(val).lower()}")
            elif isinstance(val, str):
                escaped = val.replace("\\", "\\\\").replace('"', '\\"')
                lines.append(f'{key} = "{escaped}"')
            else:
                lines.append(f"{key} = {val}")
        lines.append("")

    config_path.write_text("\n".join(lines))


def write_config(config: dict[str, Any], path: Path | None = None) -> None:
    """Write a config dict to disk as TOML (no extra dependencies).

    Parameters
    ----------
    config : dict
        Config dict to persist.
    path : Path, optional
        Destination path. Defaults to ``~/.atbfetcher/config.toml``.
    """
    _save_config_manual(config, path)


def default_config() -> dict[str, Any]:
    """Return a fresh default config dict."""
    return {
        "defaults": {
            "cache_dir": str(CONFIG_DIR),
            "output_format": "auto",
        }
    }
