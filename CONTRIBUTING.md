# Contributing to atbfetcher

Thanks for your interest in contributing to atbfetcher! This document explains how to get started.

## Development Setup

1. Install [Pixi](https://pixi.sh):
   ```bash
   curl -fsSL https://pixi.sh/install.sh | bash
   ```

2. Clone and install:
   ```bash
   git clone https://github.com/happykhan/atbfetcher.git
   cd atbfetcher
   pixi install
   ```

3. Install the test and lint environments:
   ```bash
   pixi install -e test
   pixi install -e lint
   ```

## Running Tests

```bash
pixi run -e test test
```

## Linting and Formatting

We use [Ruff](https://docs.astral.sh/ruff/) for linting and formatting:

```bash
# Check for lint errors
pixi run -e lint lint

# Check formatting
pixi run -e lint format-check

# Auto-format
pixi run -e lint format
```

## Pull Request Process

1. Fork the repository and create a feature branch from `main`
2. Write tests for any new functionality
3. Ensure all tests pass and linting is clean
4. Submit a pull request with a clear description of the changes

## Code Style

- Follow the existing code patterns in the project
- Use type hints for function signatures
- Write NumPy-style docstrings for public functions
- Keep functions focused and modular

## Reporting Issues

Please use the [GitHub issue tracker](https://github.com/happykhan/atbfetcher/issues) to report bugs or request features. Include:

- Steps to reproduce (for bugs)
- Expected vs actual behaviour
- Python version and OS
