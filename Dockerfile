FROM ghcr.io/prefix-dev/pixi:latest

WORKDIR /app

# Copy project files
COPY pyproject.toml pixi.lock ./
COPY src/ src/
COPY data/ data/
COPY tests/ tests/
COPY README.md LICENSE ./

# Install dependencies (default environment)
RUN pixi install

# Verify installation
RUN pixi run atbfetcher --help

ENTRYPOINT ["pixi", "run"]
CMD ["atbfetcher", "--help"]
