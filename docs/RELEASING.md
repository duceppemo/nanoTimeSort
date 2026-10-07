# Releasing a new version

The version lives in **one place**: `nanotimesort/__init__.py` (`__version__`).
`pyproject.toml` reads it dynamically.

```bash
# 1. Bump the version
$EDITOR nanotimesort/__init__.py        # e.g. __version__ = "1.1.0"
git commit -am "Bump version to 1.1.0"
git push

# 2. Tag and push — this alone publishes to PyPI within ~1 minute
git tag -a v1.1.0 -m "nanoTimeSort v1.1.0"
git push origin v1.1.0

# 3. (Optional) add release notes on GitHub
gh release create v1.1.0 --title "nanoTimeSort v1.1.0" --notes "..."
```

Notes:

- The publish workflow (`.github/workflows/publish.yml`) triggers on the tag push *and* on the
  GitHub release; duplicate uploads are harmless (`skip-existing: true`).
- The workflow fails fast if the tag does not match `__version__`, so a forgotten bump can't
  publish a wrong version.
- Publishing uses PyPI trusted publishing (OIDC) — no API token to rotate.
- Zenodo archives every GitHub release automatically and mints a new versioned DOI under the
  same concept DOI; the README badge always resolves to the latest one. No per-release step
  is needed (requires the repository to be enabled once at
  <https://zenodo.org/account/settings/github/>).
