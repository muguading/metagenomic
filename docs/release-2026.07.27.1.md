# Release 2026.07.27.1

This release is identified by Git tag `release-2026.07.27.1`. Deploy the tag,
not a moving branch.

## Frozen inputs

- Application source: the annotated Git tag.
- Python packages: `requirements-release.lock`.
- Portal runtime configuration: a server-local `deployment/portal.env` copied
  from `deployment/portal.env.example`; do not commit it.
- Conda environments and reference databases: record their exact asset version
  and checksum in the delivery manifest before installation.

## Verification

```bash
git checkout release-2026.07.27.1
python -m venv .venv_web
./.venv_web/bin/python -m pip install -r requirements-release.lock
./.venv_web/bin/python -m pytest tests
./.venv_web/bin/python scripts/check_deployment.py \
  --env-file deployment/portal.env --profile full-analysis --strict-assets
```

The last command requires a real server configuration and installed data
assets; the template file is intentionally not a passing deployment input.
