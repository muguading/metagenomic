<!-- headroom:rtk-instructions -->
# RTK (Rust Token Killer) - Token-Optimized Commands

When running shell commands, **always prefix with `rtk`**. This reduces context
usage by 60-90% with zero behavior change. If rtk has no filter for a command,
it passes through unchanged — so it is always safe to use.

## Key Commands
```bash
# Git (59-80% savings)
rtk git status          rtk git diff            rtk git log

# Files & Search (60-75% savings)
rtk ls <path>           rtk read <file>         rtk grep <pattern>
rtk find <pattern>      rtk diff <file>

# Test (90-99% savings) — shows failures only
rtk pytest tests/       rtk cargo test          rtk test <cmd>

# Build & Lint (80-90% savings) — shows errors only
rtk tsc                 rtk lint                rtk cargo build
rtk prettier --check    rtk mypy                rtk ruff check

# Analysis (70-90% savings)
rtk err <cmd>           rtk log <file>          rtk json <file>
rtk summary <cmd>       rtk deps                rtk env

# GitHub (26-87% savings)
rtk gh pr view <n>      rtk gh run list         rtk gh issue list

# Infrastructure (85% savings)
rtk docker ps           rtk kubectl get         rtk docker logs <c>

# Package managers (70-90% savings)
rtk pip list            rtk pnpm install        rtk npm run <script>
```

## Rules
- In command chains, prefix each segment: `rtk git add . && rtk git commit -m "msg"`
- For debugging, use raw command without rtk prefix
- `rtk proxy <cmd>` runs command without filtering but tracks usage
<!-- /headroom:rtk-instructions -->

## Test server for build and packaging

Use this server only for project-scoped environment setup, builds, packaging, and
non-production verification:

- Host: `111.231.101.166`
- User: `ubuntu`
- Platform: Ubuntu (confirm the exact version and CPU architecture after login)

### Credential handling

- Do **not** store passwords, private keys, access tokens, or database secrets in
  this repository, commits, generated reports, or command output.
- Before any password-based SSH operation, obtain the password from the current
  secure session and pass it through an ephemeral secret mechanism (for example,
  a terminal prompt or `SSH_PASSWORD` supplied only to the current process).
- Prefer a dedicated SSH key or a short-lived credential for repeated work.

### Operating rules

- Treat this host as a test server, never as production.
- Inspect disk space, OS release, architecture, existing Conda/Docker tooling,
  and active services before installing or changing anything.
- Keep project assets under a dedicated project directory; do not modify unrelated
  users, services, Conda environments, Docker images, or databases.
- Ask for confirmation before destructive actions, public network exposure,
  system-wide package changes, service restarts, or deleting remote assets.
- Record installed versions, environment names, database locations, and validation
  results in a project-local deployment report that contains no secrets.
