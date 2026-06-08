#!/bin/zsh
set -euo pipefail

ROOT="/Users/wuhhh/Desktop/徐老师/代码/metagenomic"
PYTHON="$ROOT/.venv_web/bin/python"
SCRIPT="$ROOT/scripts/run_demo_portal.py"
URL="http://127.0.0.1:5088/login"

cd "$ROOT"
echo "Starting demo portal..."
echo "URL: $URL"
echo ""
echo "Keep this window open while demonstrating."
echo ""

"$PYTHON" "$SCRIPT" &
SERVER_PID=$!

cleanup() {
  if kill -0 "$SERVER_PID" >/dev/null 2>&1; then
    kill "$SERVER_PID" >/dev/null 2>&1 || true
  fi
}

trap cleanup EXIT INT TERM

for _ in {1..30}; do
  if curl -sf "$URL" >/dev/null 2>&1; then
    open "$URL"
    wait "$SERVER_PID"
    exit 0
  fi
  sleep 1
done

echo "Demo portal did not become ready in time."
exit 1
