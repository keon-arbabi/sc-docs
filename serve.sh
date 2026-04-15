#!/usr/bin/env bash
set -euo pipefail

DOCS_DIR="$(cd "$(dirname "$0")/build/html" && pwd)"
PORT=8000
LOG=$(mktemp)

cleanup() {
    printf "\nShutting down...\n"
    kill "$HTTP_PID" "$TUNNEL_PID" 2>/dev/null
    wait "$HTTP_PID" "$TUNNEL_PID" 2>/dev/null
    rm -f "$LOG"
    exit 0
}
trap cleanup INT TERM

# Download cloudflared if missing
if [[ ! -x ~/bin/cloudflared ]]; then
    echo "Installing cloudflared..."
    mkdir -p ~/bin
    wget -q https://github.com/cloudflare/cloudflared/releases/latest/download/cloudflared-linux-amd64 -O ~/bin/cloudflared
    chmod +x ~/bin/cloudflared
fi

python -m http.server "$PORT" -d "$DOCS_DIR" &>/dev/null &
HTTP_PID=$!

~/bin/cloudflared tunnel --url "http://localhost:$PORT" >"$LOG" 2>&1 &
TUNNEL_PID=$!

# Wait for the URL to appear (timeout after 15s)
for i in $(seq 1 15); do
    if URL=$(grep -o 'https://[^|[:space:]]*trycloudflare.com' "$LOG" 2>/dev/null | head -1) && [[ -n "$URL" ]]; then
        break
    fi
    sleep 1
done

if [[ -n "${URL:-}" ]]; then
    printf "\n  %s\n\n  Ctrl+C to stop\n\n" "$URL"
else
    echo "Failed to get tunnel URL. Log:"
    cat "$LOG"
    cleanup
fi

# Keep running until Ctrl+C
while kill -0 "$TUNNEL_PID" 2>/dev/null; do
    wait "$TUNNEL_PID" 2>/dev/null || true
done
