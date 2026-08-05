#!/bin/bash
# Push and VERIFY by comparing the remote ref against local HEAD, rather than
# trusting the push command's output: `git push ... 2>&1 | tail -1` both discards
# the exit status and surfaces a macOS keychain warning ("failed to store: -50")
# that reads like a git failure but is not.
set -uo pipefail
REPO="${REPO:-PaulSchnacke/ubyw_charged_e2_af3}"
BRANCH="${BRANCH:-main}"
URL="https://x-access-token:${GITHUB_TOKEN}@github.com/${REPO}.git"
git push "$URL" "$BRANCH"; push_status=$?
local_head=$(git rev-parse HEAD)
remote_head=$(git ls-remote "$URL" "refs/heads/$BRANCH" | awk '{print $1}')
if [ "$local_head" = "$remote_head" ]; then
  echo "VERIFIED: remote $BRANCH == local HEAD ${local_head:0:7}"; exit 0
fi
echo "PUSH NOT CONFIRMED (git push exit=$push_status)" >&2
echo "  local  ${local_head:0:7}" >&2
echo "  remote ${remote_head:0:7}" >&2
exit 1
