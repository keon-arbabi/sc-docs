#!/usr/bin/env bash
# Build the docs and deploy them to Cloudflare Pages.
#
# Stable URL (never changes across deploys):  https://brisc-docs.pages.dev
#
# One-time auth (pick ONE), then `bash deploy.sh` for every update:
#   A) interactive:  npx wrangler login
#   B) token (HPC):  export CLOUDFLARE_API_TOKEN=...   (token with "Cloudflare Pages: Edit")
#                    export CLOUDFLARE_ACCOUNT_ID=...
set -eo pipefail
cd "$(dirname "$0")"

PROJECT="brisc-docs"

eval "$(~/miniforge3/bin/conda shell.bash hook)"
conda activate brisc
make clean html

# Put the no-index files at the site root so the preview is never indexed.
cp deploy/robots.txt deploy/_headers build/html/

# Same project name -> same https://$PROJECT.pages.dev URL, updated in place.
npx --yes wrangler@latest pages deploy build/html \
  --project-name "$PROJECT" --commit-dirty=true
