# MDA public website

This repository contains the static public website for **MDA (Medical Data Analysis)**.

## Public boundary

- The site explains the project, its analysis workspaces, validation evidence, and release status.
- It does not contain the private analysis engine, statistical adapters, internal validation fixtures, or research data.
- It does not upload files or run statistical analyses in the browser.

## Local preview

From this directory:

```bash
python3 -m http.server 4173
```

Then open `http://127.0.0.1:4173`.

## Publishing

GitHub Pages publishes the `main` branch from the repository root. Public installers, once ready, belong in GitHub Releases rather than the website source tree.

Before an installer is made public, it should pass the release build, signing, notarization, clean-machine installation, and checksum checks described on the site.
