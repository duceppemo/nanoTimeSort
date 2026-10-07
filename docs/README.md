# Documentation sources

`wiki/` holds the source pages for the GitHub wiki
(<https://github.com/duceppemo/nanoTimeSort/wiki>) and is the **single source of truth**.

Every push to `master` that touches `docs/wiki/**` triggers the
[Sync wiki](../.github/workflows/sync-wiki.yml) workflow, which mirrors the folder to the wiki
(deleting wiki pages that no longer exist here). **Do not edit the wiki directly** — changes
made there will be overwritten by the next sync.

`Home.md` is the wiki landing page, `_Sidebar.md` the navigation, and `[[Page]]` links map to
the other file names.

`RELEASING.md` documents how to publish a new version to PyPI.
