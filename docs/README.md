# Documentation sources

`wiki/` holds the source pages for the GitHub wiki
(<https://github.com/duceppemo/nanoTimeSort/wiki>) and is the **single source of truth**.

Every push to `master` that touches `docs/wiki/**` triggers the
[Sync wiki](../.github/workflows/sync-wiki.yml) workflow, which mirrors the folder to the wiki
(deleting wiki pages that no longer exist here). **Do not edit the wiki directly** — changes
made there will be overwritten by the next sync. A daily scheduled run re-drives the mirror in
case a sync ever fails.

One-time prerequisite: the wiki repository must exist, which GitHub only creates after the
first page is saved via the repository's Wiki tab (already done for this repository).

`Home.md` is the wiki landing page, `_Sidebar.md` the navigation, and `[[Page]]` links map to
the other file names.

`RELEASING.md` documents how to publish a new version to PyPI.
