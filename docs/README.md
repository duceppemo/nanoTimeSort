# Documentation sources

`wiki/` holds the source pages for the GitHub wiki
(<https://github.com/duceppemo/nanoTimeSort/wiki>). GitHub wikis live in a separate repository,
so after editing these pages, publish them with:

```bash
# One-time: create the wiki by visiting the Wiki tab on GitHub and saving the initial page,
# then:
git clone https://github.com/duceppemo/nanoTimeSort.wiki.git /tmp/nanoTimeSort.wiki
cp wiki/*.md /tmp/nanoTimeSort.wiki/
cd /tmp/nanoTimeSort.wiki
git add -A && git commit -m "Update wiki" && git push
```

`Home.md` is the wiki landing page; `[[Page]]` links map to the other file names.
