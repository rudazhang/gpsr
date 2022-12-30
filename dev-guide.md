# Development Guide

Release:
- `git tag <tagname>`, create tag;
- update version in `DESCRIPTION`;
- `devtools::build()`, build a bundled package (.tar.gz) to parent directory.
- `devtools::build(binary = TRUE)`, build a binary package (.tar.gz) to parent directory.
- Create release on GitHub, upload bundled and binary packages.
- Finish. GitHub automatically create source code releases (.zip, .tar.gz).

Install:
- `devtools::install()`, build a bundled package (to `/tmp`) and (re-)install.
