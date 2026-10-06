# dbSEABED raw grids (tracked copy)

Four percent-composition grids for the northern Gulf of Mexico from the dbSEABED project
("Data for Modellers"), at their native 1.2 arc-minute (0.02 degree) resolution, 891 x 383
cells, lower-left corner -98.19, 23.36. The header says `NODATA_value -9999` but missing
cells are coded `-99`; every reader in this repo treats values below 0 as missing. Read by
stage 1 (`fn.dbseabed_layers`) and stage 4c (`fn.load_dbseabed_raw`).

| Folder | File | Variable | MD5 (LF form) |
|---|---|---|---|
| `Gmf_RCK/` | `gmf_RCK_val.asc` | percent rock outcrop (an areal fraction, separate from the grain-size triangle) | `e1332f6d7983bcae79a5e2b134021d21` |
| `Gmf_GVL/` | `gmf_GVL_val.asc` | percent gravel | `1242a3351fc554895c5a5c413d8cf5d8` |
| `Gmf_SND/` | `gmf_SND_val.asc` | percent sand | `8147b22640263554185940f712f74cf7` |
| `Gmf_MUD/` | `gmf_MUD_val.asc` | percent mud | `16e40884cbdc192cbdb03c3012d50adb` |

Source: https://csdms.colorado.edu/wiki/DBSEABED#Data_for_Modellers (zips
`https://csdms.colorado.edu/csdms_wiki/images/Gmf_<CLS>.zip`). Jenkins, C. (dbSEABED,
INSTAAR, University of Colorado). Cite dbSEABED when using these layers.

Provenance of this copy: downloaded 22 Sep 2026 by `fn.pull_dbseabed()` in the
`WFS-FEM/EcospaceBasemap` repository (`R/dbSEABED_functions.R`) and copied here unchanged on
2 Oct 2026; only the `.asc` grids, not the zips' `.jpg`/`.avl` files. EcospaceBasemap tracks
the same four files under `data/dbseabed/` since its commit `4236720` (pull request #3), with
the same MD5s in its `SOURCE.md`. Both repositories ship them because the CSDMS server was
unreachable during the October 2026 reproducibility review (issue #2), which would otherwise
block a fresh clone from running stage 1 or stage 4c.

If the folder is removed, `git checkout -- data/dbseabed` restores it. To use another copy,
set `dir.dbseabed` in `config.local.R`, or set `dir.ecospace.basemap` to an EcospaceBasemap
clone and the driver reads its `data/dbseabed/`. The download code lives only in
EcospaceBasemap (`fn.pull_dbseabed()`); it was removed from this repository under issue #3.

MD5s are of the LF form of each file, which is how git stores them and how a checkout
delivers them (`.gitattributes` pins `*.asc` to LF). The zips from CSDMS hold CRLF files,
which hash differently (RCK `46c0bd22...`, GVL `5ce2c0b6...`, SND `bf01968f...`, MUD
`a173130e...`); a freshly downloaded copy therefore matches only after `tr -d ''` or a
commit-and-checkout. Check with
`tools::md5sum(list.files('data/dbseabed', '\.asc$', recursive=TRUE, full.names=TRUE))`
after a fresh checkout.
