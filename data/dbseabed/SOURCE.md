# dbSEABED raw grids (tracked copy)

Four percent-composition grids for the northern Gulf of Mexico from the dbSEABED project
("Data for Modellers"), at their native 1.2 arc-minute (0.02 degree) resolution, 891 x 383
cells, NODATA = -9999, values < 0 to be treated as missing:

| Folder | File | Variable |
|---|---|---|
| `Gmf_RCK/` | `gmf_RCK_val.asc` | percent rock outcrop (an areal fraction, separate from the grain-size triangle) |
| `Gmf_GVL/` | `gmf_GVL_val.asc` | percent gravel |
| `Gmf_SND/` | `gmf_SND_val.asc` | percent sand |
| `Gmf_MUD/` | `gmf_MUD_val.asc` | percent mud |

Source: https://csdms.colorado.edu/wiki/DBSEABED#Data_for_Modellers (zips
`https://csdms.colorado.edu/csdms_wiki/images/Gmf_<CLS>.zip`). Jenkins, C. (dbSEABED,
INSTAAR, University of Colorado). Cite dbSEABED when using these layers.

Provenance of this copy: downloaded 22 Sep 2026 by `fn.pull_dbseabed()` in the
`WFS-FEM/EcospaceBasemap` repository and copied here unchanged on 2 Oct 2026 (byte-identical,
only the `.asc` grids; the zips' auxiliary `.jpg`/`.avl` files are omitted). The grids are
tracked in this repo because the CSDMS server was unreachable during the October 2026
reproducibility review, which would otherwise have blocked a fresh clone from running stage 1.
`fn.pull_dbseabed(dir.dbseabed)` in `R/_setup.R` re-downloads them if the folder is removed.
The duplication with EcospaceBasemap is tracked in issue #3.
