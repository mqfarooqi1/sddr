## Resubmission

This is a resubmission (0.1.1). In response to the review of 0.1.0 I have:

* Removed the invalid DOI (10.1016/0014-2921(93)90031-3) from the Quah (1993)
  reference; the reference is kept as plain text.
* Note on the flagged "possibly misspelled words": 'Quah' and 'Rey' are author
  surnames and 'ergodic' is a standard technical term; all are spelled
  correctly.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release, so there is a NOTE for "New submission".

## Test environments

* Local: Windows 11, R 4.5.2
* GitHub Actions: ubuntu-latest (R-devel, release, oldrel-1),
  macOS-latest (release), windows-latest (release)

## Notes

* The reference implementation (the Python library 'giddy', via 'reticulate')
  is used only in optional, skip-guarded tests for numerical validation. It is
  not required to install, load, use, or check the package, and those tests are
  skipped on CRAN.
* 'spdep' is used only when a user supplies spatial weights as a 'listw' or
  'nb' object; a plain weights matrix needs no extra packages.
