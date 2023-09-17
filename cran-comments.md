This resubmission addresses problems associated with vignette building and remote resources.
We now pre-compile the paper-comparison vignette prior to the formal build process following https://ropensci.org/blog/2019/12/08/precompute-vignettes/.

## Test environments

* Local: Windows 10 home install (build 19045), 64-bit R 4.3.0 (2023-04-21 ucrt)

* GitHub Actions: MacOS 12.6.8 21G725, R 4.3.1 (2023-06-16)

* GitHub Actions: Microsoft Windows Microsoft Windows Server 2022, 10.0.20348, Datacenter, x86_64-w64-mingw32 (64-bit), R 4.3.1 (2023-06-16 ucrt)
* GitHub Actions: Microsoft Windows Microsoft Windows Server 2022, 10.0.20348, Datacenter, x86_64-w64-mingw32 (64-bit), R 3.6.3 (2020-02-29) 
* GitHub Actions: Microsoft Windows Microsoft Windows Server 2022, 10.0.20348, Datacenter, x86_64-w64-mingw32 (64-bit), R 4.1.3 (2022-03-10)

* GitHub Actions: Ubuntu 22.04.3 LTS, x86_64-pc-linux-gnu (64-bit), R 4.3.1 (2023-06-16)  
* GitHub Actions: Ubuntu 22.04.3 LTS, x86_64-pc-linux-gnu (64-bit), R 4.2.3 (2023-03-15)   
* GitHub Actions: Ubuntu 22.04.3 LTS, x86_64-pc-linux-gnu (64-bit), R 4.1.3 (2022-03-10)  
* GitHub Actions: Ubuntu 22.04.3 LTS, x86_64-pc-linux-gnu (64-bit), R 4.0.5 (2021-03-31)  
* GitHub Actions: Ubuntu 22.04.3 LTS, x86_64-pc-linux-gnu (64-bit), R 3.6.3 (2020-02-29)

* win-builder: Windows Server 2022 x64 (build 20348), x86_64-w64-mingw32, R Under development (unstable) (2023-09-16 r85157 ucrt)
* win-builder: Windows Server 2022 x64 (build 20348), x86_64-w64-mingw32, R 4.3.1 (2023-06-16 ucrt)
* win-builder: Windows Server 2022 x64 (build 20348), x86_64-w64-mingw32, R 4.2.3 (2023-03-15 ucrt)

* mac-builder: macosx, macOS 13.3.1 (22E261), r-release-macosx-arm64, R 4.3.0

* R-hub builder: Windows Server 2022, R-devel, 64 bit
* R-hub builder: 
* R-hub builder:


## R CMD check results:
There were no ERRORs, WARNINGs, or substantive NOTEs

### There are spurious NOTES associated with URLs on the win-builder system:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Juniper L. Simonis <juniper.simonis@weecology.org>'

Found the following (possibly) invalid URLs:
  URL: https://www.nsf.gov/awardsearch/showAward?AWD_ID=1315138
    From: README.md
    Status: Error
    Message: Operation timed out after 60004 milliseconds with 0 bytes received

Found the following (possibly) invalid DOIs:
  DOI: 10.1002/ecy.2373
    From: DESCRIPTION
    Status: Forbidden
    Message: 403

### There are spurious NOTES associated with URLs and the manual build on the Windows Server of R hub:

* checking CRAN incoming feasibility ... [224s] NOTE
Maintainer: 'Juniper L. Simonis <juniper.simonis@weecology.org>'
Found the following (possibly) invalid URLs:
  URL: https://www.nsf.gov/awardsearch/showAward?AWD_ID=1622425
    From: README.md
    Status: Error
    Message: libcurl error code 28:
      	Operation timed out after 60001 milliseconds with 0 bytes received

* checking HTML version of manual ... [14s] NOTE
Skipping checking math rendering: package 'V8' unavailable

* checking for non-standard things in the check directory ... NOTE
Found the following files/directories:
  ''NULL''

* checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'

## Downstream dependencies
There are currently no downstream dependencies for this package.
