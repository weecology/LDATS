## Resubmission
This is a resubmission. In this version I have:

* Substantially reduced runtime of the full check by reducing iterations in tests
* Added Weecology (the copyright holder) to the DESCRIPTION file 
* Capitalized only the sentence beginnings and names in the description text 
* Elaborated the Description field in the DESCRIPTION file including references as requested 
* Added a \value entry with explanation to all functions' .Rd files with 
* Added small executable examples to all functions' .Rd files 
* Longer running examples have been enclosed in \donttest{}, and all \dontrun{} wrappings have been replaced with \donttest{} 
* Ensured re-setting of user options (set via par()) using on.exit() 
* All messages to the console are now written using message()/warning() rather than print()/cat() 

## Test environments
* local Windows 10 home install, R 3.6.1 64-bit and 32-bit 
* ubuntu 14.04.5 LTS (on travis-ci), R 3.6.1 and R-devel (2019-07-24 r76881) 
* win builder, R 3.6.1 and R-devel (2019-07-05 r76784) 
* R-hub builder, Ubuntu Linux 16.04, R-release 
* R-hub builder, Fedora Linux, R-devel 

## R CMD check results:
There were no ERRORs or WARNINGs 
There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Juniper L. Simonis <juniper.simonis@weecology.org>'

New submission

Possibly mis-spelled words in DESCRIPTION:
  Blei (26:68)
  Kleykamp (27:65)
  Venables (28:36)
  al (27:8, 29:55)
  changepoints (26:5)
  et (27:5, 29:52)


  Those words are spelled correctly.

## Downstream dependencies
There are currently no downstream dependencies for this package.
