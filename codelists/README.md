# `codelists/`

##Overview

This folder contains the code lists which were used to identify different
vaccines in the observation and prescription table of CPRD Aurum.
We created for categories of code lists:
 * A code list created from the medical dictionary listing neutral vaccine codes
and codes for administered vaccines (flagged)
 * A code list for all side-effects from the vaccine of interest from the 
 medical dictionary. This code list was not used in the final analysis as it 
 did not add specific information on vaccine timing.
 * A code list with terms of declined vaccination from the medical dictionary.
 * A code list for vaccine products from the product dictionary.
 
All code lists were created based on the CPRD Aurum dictionaries May 2022. 
Please note, terms need updating when using another release of CPRD Aurum data.

The antigen of interest is flagged in all the code lists (measles, pertussis,
and PCV).

##Search strategies

All code lists for the different categories of vaccination were derived from
searching the CPRD medical and product dictionaries. The search terms used can
be found in the folder 'Search_strategy'. All search results were manually 
screened and classified into the four categories.

#Use of the code lists

The use of the different code lists in combination with each other is explained
in the manuscript. 

##Caution

medcodeid and productcodeid should be treated as characters. Otherwise, R will
round the results. They also have to be imported as text or without any
prespecified class into excel, otherwise the ids will be automatically 
rounded  by excel. 