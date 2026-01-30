===============
Release History
===============

Initial Release (2025-04-05)
----------------------------
Provisional Release

Version 1.0.0 (2025-06-24)
--------------------------
Official Software Release

Version 1.0.5 (2025-11-17)
--------------------------
PyPi and conda-forge official releases

Version 1.1.0 (2025-11-25)
--------------------------
* Additional options for AnalysisProject
    * Add option to allow instantiation of AnalysisProject by sending
        a dictionary of properties directly rather than requiring
        an external yml file.
    * Add option to AnalysisProject to suppress writing to disk

Version 1.1.1 (2025-11-30)
--------------------------
* Debugging options for AnalysisProject
    * allowing for `write_results_to_files` even without a yml file. Results will be written to `default` as a root rather than the YML filename in this case.


Version 1.1.2 (2025-12-06)
--------------------------
* Updates for Ward-Lough from PR #69
    * change variable name for time from `t` to `time` for consistency with other solultions
    * accomodate a single scalar time at which to calculate

Version 1.2.0 (2026-01-10)
--------------------------
* Refactoring to properly name the Dudley Ward and Lough functions
    * Recognized we were erroneously referring to the Dudley Ward and Lough solution as simply Ward and Lough. Refactored to change the spelling in comments and function names.
    * solution `ward_lough_drawdown` is now `dudley_ward_lough_drawdown`
    * solution `ward_lough_depletion` is now `dudley_ward_lough_depletion`

Version 1.2.5 (2026-01-20)
--------------------------
* Performance enhancements to Hunt (2003) solutions. 
    * precalculation of some constants and tuning of numerical integral
    * vectorization of some calculations
    * no syntax changes necessary for use

Version 1.3.0 (2026-01-30)
--------------------------
* Performance enhancements to all depletion solutions. 
    * precalculation of a unit pumping response precedes multiplication rather than full depletion calculation to fill out pumping time series
    * no syntax changes necessary for use
