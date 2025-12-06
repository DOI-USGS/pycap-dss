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

