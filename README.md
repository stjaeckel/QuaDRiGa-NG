# QuaDRiGa-NG

This is a fork of the original QuaDRiGa repository (v2.6.1) hosted by the Fraunhofer HHI (https://github.com/fraunhoferhhi/QuaDRiGa). 

QuaDRiGa, short for QUAsi Deterministic RadIo channel GenerAtor, is used for generating realistic radio channel impulse responses for system-level simulations of mobile radio networks. These simulations are used to determine the performance of new digital-radio technologies in order to provide an objective indicator for the standardization process in bodies like the third generation partnership program 3GPP.

I, Stephan Jaeckel, PhD, wrote most of the ~80.000 lines of code during my time at Fraunhofer from 2011 to 2020. QuaDRiGa was the topic of my PhD thesis (https://www.db-thueringen.de/receive/dbt_mods_00032895) and numerous publications on the topic. I left Fraunhofer in 2020 to continue my career as a freelance consultant in the industry. Since then, I am no longer responsible for maintaining QuaDRiGa at Fraunhofer. I do not know what Fraunhofer's plans for QuaDRiGa are nor who will be responsible for technical support in the future. Since I offer support on QuaDRiGa as a freelance service for multiple contacts and use it for university teaching, I created this fork to maintain the code, provide bug fixes and develop new features.

For questions, bug reports, feature request, etc. please contact me at quadriga@sjc-wireless.com, through GitHub or LinkedIn (https://www.linkedin.com/in/stephan-jaeckel-phd-a66990a4). 

QuaDRiGa v2.6.1 is released under the
Software License for The QuaDRiGa Channel Model  
© Copyright 2011 - 2021 Fraunhofer-Gesellschaft zur Förderung der angewandten Forschung e.V., All rights reserved.

## Current Version

The current experimental version requires the ``arrayant-lib`` library to work. This library implements some core functionality in C++/MEX for faster execution and interfaces for future projects. The library can be found here: https://github.com/stjaeckel/arrayant-lib. After installing the library, calling ``qd_simulation_parameters`` on the MATLAB/Octave command prompt should show the QuaDRiGa version 3.0.0 and ``arrayant-lib``-version 0.1.0. If it shows an error, you might need to recompile the library.

List of changes:
* v3.0.0 - Replacing the ``qd_array.interpolate`` method with the one provided by the ``arrayant-lib`` library improves performance by up to 70% on multi-core CPUs


## Future plans

I plan to extend the existing QuaDriGa code in several ways. However, this depends on the requirements from projects, contracts and my available time. The first extensions will be:

* Array antenna processing functions written in C++ with MEX and Python wrappers
* Sum-of-sinusoids random number generator written in C++ with MEX and python wrappers

These functions contain the most computing intense parts of QuaDRiGa and may use up to 80% of the computing time. Providing native C++ code with potential (optional) NIVIDIA-CUDA and multi-core acceleration trough OpenMP can provide enormous performance gains. Also, these parts can be valuable components for other software projects unrelated to QuaDRiGa (such as RayTracing or DataScience applications). Hence, the new functionality will be provided by a set of separate libraries released under a Free-Software license. Integration with the existing QuaDRiGa code will be realized through a MEX wrapper. Native python wrappers will also made available. 

Further ideas are:

* QuaDRiGa support library (written in C++) for common functions such as parameter generation, channel generation, data import/export, channel interpolation
* Python implementation of the whole framework
* Integration of Ray-Tracing functions for site-specific simulations
* Support for latest 3GPP channel models
* Support for MATLAB, Octave and Python on Windows and Linux platforms
* WiKi containing the documentation
* Community contributions through GitHUB

I am currently working on this in my free time with limited funding and can therefore not provide a definitive timeline. Please contact me if you have questions or comments.
