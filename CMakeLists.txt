cmake_minimum_required(VERSION 3.16.3...3.19.7 FATAL_ERROR)

project(CoronaryCaScore)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://www.slicer.org/wiki/Documentation/Nightly/Extensions/CoronaryCaScore")
set(EXTENSION_CONTRIBUTORS "Vittorio Censullo (AITeRTC)")
set(EXTENSION_DESCRIPTION "Coronary Artery Calcium Score - Automated CAC scoring with per-vessel analysis and MESA percentiles")
set(EXTENSION_ICONURL "https://www.example.com/Slicer/Extensions/CoronaryCaScore.png")
set(EXTENSION_SCREENSHOTURLS "https://www.example.com/Slicer/Extensions/CoronaryCaScore/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a list or "NA" if no dependencies

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(CoronaryCaScore)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
