include_directories(. 
		    ${PROJECT_SOURCE_DIR}/util)
ADD_LIBRARY (detection STATIC
		      DecorrelatorDetector
		      LinearDetector
		      MMSEDetector
		      RMMSEDetector)
add_subdirectory(trellis)