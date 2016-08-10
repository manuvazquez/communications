include_directories(. 
		    ${PROJECT_SOURCE_DIR}/util
		    ${PROJECT_SOURCE_DIR}/filters
		    ${PROJECT_SOURCE_DIR}/statistics
		    ${PROJECT_SOURCE_DIR}/random
		    ${PROJECT_SOURCE_DIR}/transmission
		    ${PROJECT_SOURCE_DIR}/estimation)
ADD_LIBRARY (detection STATIC
		      DecorrelatorDetector
		      LinearDetector
		      MMSEDetector
		      SOSMMSEDetector
		      EmbeddedICSOSMMSEDetector
		      RMMSEDetector)
add_subdirectory(trellis)