include_directories(. 
		    ${PROJECT_SOURCE_DIR}/estimation
		    ${PROJECT_SOURCE_DIR}/transmission)
ADD_LIBRARY (trellis STATIC
		    PSPPath
		    Trellis
		    ViterbiPath)
