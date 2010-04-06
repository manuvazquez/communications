include_directories(.
		    ${PROJECT_SOURCE_DIR}/util
		    ${PROJECT_SOURCE_DIR}/random)
ADD_LIBRARY (statistics STATIC
			ARprocess
			StatUtil
			UsersActivityDistribution)
TARGET_LINK_LIBRARIES (statistics random)