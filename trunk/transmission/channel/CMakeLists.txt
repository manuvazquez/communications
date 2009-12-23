include_directories(${PROJECT_SOURCE_DIR}/estimation)
ADD_LIBRARY (channel STATIC 
		      ARchannel
		      ARMultiuserCDMAchannel
		      BesselChannel
		      EstimatedMIMOChannel
		      MIMOChannel
		      MultiuserCDMAchannel
		      StillMemoryMIMOChannel
		      TimeInvariantChannel)
add_subdirectory(delayPowerProfile)