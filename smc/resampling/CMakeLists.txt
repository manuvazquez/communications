include_directories(.)
ADD_LIBRARY (resampling STATIC
		BestParticlesResamplingAlgorithm
		MultinomialResamplingAlgorithm
		ResamplingAlgorithm
		ResidualResamplingAlgorithm
		WithoutReplacementResamplingAlgorithm
		WithThresholdResamplingAlgorithmWrapper)