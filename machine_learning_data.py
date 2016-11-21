class Data(object):
	"""
	This data class stores the input features and labels for the machine learning algorithms

	It consists of three data sets (training, validation, test)
	"""

	def __init__(self, data):
		super(ClassName, self).__init__()
		
		training,test,validation = split_data(data)

		self.training = DataSet(training)
		self.test = DataSet(test)
		self.validation = DataSet(validation)

class DataSet(object):
	"""
	Each data set consists of two numpy arrays, one for inputs and one for labels.

	Initially, the inputs will be an array of arrays of transcriptome product values
	and the labels will be an array of arrays of mri intensities.
	"""
	
	def __init__(self, data):
		super(ClassName, self).__init__()
		self.data = munge_data(data)
