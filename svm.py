from sklearn import svm
from machine_learning_data import Data, DataSet

data = Data(ml_data)

# This is a regression problem
# We want to predict MRI pixel intensity from the transcript profiles at that coordinate.

training_X = data.training.input_vector
training_Y = data.training.label

validation_X = data.validation.input_vector
validation_Y = data.validation.label

clf = svm.SVC()
clf.fit(X,Y)

test_X = data.test.input_vector
test_Y = data.test.label