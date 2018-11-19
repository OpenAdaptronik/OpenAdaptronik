# Random Decrement Analysis with Python 
## Abstract
These files enable a spectral analysis of time domain data with Python using the Random Decrement method. Using a time series of 
(usually: vibration) data as input, the RD method calculates an estimation of the autocorrelation function. 
## Environment
* Python 3.7
* Numpy
* Scipy
* PyQT5
* PyQTGraph
* Matplotlib
## Further needed
Time data! <br /> 
Either acquired with [Phyphox](https://phyphox.org/)  (great package for DIY signal acquisition with smartphones) or our 
OpenAdaptronik DIY vibration data logger (sse under "hardware")

## Files
* DIY_RDAnalysis.py: GUI 
* rd_estim.py: function to calculate the RD signature from a time series
* decimate_filter: function for a filtered downsampling of the time series
## Usage
DIY_RDAnalysis requires the time series input in a CSV file. Pre defined formats are the CSV exported by the Phyphox app or the OA datalogger

The data analysis is based on averaging triggered frames from the time series. Choose a suitable part of the time series in the upper plot. 
Then adjust the trigger level to the signal amplitude, push GO! and take a look at the resulting spectrum. 
Fiddle with the parameters (downsampling, trigger level, trigger condition) to enhance the result.

