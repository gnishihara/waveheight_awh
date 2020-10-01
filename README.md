# Estimate wave heights using a pressure data logger.

Calculate wave height using the zero-crossing method for the AWH data loggers.
The code to calculate wave heights with the zero-crossing method is heavily
borrowed from the oceanwaves R package.
I just modified the wave calculating function so that it does not calculate the
statistics internally. 
This is because of the way the AWH logger was setup to record pressure.

Note: I think the AWH logger pressure recordings need to be changed.
For example the data analyzed was from a logger that was setup with the following settings:

* Sampling interval: 0.1 seconds
* Sampling period: 10 seconds.
* Burst interval: 60 seconds

So every minute a 10 seconds of data will be taken at 10 Hz. 
For a total of 100 data points.
This setup will not measure waves with periods longer than 10 seconds.
Perhaps, the sampling period should be increased to 60 seconds.
Then the burst interval should be increased to 300 seconds (5 minutes) or more to save power.

The wave periods in coastal areas are generally less than 20 seconds.

References
1. <http://www.coastalwiki.org/wiki/Waves>
2. <https://pubs.usgs.gov/of/2002/of02-206/phy-environment/n-wave-climate.html>
