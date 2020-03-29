#include "Fit.h"


int main(int argc, char* argv[])
{

   std::cout << "fitting a time series" << std::endl;

   str inDs = "/u/k-data/fattahi/Kurdistan/TimeSeries/72/Sequential/TimeSeries_Post/timeSeries.vrt"; 

   //str tempCohDs = "/u/k-data/fattahi/Kurdistan/TimeSeries/72/Sequential/Datum_connection/EVD/tcorr.bin";
   str tempCohDs = "/u/k-data/fattahi/Kurdistan/TimeSeries/72/Sequential/TimeSeries_Post/fit/tcorr.bin";
   //str inDs = "/u/k-data/fattahi/Kurdistan/TimeSeries/72/Sequential/TimeSeries_Post/fit/54_18_looks/timeSeries_54r_18az.vrt";

   str outDs = "/u/k-data/fattahi/Kurdistan/TimeSeries/72/Sequential/TimeSeries_Post/fit/fitParameters.envi";

   str tauDs = "/u/k-data/fattahi/Kurdistan/TimeSeries/72/Sequential/TimeSeries_Post/fit/EstimatedTau.bil";   

   double tEq = 2017.8624229979466;
   str modelType = "linear_coseismic_postseismic";
   int lines_per_block=256;
   int Npar = 4 ;
   float tcorThresh = 0.5;

   Fit fitObj;
   fitObj.OpenDataset(inDs, tempCohDs);
   fitObj.get_time();
   fitObj.set_Npar(Npar);

   fitObj.set_blocksize(lines_per_block);   
   fitObj.set_outputDataset(outDs);
   fitObj.set_earthquake_time(tEq);
   fitObj.set_temporal_coherence_threshold(tcorThresh);
   fitObj.set_tauDataset(tauDs);

   // estimate Tau from a range of Taus
   fitObj.estimate_Tau(modelType);

   /*
   // create an ouput dataset
   fitObj.CreateDataset();

   // estimate other parameters with Tau fixed from previous step
   fitObj.fit_timeseries(modelType);
   */
   fitObj.CloseDataset();

}

