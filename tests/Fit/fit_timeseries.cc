#include "Fit.h"


int main(int argc, char* argv[])
{

   std::cout << "fitting a time series" << std::endl;

   str inDs = "/u/k-data/fattahi/Kurdistan/TimeSeries/72/Sequential/TimeSeries_Post/timeSeries.vrt"; 
   //str inDs = "/u/k-data/fattahi/Kurdistan/TimeSeries/72/Sequential/TimeSeries_Post/test_fit_2/54_18_looks/timeSeries_54r_18az.vrt";

   str outDs = "/u/k-data/fattahi/Kurdistan/TimeSeries/72/Sequential/TimeSeries_Post/test_fit/fitParameters_lks.envi";

   

   double tEq = 2017.8624229979466;
   str modelType = "linear_coseismic_postseismic";
   double Tau = 0.027; 
   int lines_per_block=25;
   int Npar = 4 ;


   Fit fitObj;
   std::cout << "Open Dataset" << std::endl;
   fitObj.OpenDataset(inDs);
   std::cout << "get time " << std::endl;
   fitObj.get_time();
   std::cout << "Set Npar " << std::endl;
   fitObj.set_Npar(Npar);

   fitObj.set_blocksize(lines_per_block);   
   fitObj.set_outputDataset(outDs);
   fitObj.set_earthquake_time(tEq);
   fitObj.set_Tau(Tau);

   std::cout << "Create Dataset " << std::endl;
   fitObj.CreateDataset();
   std::cout << "fit timeseries " << std::endl;
   fitObj.fit_timeseries(modelType);
   std::cout << " Closing dataset" << std::endl;
   fitObj.CloseDataset();
   std::cout << "That was it!" << std::endl;
}

