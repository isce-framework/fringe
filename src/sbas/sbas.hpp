#ifndef FRINGE_SBAS_H
#define FRINGE_SBAS_H

#include <iostream>
#include <vector>
#include <armadillo>


struct Scene
{
    std::string date;
    std::string bperpName;

    Scene(){};
    ~Scene(){};

    void print()
    {
        std::cout<< "Date: " << date << "\n";
        std::cout<< "Bperp: " << bperpName << "\n";
    }
};


struct Pair
{
    //Master date
    std::string masterDate;

    //Slave date
    std::string slaveDate;

    //Filename for interferogram / offsets
    std::string ifgName;

    //Filename for coherence / quality
    std::string cohName;

    //Threshold to apply to coherence / quality
    float threshold;

    //Scale factor to multiply input data with
    float scale;

    //Values that will be determined in C++ code
    //Reference offset to take out of observations before scaling
    float referenceOffset;

    //Time difference in years
    float deltaT;

    //Constructors / destructors
    Pair():threshold(0.0f), scale(1.0f){};
    ~Pair(){};

    void print()
    {
        std::cout << "Pair: " << masterDate << " - " << slaveDate << std::endl;
        std::cout << "Observation: " << ifgName << std::endl;
        
        if (!cohName.empty()) 
        {
            std::cout << "Quality:  " << cohName << std::endl;
            std::cout << "Threshold: " << threshold << std::endl;
        }

    };

};


struct sbasOptions
{
    //Number of pairs
    int nPairs;

    //All the pairs
    std::vector<Pair> pairs;

    //All the scenes
    std::vector<Scene> scenes;

    //Reference box
    std::vector<int> refbox;

    //Valid region
    std::vector<int> bbox;

    //Output directory
    std::string outDir;

    //Reference date
    std::string refDate;

    //Values that will be determined in C++ code
    //Number of unique images
    int nSAR;

    //Dates
    std::vector<std::string> dates;

    //DEM error related stuff
    bool estimateDEMError;
    double startingRange;
    double rangeSpacing;
    double wavelength;
    std::string incAngleFile; 


    //Constructors/destructors
    sbasOptions(): refbox(4), bbox(4), estimateDEMError(false) {};
    ~sbasOptions(){};

    void print()
    {
        std::cout << "Number of pairs : " << nPairs << std::endl;
        std::cout << "Reference box: " << refbox[0] << " " << refbox[1] 
                  << " " << refbox[2] << " " << refbox[3] << std::endl;

        std::cout << "Valid region: " << bbox[0] << " " << bbox[1] << " "
                  << bbox[2] << " " << bbox[3] << std::endl;

        std::cout << "Reference date: " << refDate << std::endl;

        std::cout << "Output directory: " << outDir << std::endl;
    }
};


struct Solution
{
    arma::fvec tsest;
    double temporalCorrelation;
    double zerr;
    double velocity;
    bool rewrap;
};

struct Inverter
{
    int nPairs;
    int nSAR;
    int refIndex;
    bool fullRank;

    //Connecitivity matrix related
    arma::fmat Jmat;
    arma::fmat Jmatinv;
    arma::fvec epochs;

    //DEM error / model fit related
    arma::fmat Gmat;
    arma::fmat Gmatinv;
    arma::fmat GmatinvTrans;
    arma::fmat Pmat;

    void setup(int nPair, int nsar, int order)
    {
        Jmat.resize(nPair, nsar);
        Jmat.zeros();
        
        epochs.resize(nsar);
        Gmat.resize(nsar, order+1);

        nPairs = nPair;
        nSAR = nsar;
    }

    void prepare()
    {
        //This should only be called once in the code

        //Setting up the SBAS matrix
        Jmat.shed_col(refIndex);
        fullRank = (arma::rank(Jmat) == (nSAR-1));
        
        if (!fullRank)
        {
            std::cout << "Warning :: Network graph is not connected \n";
            throw std::invalid_argument("Network graph is not connected");
        }
        else
        {
            std::cout << "Network graph is connected. \n";
            arma::mat Jmatdbl  = arma::conv_to<arma::mat>::from(Jmat);
            arma::mat Jmatdblinv = arma::pinv(Jmatdbl);

            Jmatinv = arma::conv_to<arma::fmat>::from(Jmatdblinv);
        }


        //Setting up Bperp matrix for DEM error
        Gmat.col(0).fill(1);
        for (int ii=1; ii<Gmat.n_cols; ii++)
        {
            Gmat.col(ii) = Gmat.col(ii-1) % epochs;
        }

        Gmat.shed_row(refIndex);

        if (arma::rank(Gmat) != Gmat.n_cols)
        {
            std::cout << "Warning :: Design matrix is not full rank for dem error inversion. \n";
            throw std::invalid_argument("Poor design matrix for DEM error estimation");
        }
        {
            arma::mat Dmat = arma::conv_to<arma::mat>::from(Gmat);
            arma::mat Dmatinv = arma::pinv(Dmat);
            Gmatinv = arma::conv_to<arma::fmat>::from(Dmatinv);
            GmatinvTrans = arma::trans(Gmatinv);
            Pmat = arma::conv_to<arma::fmat>::from(Dmat * Dmatinv);
        }

    };


    void solve(const arma::fvec &inphs, Solution* worker) const
    {
        worker->tsest = Jmatinv * inphs;
        
        std::complex<float> unit(0., 1.);
        std::complex<float> tcorr = arma::sum(arma::exp(unit*(inphs - Jmat * worker->tsest)));

        worker->temporalCorrelation = std::abs(tcorr)/(1.0*nPairs);
        worker->rewrap=false;

    };


    void estimateDEMError(const arma::fvec &bperpfact, Solution* worker) const
    {
        //This is the pseudo-inverse update part
        //Computation of pseudo-inverse without repeated matrix operations.
        arma::fvec Ax = Gmatinv * bperpfact;
        arma::fvec Px = Pmat * bperpfact;
        arma::fvec b(Gmatinv.n_cols);

        float alpha = arma::dot(bperpfact, Px);
        if (alpha < 1.0e-9)
        {
            float eta = arma::dot(Ax, Ax);
            b = GmatinvTrans * (Ax / (1.0 + eta));
        }
        else
        {
            b = Px/alpha;
        }

        worker->zerr = arma::dot(b, worker->tsest);
        arma::fvec model = Gmatinv * worker->tsest -  Ax * worker->zerr;
        worker->velocity = model(1);

        //Return cleaned up time-series
        worker->tsest = worker->tsest - bperpfact * worker->zerr;
    };


    void estimateVelocity(Solution *worker) const
    {
        arma::fvec model =Gmatinv * worker->tsest;
        worker->velocity = model(1);
    };

    Inverter(){};
    ~Inverter(){};
};




#endif
