#include <LatAnalyze/Core/Math.hpp>
#include <LatAnalyze/Core/OptParser.hpp>
#include <LatAnalyze/Io/Io.hpp>
#include <LatAnalyze/Statistics/Dataset.hpp>

using namespace std;
using namespace Latan;

#ifndef DEF_NSAMPLE
#define DEF_NSAMPLE 2000
#endif


int main(int argc, char *argv[])
{
  // parse command line/options ////////////////////////////////////////////////
  OptParser              opt;
  bool                   parsed;
  string                 g, L, m2;
  double                 N;
  
  opt.addOption("" , "help"      , OptParser::OptType::trigger, true,
                "show this help message and exit");
  parsed = opt.parse(argc, argv);
    if (!parsed or (opt.getArgs().size() != 3) or opt.gotOption("help"))
    {
        cerr << "Usage: " << argv[0] << " <g> <L> <m^2> " << endl;

        cerr << endl
             << "Possible options:" << endl
             << opt << endl;

        return EXIT_FAILURE;
    }
  
    g  = opt.getArgs()[0];
    L  = opt.getArgs()[1];
    m2 = opt.getArgs()[2];
    N  = 2.0;

    // Data prefix ///////////////////////////////////////////////////////
    string            corrDir, infile_prefix, outfile_prefix, outfile_dir;
    vector<string>    flowtime_list;
    vector<Index>     q_list;
    
    corrDir           = "data/processed_data/fft_correlator/";
    infile_prefix     = corrDir+"fft_g"+g+"_L"+L+"_m2"+m2+"_wilson_twopt_";
    outfile_prefix    = "c3_g"+g+"_L"+L+"_m2"+m2;
    outfile_dir       = "data/processed_data/flowtime_fit_data";
    
    // specify range of flowtime and momentum for ensembles////////////////

    if (g == "0.1")
    {
        flowtime_list = {"1.0", "2.0", "3.0", "4.0", "5.0", "6.0", "7.0", "8.0", "9.0"};
    } else if (g == "0.2")
    {
        flowtime_list = {"1.0", "1.4", "1.8", "2.2", "2.6", "3.0", "3.4", "3.8", "4.2"};
    } else if (g == "0.3")
    {
        flowtime_list = {"1.0", "1.3", "1.6", "1.9", "2.2", "2.5", "2.8", "3.1", "3.4"};
    }

    if (L == "64")
    {
        q_list = {1};
    } else if (L == "128")
    {
        q_list = {2};
    } else if (L == "256")
    {
        q_list = {4};
    }
    
    // read and re-arrange data//////////////////////////////////////////
    Index                   nSample = 2000;
    DMatSample              out(nSample, flowtime_list.size(), 1), emtc, trphi;
    vector<DMatSample>      out_dir_avg(q_list.size(), out);

    // direction 0
    for (Index q=0; q < q_list.size(); ++q){
        for (Index t=0; t < flowtime_list.size(); ++t)
        {
            emtc  = Io::load<DMatSample>(infile_prefix+"emtc_0_0_wilson_"+flowtime_list[t]+"_trphi_2_"+to_string(q_list[q])+"_0.h5");
            trphi = Io::load<DMatSample>(infile_prefix+"trphi_2_wilson_"+flowtime_list[t]+"_trphi_2_"+to_string(q_list[q])+"_0.h5");
            for (Index s = central; s < nSample; ++s)
            {
                out[s](t) = emtc[s](0)/(N*trphi[s](0));
                out_dir_avg[q][s](t) = out[s](t)/3;
            }
        }
        ;
    }

    // direction 1
    for (Index q=0; q < q_list.size(); ++q){
        for (Index t=0; t < flowtime_list.size(); ++t)
        {
            emtc  = Io::load<DMatSample>(infile_prefix+"emtc_1_1_wilson_"+flowtime_list[t]+"_trphi_2_0_"+to_string(q_list[q])+".h5");
            trphi = Io::load<DMatSample>(infile_prefix+"trphi_2_wilson_"+flowtime_list[t]+"_trphi_2_0_"+to_string(q_list[q])+".h5");
            for (Index s = central; s < nSample; ++s)
            {
                out[s](t) = emtc[s](0)/(N*trphi[s](0));
                out_dir_avg[q][s](t) += out[s](t)/3;
            }
        }
    }

    // direction 2
    for (Index q=0; q < q_list.size(); ++q){
        for (Index t=0; t < flowtime_list.size(); ++t)
        {
            emtc  = Io::load<DMatSample>(infile_prefix+"emtc_2_2_wilson_"+flowtime_list[t]+"_trphi_2_0_0.h5");
            trphi = Io::load<DMatSample>(infile_prefix+"trphi_2_wilson_"+flowtime_list[t]+"_trphi_2_0_0.h5");
            for (Index s = central; s < nSample; ++s)
            {
                out[s](t) = emtc[s](q_list[q])/(N*trphi[s](q_list[q]));
                out_dir_avg[q][s](t) += out[s](t)/3;
            }
        }

    }

    // direction average
    for (Index q=0; q < q_list.size(); ++q){
        Io::save<DMatSample>(out_dir_avg[q], outfile_dir+"/"+outfile_prefix+"_flowtime.h5");
    }


  return EXIT_SUCCESS;
}
