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
  string                 manFilename, g, L, m2, corrDir;

  
  opt.addOption("" , "help"      , OptParser::OptType::trigger, true,
                "show this help message and exit");
  parsed = opt.parse(argc, argv);
    if (!parsed or (opt.getArgs().size() != 3) or opt.gotOption("help"))
    {
        cerr << "NEED TO UPDATE: usage: " << argv[0] << " g, L, m2 " << endl;
        cerr << "NO Parameter file fit symbols: "
             << "NO  x_0, m2, L, t (flowtime), mcrit " << endl;
        cerr << endl
             << "Possible options:" << endl
             << opt << endl;

        return EXIT_FAILURE;
    }
  
    g = opt.getArgs()[0];
    L = opt.getArgs()[1];
    m2 = opt.getArgs()[2];


  // load and trace data ///////////////////////////////////////////////////////
    //corrDir = "g"+g+"/L"+L+"/m2"+m2+"/";
    corrDir = "data/processed_data/fft_correlator/";
    string infile_prefix = corrDir+"fft_g"+g+"_L"+L+"_m2"+m2+"_wilson_twopt_";
    string outfile_prefix = "c3_g"+g+"_L"+L+"_m2"+m2;
    string outfile_dir = "data/processed_data/flowtime_fit_data";
    std::vector<string> flowtime_list;
    std::vector<Index>  q_list;
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
    
    cout << "g: " << g << " |    L: " << L << " |    m2: " << m2 << endl;
    Index nSample = 2000;
    DMatSample out(nSample, flowtime_list.size(), 1);
    std::vector<DMatSample> out_dir_avg(q_list.size(), out);
    DMatSample emtc;
    DMatSample trphi; 
    Index counter;

    //dir0
    for (Index q=0; q < q_list.size(); ++q){
        for (Index t=0; t < flowtime_list.size(); ++t)
        {
            emtc = Io::load<DMatSample>(infile_prefix+"emtc_0_0_wilson_"+flowtime_list[t]+"_trphi_2_"+to_string(q_list[q])+"_0.h5");
            trphi = Io::load<DMatSample>(infile_prefix+"trphi_2_wilson_"+flowtime_list[t]+"_trphi_2_"+to_string(q_list[q])+"_0.h5");
            for (Index s = central; s < nSample; ++s)
            {
                out[s](t) = emtc[s](0)/(2.0*trphi[s](0));
                out_dir_avg[q][s](t) = out[s](t)/3;
            }
        }
        ;
    }
    cout << "nsample/size = " << emtc.size() << endl;
    cout << "Complete dir 0" << endl;
    //dir1
    for (Index q=0; q < q_list.size(); ++q){
        for (Index t=0; t < flowtime_list.size(); ++t)
        {
            emtc = Io::load<DMatSample>(infile_prefix+"emtc_1_1_wilson_"+flowtime_list[t]+"_trphi_2_0_"+to_string(q_list[q])+".h5");
            trphi = Io::load<DMatSample>(infile_prefix+"trphi_2_wilson_"+flowtime_list[t]+"_trphi_2_0_"+to_string(q_list[q])+".h5");
            for (Index s = central; s < nSample; ++s)
            {
                out[s](t) = emtc[s](0)/(2.0*trphi[s](0));
                out_dir_avg[q][s](t) += out[s](t)/3;
            }
        }
    }
    cout << "Complete dir 1" << endl;

    //dir2
    for (Index q=0; q < q_list.size(); ++q){
        for (Index t=0; t < flowtime_list.size(); ++t)
        {
            emtc = Io::load<DMatSample>(infile_prefix+"emtc_2_2_wilson_"+flowtime_list[t]+"_trphi_2_0_0.h5");
            trphi = Io::load<DMatSample>(infile_prefix+"trphi_2_wilson_"+flowtime_list[t]+"_trphi_2_0_0.h5");
            for (Index s = central; s < nSample; ++s)
            {
                out[s](t) = emtc[s](q_list[q])/(2.0*trphi[s](q_list[q]));
                out_dir_avg[q][s](t) += out[s](t)/3;
            }
        }

    }
    cout << "Complete dir 2" << endl;

    for (Index q=0; q < q_list.size(); ++q){

        Io::save<DMatSample>(out_dir_avg[q], outfile_dir+"/"+outfile_prefix+"_mom_dir_avg_flowtime.h5");

    }
    cout << "Complete avg" << endl;

  return EXIT_SUCCESS;
}
