#include <LatAnalyze/Core/OptParser.hpp>
#include <LatAnalyze/Functional/CompiledModel.hpp>
#include <LatAnalyze/Io/Io.hpp>
#include <LatAnalyze/Statistics/MatSample.hpp>
#include <LatAnalyze/Core/Math.hpp>
#include <LatAnalyze/Numerical/MinuitMinimizer.hpp>
#include <LatAnalyze/Numerical/NloptMinimizer.hpp>
#include <LatAnalyze/Core/Plot.hpp>
#include <LatAnalyze/Statistics/XYSampleData.hpp>
#include <LatAnalyze/Io/XmlReader.hpp>


using namespace std;
using namespace Latan;



void setFitRange(XYSampleData &data, const Index ti, const Index tf, const Index nt)
{
    for (Index t = 0; t < nt; ++t)
    {
        data.fitPoint((t >= ti) and (t <= tf), t);
    }
}


int main(int argc, char *argv[])
{
    // parse arguments /////////////////////////////////////////////////////////
    OptParser            opt;
    bool                 parsed, doPlot, doCorr, doPaperFit, savePlot;
    string               g, L, m2, paramFileName, corrFileName, saveFileName, saveFilePrefix, savePlotPrefix, outFmt;
    Index                shift, nPar, thinning;
    double               ti, tf;
    double               svdTol;
    std::vector<double> corrFlowtimes;
    Minimizer::Verbosity verbosity;
    

    opt.addOption("d" , "default", OptParser::OptType::trigger  , true,
                  "default fit range (for paper)");    
    opt.addOption("" , "svd"      , OptParser::OptType::value  , true,
                  "singular value elimination threshold", "0.");
    opt.addOption("v", "verbosity", OptParser::OptType::value  , true,
                  "minimizer verbosity level (0|1|2)", "0");
    opt.addOption("o", "output", OptParser::OptType::trigger  , true,
                  "output file", "");
    opt.addOption("" , "uncorr"   , OptParser::OptType::trigger, true,
                  "only do the uncorrelated fit");
    opt.addOption("s", "save-plot", OptParser::OptType::trigger, true,
                    "saves the source and .pdf", "");
    opt.addOption("", "help"      , OptParser::OptType::trigger, true,
                  "show this help message and exit");
    parsed = opt.parse(argc, argv);
    if (!parsed or (opt.getArgs().size() != 1) or opt.gotOption("help"))
    {
        cerr << "usage: " << argv[0] << " <param_xml> <options> <correlator file>" << endl;
        cerr << endl << "Possible options:" << endl << opt << endl;
        
        return EXIT_FAILURE;
    }
    paramFileName = opt.getArgs().front();
    doPaperFit   = opt.gotOption("default");
    svdTol       = opt.optionValue<double>("svd");
    doCorr       = !opt.gotOption("uncorr");
    savePlot     = opt.gotOption("save-plot");
    switch (opt.optionValue<unsigned int>("v"))
    {
        case 0:
            verbosity = Minimizer::Verbosity::Silent;
            break;
        case 1:
            verbosity = Minimizer::Verbosity::Normal;
            break;
        case 2:
            verbosity = Minimizer::Verbosity::Debug;
            break;
        default:
            cerr << "error: wrong verbosity level" << endl;
            return EXIT_FAILURE;
    }
    Latan::XmlReader paramFile(paramFileName);
    corrFileName = paramFile.getFirstValue<string>("filenames", "file");
    corrFlowtimes = paramFile.getAllValues<double>("flowtimes", "t");
    if (doPaperFit)
    {
      ti = paramFile.getFirstValue<double>("defaultfitrange", "min");
      tf = paramFile.getFirstValue<double>("defaultfitrange", "max");
    }
    else
    {
      ti = paramFile.getFirstValue<double>("customfitrange", "min");
      tf = paramFile.getFirstValue<double>("customfitrange", "max");
    }
    g = paramFile.getFirstValue<string>("spacing", "ag");
    L = paramFile.getFirstValue<string>("volume", "L");
    m2 = paramFile.getFirstValue<string>("mass", "m2");
    saveFilePrefix = paramFile.getFirstValue<string>("savefiledir", "dir");
    savePlotPrefix = paramFile.getFirstValue<string>("saveplotdir", "dir");
    saveFileName = "c3-flowtime-fit-value_g"+g+"_L"+L+"_m2"+m2;
    cout << "g: " << g << " | L: " << L << " | m2 = " << m2 << " | ti: " << ti<< " | tf: " << tf  << endl;
    cout << "corrrelator file: " << corrFileName << endl;

    // load correlator /////////////////////////////////////////////////////////
    DMatSample tmp, corr;
    Index      nSample, nt;
    
    tmp     = Io::load<DMatSample>(corrFileName);
    nSample = tmp.size();
    nt      = tmp[central].rows();
    tmp     = tmp.block(0, 0, nt, 1);
    corr    = tmp;


    
    // make models /////////////////////////////////////////////////////////////
    DoubleModel mod;
    
    nPar        = 2;
    mod.setFunction([](const double *x, const double *p)
                    {
                        return p[0] + p[1]*x[0];
                    }, 1, nPar);

    
    // fit /////////////////////////////////////////////////////////////////////
    DMatSample          tvec(nSample);
    XYSampleData        data(nSample);
    DVec                init(nPar);
    NloptMinimizer      globMin(NloptMinimizer::Algorithm::GN_CRS2_LM);
    MinuitMinimizer     locMin;
    vector<Minimizer *> unCorrMin{&globMin, &locMin};
    std::vector<double> timeshift(nt);

    timeshift = corrFlowtimes;
    FOR_STAT_ARRAY(tvec, s)
    {
        //tvec[s] = DVec::LinSpaced(nt, 0, nt - 1);
        tvec[s] =  Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(timeshift.data(), timeshift.size());

    }
    data.addXDim(nt, "t/a", true);
    data.addYDim("C(t)");
    data.setUnidimData(tvec, corr);
    // set parameter name ******************************************************

    mod.parName().setName(0, "C_3");
    mod.parName().setName(1, "Omega'");   

    // set initial values ******************************************************

    init(0) = 0.05;
    init(1) = 1;
    
    // set limits for minimisers ***********************************************

    globMin.setLowLimit(0, 0.);
    globMin.setHighLimit(0, 10. * init(0));

    globMin.setLowLimit(1, 0.);
    globMin.setHighLimit(1, 10. * init(1));

    globMin.setPrecision(0.001);
    globMin.setMaxIteration(100000);
    globMin.setVerbosity(verbosity);
    locMin.setMaxIteration(1000000);
    locMin.setVerbosity(verbosity);

    // fit /////////////////////////////////////////////////////////////////////
    
      SampleFitResult fit;

      setFitRange(data, ti, tf, nt);
      if (doCorr)
      {
        cout << "-- uncorrelated fit..." << endl;
      }
      cout << "using linear model" << endl;
      data.setSvdTolerance(svdTol);
      data.assumeYYCorrelated(false, 0, 0);
      fit = data.fit(unCorrMin, init, mod);
      if (doCorr)
      {
        fit.print();
      }
      if (doCorr)
      {
        cout << "-- correlated fit..." << endl;
        cout << "using linear model" << endl;
        init = fit[central];
        data.assumeYYCorrelated(true, 0, 0);
        fit = data.fit(locMin, init, mod);
        fit.print();
      }
      if (!saveFileName.empty())
      {
        cout << "Saving c3 fit value to: " << saveFilePrefix << saveFileName <<".h5"<< endl;
        Io::save(fit, saveFilePrefix + saveFileName+".h5");
      }
      // plots ***************************************************************

      if (savePlot)
      {
        Plot p;

        //p << PlotRange(Axis::x, 0, nt - 1);
        p << Label("1/(g*sqrt(t))", Axis::x);
        p << Label("c_3 + f_g(g*sqrt(t))", Axis::y);
        p << Title("ag=" + g + "  N_L=" + L + "  m^2=" + m2);
        p << PlotRange(Axis::x, 0, corrFlowtimes.back() + 1);

        p << Color("rgb 'blue'") << PlotPredBand(fit.getModel(_), 0, corrFlowtimes.back() + 1);
        p << Color("rgb 'blue'") << PlotFunction(fit.getModel(), 0, corrFlowtimes.back() + 1);
        p << Color("rgb 'red'") << PlotData(data.getData());
        cout << "Saving plot to: " << savePlotPrefix << saveFileName << endl;
        p.save(savePlotPrefix + saveFileName);
      }

    return EXIT_SUCCESS;
}
