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



int main(int argc, char *argv[])
{
    // parse arguments /////////////////////////////////////////////////////////
    OptParser              opt;
    bool                   parsed, doCorr, doGlob, savePlot;
    string                 paramFileName, savePlotPrefix;
    double                 svdTol;
    Minimizer::Verbosity   verbosity;

    opt.addOption("", "svd", OptParser::OptType::value, true,
                  "singular value elimination threshold", "0.");
    opt.addOption("v", "verbosity", OptParser::OptType::value, true,
                  "minimizer verbosity level (0|1|2)", "0");
    opt.addOption("", "uncorr", OptParser::OptType::trigger, true,
                  "only do the uncorrelated fit");
    opt.addOption("g", "glob", OptParser::OptType::trigger, true,
                  "use global minimiser first");
    opt.addOption("s", "save-plot", OptParser::OptType::trigger, true,
                  "saves the source and .pdf", "");
    opt.addOption("", "help", OptParser::OptType::trigger, true,
                  "show this help message and exit");
    parsed = opt.parse(argc, argv);
    if (!parsed or (opt.getArgs().size() != 1) or opt.gotOption("help"))
    {
        cerr << "usage: " << argv[0] << " <parameter file> <options> " << endl;
        cerr << endl
             << "Possible options:" << endl
             << opt << endl;

        return EXIT_FAILURE;
    }

    paramFileName = opt.getArgs()[0];
    svdTol        = opt.optionValue<double>("svd");
    doCorr        = !opt.gotOption("uncorr");
    doGlob        = opt.gotOption("g");
    savePlot      = opt.gotOption("s");

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

    // read parameter file ////////////////////////////////////////////////////////
    Index               nSample, filenum;
    Latan::XmlReader    paramFile(paramFileName);
    string              fitfunction, modelprefix, model;
    vector<string>      corrFileNames, massFileNames;
    vector<double>      corrSpacings;
    vector<int>         corrVolumes;

    model             = paramFile.getFirstValue<string>("fitfunction", "model");
    modelprefix       = "model_" + model;
    savePlotPrefix    = paramFile.getFirstValue<string>("saveplotdir", "dir");
    corrFileNames     = paramFile.getAllValues<string>("filenames", "file");
    massFileNames     = paramFile.getAllValues<string>("masses", "massfile");
    corrSpacings      = paramFile.getAllValues<double>("spacings","ag");
    corrVolumes       = paramFile.getAllValues<int>("volumes","gL");
    filenum           = corrFileNames.size();
    if (corrFileNames.size()!=massFileNames.size() && massFileNames.size()!=corrSpacings.size() && corrSpacings.size()!=corrVolumes.size())
    {
      cerr << "Wrong file number in parameter file." << endl;
      return EXIT_FAILURE;
    }
    
    // load correlator files ////////////////////////////////////////////////////////

    vector<DMatSample>      corrSample(filenum);
    vector<DSample>         massSample(filenum);
    for (int i = 0; i < corrFileNames.size(); i++)
    {
      corrSample[i] =  Io::load<DMatSample>(corrFileNames[i]);
      massSample[i] =  Io::load<DSample>(massFileNames[i]);
    }
    nSample         = corrSample[0].size();
    

    // make model //////////////////////////////////////////////////////////////
    DoubleModel modGlobal;
    Index       nPar, nArg;
    if (model == "1")
    {
      nPar = 1;
      fitfunction = "return p_0;";
    } else if (model == "2")
    {
      nPar = 2;
      fitfunction = "return p_0+p_1*x_0;";
    } else if (model == "3")
    {
      nPar = 2;
      fitfunction = "return p_0+p_1*x_1;";
    } else if (model == "4")
    {
      nPar = 2;
      fitfunction = "return p_0+p_1*x_2;";
    } else {
      cerr << "model not defined" << endl;
      return EXIT_FAILURE;
    }
    nArg        = 3; 
    modGlobal   = compile(fitfunction, nArg, nPar);

    // fit /////////////////////////////////////////////////////////////////////
    XYSampleData         data(nSample);
    SampleFitResult      fit;
    DVec                 init(nPar);
    NloptMinimizer       globMin(NloptMinimizer::Algorithm::GN_CRS2_LM);
    MinuitMinimizer      locMin;
    vector<Minimizer *>  unCorrMin{&globMin, &locMin};

    //add fit dimensions and range///////////////////////////////////////////////
    data.addXDim(filenum, "m2", false);
    data.addXDim(filenum, "1/gL", true);
    data.addXDim(filenum, "ag", true);
    data.addYDim();
    for (Index s = central; s < nSample; ++s)
    {
        for (Index corridx = 0; corridx < filenum; corridx++)
        {

            data.x(corridx, 0)[s] = massSample[corridx][s];
            data.x(corridx, 1)[s] = pow(corrVolumes[corridx]*corrSpacings[corridx],-1);
            data.x(corridx, 2)[s] = corrSpacings[corridx];
            data.y(data.dataIndex(corridx, corridx, corridx), 0)[s] = corrSample[corridx][s](0);
        }
    }

    // set parameter name /////////////
    for (Index p = 0; p < nPar; p++)
    {
        modGlobal.parName().setName(p, "p_" + strFrom(p));
    }

    // set limits for minimiser //////////////
    for (Index p = 0; p < nPar; p++)
    {
        init(p) = 0.1;
        locMin.setLowLimit(p, -2.0);
        locMin.setHighLimit(p, 2.0);
        if (doGlob)
        {
            globMin.setLowLimit(p, -2.0);
            globMin.setHighLimit(p, 2.0);
        }
    }
    globMin.setPrecision(0.001);
    globMin.setMaxIteration(100000);
    globMin.setVerbosity(verbosity);
    locMin.setPrecision(0.001);
    locMin.setMaxIteration(1000000);
    locMin.setVerbosity(verbosity);

    // fit /////////////////////////////////

        for (Index corridx = 0; corridx < filenum; corridx++)
        {

            data.fitPoint(true, data.dataIndex(corridx,corridx,corridx), 0);

        }
        cout << "model prefix: " << modelprefix << endl;
        if (!doCorr)
        {
            cout << "-- uncorrelated fit..." << endl;
        }
        data.setSvdTolerance(svdTol);
        data.assumeYYCorrelated(false, 0, 0);
        data.assumeXXCorrelated(false, 0, 0);
        data.assumeXXCorrelated(false, 0, 1);
        data.assumeXXCorrelated(false, 0, 2);
        data.assumeXXCorrelated(false, 1, 1);
        data.assumeXXCorrelated(false, 1, 2);
        data.assumeXXCorrelated(false, 2, 2);
        data.assumeXYCorrelated(false, 0, 0);
        data.assumeXYCorrelated(false, 1, 0);
        data.assumeXYCorrelated(false, 2, 0);

        fit = doGlob ? data.fit(unCorrMin, init, modGlobal) : fit = data.fit(locMin, init, modGlobal);
        if (!doCorr)
        {
            fit.print();
        }

        if (doCorr)
        {
            cout << "-- correlated fit..." << endl;
            data.assumeYYCorrelated(true, 0, 0);
            fit = data.fit(locMin, init, modGlobal);
            fit.print();
        }
    

    // plots ///////////////////////////////////////////////////////////////////
    if (savePlot)
    {

        Plot plot0, plot1, plot2;
        DVec ref(3);
        XYStatData res0, res1, res2;
        ref(0) = 0.;
        ref(1) = 0.;
        ref(2) = 0.;
        //axis 0
        res0 = data.getPartialResiduals(fit, ref, 0).getData();
        plot0 << PlotRange(Axis::x, 0., 0.1);
        plot0 << Color("rgb 'blue'");
        plot0 << PlotPredBand(fit.getModel(_).bind(0, ref), 0., 0.1);
        plot0 << Color("rgb 'blue'");
        plot0 << PlotFunction(fit.getModel().bind(0, ref), 0., 0.1);
        plot0 << Color("rgb 'red'");
        plot0 << PlotData(res0, 0, 0);
        plot0 << Label("m2bar", Axis::x);
        plot0 << Label("fit", Axis::y);
        plot0 << Caption("versus m2bar");
        plot0.save(savePlotPrefix + modelprefix + "_global_fit_residual_m2bar");
        //axis 1
        res1 = data.getPartialResiduals(fit, ref, 1).getData();
        plot1 << PlotRange(Axis::x, 0., 0.2);
        plot1 << Color("rgb 'blue'");
        plot1 << PlotPredBand(fit.getModel(_).bind(1, ref), 0., 0.2);
        plot1 << Color("rgb 'blue'");
        plot1 << PlotFunction(fit.getModel().bind(1, ref), 0., 0.2);
        plot1 << Color("rgb 'red'");
        plot1 << PlotData(res1, 1, 0);
        plot1 << Label("1/gL", Axis::x);
        plot1 << Label("fit", Axis::y);
        plot1 << Caption("versus 1/gL");
        plot1.save(savePlotPrefix + modelprefix + "_global_fit_residual_gL");

        //axis 2
        res2 = data.getPartialResiduals(fit, ref, 1).getData();
        plot2 << PlotRange(Axis::x, 0., 0.4);
        plot2 << Color("rgb 'blue'");
        plot2 << PlotPredBand(fit.getModel(_).bind(2, ref), 0., 0.4);
        plot2 << Color("rgb 'blue'");
        plot2 << PlotFunction(fit.getModel().bind(2, ref), 0., 0.4);
        plot2 << Color("rgb 'red'");
        plot2 << PlotData(res2, 2, 0);
        plot2 << Label("ag", Axis::x);
        plot2 << Label("fit", Axis::y);
        plot2 << Caption("versus ag");
        plot2.save(savePlotPrefix + modelprefix + "_global_fit_residual_ag");
    }


    return EXIT_SUCCESS;
}
