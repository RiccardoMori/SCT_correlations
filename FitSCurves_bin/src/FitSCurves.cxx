#include "TFile.h"
#include "TApplication.h"
#include "TBrowser.h"




#include <iostream>
#include "TCanvas.h"
#include "TMath.h"
#include "TH2.h"
#include <vector>
#include "TProfile.h"
#include "TFitResult.h"
#include "TF1.h"
#include "TBenchmark.h"
#include "plane_def.h"
#include "s_plot_collection.h"
#include "s_DrawOption.h"
#include "s_cuts.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "processorBase.h"

#include "landgausFit.h"
#include <algorithm>
#include <fstream>
#include "tclap/CmdLine.h"
#include "tclap/ArgException.h"
#include <string>
#include "sct_events/rootEventBase.hh"
#include "sct_events/rootEventRunOutput.hh"
#include "treeCollection.h"
#include "TTree.h"
#include "SCT_helpers.h"
#include "scurve_fit_collection.hh"
#include "internal/exceptions.hh"

std::string trim(const std::string & s) {
  static const std::string spaces = " \t\n\r\v";
  size_t b = s.find_first_not_of(spaces);
  size_t e = s.find_last_not_of(spaces);
  if (b == std::string::npos || e == std::string::npos) {
    return "";
  }
  return std::string(s, b, e - b + 1);
}



std::vector<std::string> split(const std::string & str, const std::string & delim, bool dotrim) {
  std::string s(str);
  std::vector<std::string> result;
  if (str == "") return result;
  size_t i;
  while ((i = s.find_first_of(delim)) != std::string::npos) {
    result.push_back(dotrim ? trim(s.substr(0, i)) : s.substr(0, i));
    s = s.substr(i + 1);
  }
  result.push_back(s);
  return result;
}
std::vector<std::string> split(const std::string & str, const std::string & delim) {
  return split(str, delim, false);
}
std::vector<unsigned> parsenumbers(const std::string & s) {
  std::vector<unsigned> result;
  std::vector<std::string> ind = split(s, ",");
  for (auto&e : ind) {
    result.push_back(atoi(e.c_str()));
  }
  return result;
}



using namespace std;
using namespace TCLAP;
int main(int argc, char **argv) {
#ifdef _DEBUG
  int argc_asdasddsa = 1;
  char **argv_sadsda = new char*[1];
  argv_sadsda[0] = "dummy";
  TApplication theApp("App", &argc_asdasddsa, argv_sadsda);

#endif // _DEBUG


  try {

    CmdLine cmd("SCurve Fit", ' ', "0.1");
    ValueArg<std::string> FileNameArg("f", "inFile", "filename", true, "", "string");
    cmd.add(FileNameArg);
    ValueArg<std::string> outFilename("o", "outFile", "output filename", false, "", "string");
    cmd.add(outFilename);
    ValueArg<double> MPV_arg("m", "mpv", "start value for the landau gauss fit", false, 80, "double");
    cmd.add(MPV_arg);
    ValueArg<double> N_arg("n", "NumberOfElements", "The number of Events needed for the individual strips", false, 1000, "double");
    cmd.add(N_arg);
    ValueArg<double> MPV_arg_lower_cut("c", "mpv_lower_cut", "lower cut for mpv value", false, 280, "double");
    cmd.add(MPV_arg_lower_cut);

    ValueArg<std::string> ignore_elements("i", "ignorsThreshold", "ignores threshold by index", false, "", "string");
    cmd.add(ignore_elements);
    cmd.parse(argc, argv);
    auto index = parsenumbers(ignore_elements.getValue());


    std::string inFile = FileNameArg.getValue();

    std::string outfile;
    if (!outFilename.isSet()) {
      outfile = inFile.substr(0, inFile.size() - 5);
      outfile += "_scurve.txt";
    } else {
      outfile = outFilename.getValue();
    }
    std::cout << outfile << std::endl;

    TFile _file0(inFile.c_str());
    if (!_file0.IsOpen()) {
      std::cout << "Unable to Open the File    \n ";
      return -3;
    }

    TTree* noise = (TTree*)_file0.Get("hitmap");
    TTree* Hits = (TTree*)_file0.Get("out");

    if (!noise&&!Hits) {
      std::cout << "TTree not found " << std::endl;
      return -2;

    }
    scurve_fit_collection p(outfile.c_str());


    if (noise) {
      p.setNoiseRun(noise);
    } else if (Hits) {
      p.setBeamRuns(Hits);
    }
    p.set_mpv_start_low_cut(MPV_arg_lower_cut.getValue());
    p.set_number_of_min_events(N_arg.getValue());
    p.setStartMPV(MPV_arg.getValue());
    gErrorIgnoreLevel = kError;  // ignoring root printouts (replace of histograms) 
    p.processTotal("total_efficiency:Threshold:error_efficiency");
    p.processStrip("Occupancy:Threshold:Occupancy_error", index);
#ifdef _DEBUG
    //   theApp.Run();
#endif
    return 0;
  } catch (ArgException &e)  // catch any exceptions
  {

    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
    return -1;
  }

  return 0;
}
