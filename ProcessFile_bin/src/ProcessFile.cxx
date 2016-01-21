
#include "processorBase.h"

#include "TError.h"

#include "TApplication.h"
#include "TBrowser.h"
#include "s_cuts.h"
#include "TF1.h"
#include "TProfile.h"
#include "xml_helpers/xml_fileList.hh"

#include "tclap/CmdLine.h"
#include <iostream>
#include "TCanvas.h"
#include <thread>

#include "TSystem.h"
#include "internal/exceptions.hh"

using namespace xml_util;
using namespace TCLAP;
using namespace std;







std::string trim(const std::string &s) {
  static const std::string spaces = " \t\n\r\v";
  size_t b = s.find_first_not_of(spaces);
  size_t e = s.find_last_not_of(spaces);
  if (b == std::string::npos || e == std::string::npos) {
    return "";
  }
  return std::string(s, b, e - b + 1);
}

std::vector<std::string> split(const std::string &str,
  const std::string &delim, bool dotrim) {
  std::string s(str);
  std::vector<std::string> result;
  if (str == "")
    return result;
  size_t i;
  while ((i = s.find_first_of(delim)) != std::string::npos) {
    result.push_back(dotrim ? trim(s.substr(0, i)) : s.substr(0, i));
    s = s.substr(i + 1);
  }
  result.push_back(s);
  return result;
}
std::unique_ptr<xmlImputFiles::MinMaxRange<double>> make_range(const std::string& r_) {
  auto r = split(r_, ":", true);

  if (r.size() == 3) {
    return std::unique_ptr<xmlImputFiles::MinMaxRange<double>>(new xmlImputFiles::MinMaxRange<double>(stof(r[0]), stof(r[1]), stof(r[2])));
  }
  else if (r.size() == 2) {
    return  std::unique_ptr<xmlImputFiles::MinMaxRange<double>>( new xmlImputFiles::MinMaxRange<double>(stof(r[0]), stof(r[1])));
  }

  std::cout << "error processing range " << r_ << std::endl;
  return   std::unique_ptr<xmlImputFiles::MinMaxRange<double>>(); 

}


struct  inParam {
  int argc;
  char **argv;
};
int asyncMain(void *arg) {
  try {
    inParam* para = static_cast<inParam *>(arg);
    int argc = para->argc;
    char **argv = para->argv;

    CmdLine cmd("ProcessFile", ' ', "0.1");
    ValueArg<std::string> FileNameArg("i", "inFile", "xml filename", true, "", "string");
    cmd.add(FileNameArg);
    ValueArg<std::string>  inPath("p", "inPath", "path to the root files", true, "", "string");
    cmd.add(inPath);
    ValueArg<std::string>  output_path("o", "outPath", "output path", false, "dummy.root", "string");
    cmd.add(output_path);
    ValueArg<std::string>  processor_type("s", "Processor", "which processor to use Standard or Modulo", false, "Standard", "string");
    cmd.add(processor_type);
    ValueArg<int> element("e", "element", "element of interest  in the XML file", true, 1, "int");
    cmd.add(element);

    ValueArg<std::string> residualRange("c", "residualCut", "range for the residual -r -10:10 => min =-10 , max = 10  ", false, "-10:10", "range");

    cmd.add(residualRange);


    cmd.parse(argc, argv);  //terminates on error

    TFile * __file1 = new TFile(output_path.getValue().c_str(), "recreate");


    auto p= create_processor(processor_type.getValue());

    p->setPrintout(true);

    gErrorIgnoreLevel = kError;  // ignoring root printouts (replace of histograms) 



    p->Add_XML_RunList(FileNameArg.getValue(), inPath.getValue(), ".", element.getValue());

#ifdef _DEBUG
    TApplication theApp("App", &argc, argv);
#endif // _DEBUG

    p->process();

    auto  r = make_range(residualRange.getValue());
    p->saveHistograms(__file1, r.get());

    __file1->Write();
    new TBrowser();
#ifdef _DEBUG

    theApp.Run();
#endif // _DEBUG

  }
  catch (...) {
    display_registered_processors(std::cout);
    return sct_corr::handleExceptions();
  }
  return 0;
}

int main(int argc, char **argv) {

  try {
    inParam para;
    para.argc = argc;
    para.argv = argv;

    asyncMain(&para);
  } catch (...) {
    std::cout << "error on exit" << std::endl;
  }

  

  return 0;


}
