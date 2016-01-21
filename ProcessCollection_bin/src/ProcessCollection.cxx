
#include "processorBase.h"
#include "TError.h"


#include <iostream>
#include <sstream>


#include "tclap/CmdLine.h"
#include "xml_helpers/xml_fileList.hh"
#include <thread>
#include "internal/exceptions.hh"


using namespace xml_util;

using namespace rapidxml;


using namespace std;
using namespace TCLAP;

void remove_root_printouts() {
  gErrorIgnoreLevel = kError;  // ignoring root printouts (replace of histograms) 
}

struct  inParam {
  int argc;
  char **argv;
};


int asyncMain(void *arg) {

  inParam* para = static_cast<inParam *>(arg);
  int argc = para->argc;
  char **argv = para->argv;

  try {

    CmdLine cmd("ProcessCollection", ' ', "0.1");
    
    ValueArg<std::string> FileNameArg("i", "inFile", "xml filename", true, "", "string");
    cmd.add(FileNameArg);

    ValueArg<std::string>  inPath("p", "inPath", "path to the root files", true, "", "string");
    cmd.add(inPath);
    
    ValueArg<std::string>  outpath("o", "outPath", "output path", false, ".", "string");
    cmd.add(outpath);
    ValueArg<std::string>  processor_type("s", "Processor", "which processor to use Standard or Modulo", false, "Standard", "string");
    cmd.add(processor_type);
    TCLAP::SwitchArg regprocessor_type("r", "registeredProcessors", "shows the processors registered");
    cmd.add(regprocessor_type);
#ifdef _DEBUG
    cmd.setExceptionHandling(false);
#endif // _DEBUG

    cmd.parse(argc, argv);
    if (regprocessor_type.getValue())  {

      display_registered_processors(std::cout);
      return 0;
    }

    auto p = create_processor(processor_type.getValue());

    p->Add_XML_RunList(FileNameArg.getValue(), inPath.getValue(), outpath.getValue());
    p->setPrintout(true);
    p->process();
  } catch (ArgException &e)  // catch any exceptions
  {
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
    return -1;
  } catch (...) {
    display_registered_processors(std::cout);
    return sct_corr::handleExceptions();
  }

  return 0;
}
int main(int argc, char **argv) {

  remove_root_printouts();
  
  




  inParam para;
  para.argc = argc;
  para.argv = argv;


  asyncMain(&para);
  



  return 0;

}

