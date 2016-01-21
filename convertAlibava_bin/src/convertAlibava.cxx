#include "TApplication.h"

#include "../inc/treeDUT.h"
#include <iostream>

#include <algorithm>
#include "TH2.h"
#include "tclap/CmdLine.h"
#include "tclap/ArgException.h"

using namespace std;
using namespace TCLAP;
int main(int argc, char **argv) {

  CmdLine cmd("covert Alibava files to sct files", ' ', "0.1");
  ValueArg<std::string> FileNameArg("f", "inFile", "filename", true, "", "string");
  cmd.add(FileNameArg);
  ValueArg<std::string> outFilename("o", "outFile", "output filename", true, "", "string");
  cmd.add(outFilename);
  cmd.parse(argc, argv);
 
  TFile _file(FileNameArg.getValue().c_str());
  TTree *tree = (TTree *)_file.Get("treeDUT");
  treeDUT t(tree,outFilename.getValue().c_str());
  t.Loop();

 


  return 0;
}
