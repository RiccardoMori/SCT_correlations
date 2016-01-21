#include "TFile.h"
#include "TApplication.h"
#include "TBrowser.h"



#include "../sct_plot/inc/sct_plots.h"
#include <iostream>
#include "TCanvas.h"
#include "TMath.h"
#include "TH2.h"
#include <vector>



#include "s_plot_collection.h"
#include "s_plane.h"
#include "s_DrawOption.h"


int main(int argc, char **argv) {
  TApplication theApp("App", &argc, argv);


 
  theApp.Run();
}

