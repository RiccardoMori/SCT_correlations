#ifndef plotsBase_h__
#define plotsBase_h__


#include "plot_save_option_def.h"
#include "s_plot_prob.h"
#include "s_plot_collection.h"
#include "Rtypes.h"
#include "s_DrawOption.h"
#include "s_plane.h"
#include "plane_def.h"



namespace sct_corr{
enum returnProcessEvent :int{
  FILL_OK,
  FILL_DONE,
  FILL_SKIP
};
class sct_event_buffer;
  class plot{
  public:

    plot(const s_plot_prob& = "");
    virtual const char* getType() const = 0;
    virtual bool MakeReadyForData(sct_event_buffer* outputBuffer) = 0;
    virtual void pushAxis(const axis_ref* axis) = 0;
    virtual void pushPlane(S_plane* axis) = 0;
    virtual returnProcessEvent  ProcessCurrentEvent() = 0;
    virtual void fill() = 0;
    virtual Long64_t Draw(const char* options, const char* cuts = "", const char* axis = "y:x") = 0;
    virtual Long64_t Draw(const S_DrawOption&) = 0;
    virtual void setParameter(const char* tag, const char * value);
    virtual s_plane_collection getOutputcollection();

    virtual sct_type::collectionName_t getOutputName()  const = 0;
    sct_type::collectionName_t getName() const;
    bool getSave2disk() const;
  private:
    s_plot_prob m_prob;
    mutable  sct_type::collectionName_t m_name;
  };
  using S_plot = std::shared_ptr<sct_corr::plot>;
}

#endif // plotsBase_h__
