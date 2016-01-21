#ifndef plot2d_h__
#define plot2d_h__


#include "internal/plotsBase.hh"
#include "plane_def.h"
#include "sct_events/rootEvent_X_Y_hits.hh"
#include "sct_types.h"

namespace sct_corr{
  class treeCollection_ouput;
  class axis_ref;
  class plot2d :public plot{

  public:
    plot2d(const s_plot_prob& = "");
    virtual Long64_t Draw(const char* options, const char* cuts = "", const char* axis = "y:x") override;
    virtual Long64_t Draw(const S_DrawOption& opt) override;
    virtual void ProcessEvent() = 0;
    virtual returnProcessEvent ProcessCurrentEvent() override final;
    virtual void fill() override final;
    void pushHit(Double_t x, Double_t y);
    void pushHit(Double_t x, Double_t y, Double_t ID);
    void pushHit(const plane_hit& hit, Double_t ID);
    virtual bool MakeReadyForData(sct_event_buffer* outputBuffer) override;
    virtual void pushAxis(const axis_ref* axis);
    virtual void pushPlane(S_plane* axis);
    virtual sct_type::collectionName_t getOutputName()  const;
    Int_t get_current_event() const;
  protected:
    const axis_ref* m_x = nullptr;
    const axis_ref* m_y = nullptr;

    rootEvent_X_Y_hits m_outputEvent;

    std::shared_ptr<treeCollection_ouput> m_outTree;
    Int_t m_current = 0;
  };
}
#endif // plot2d_h__
