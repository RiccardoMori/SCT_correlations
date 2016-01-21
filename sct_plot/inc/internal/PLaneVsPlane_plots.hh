#ifndef PLaneVsPlane_plots_h__
#define PLaneVsPlane_plots_h__


#include "internal/plotsBase.hh"
#include <memory>
#include "sct_events/rootEventBase.hh"
#include "sct_events/rootEvent_X_Y_hits.hh"
namespace sct_corr{
  class treeCollection_ouput;


  class plotPlaneVsPlane :public plot{
  public:


    plotPlaneVsPlane(const s_plot_prob& = "");
    virtual bool MakeReadyForData(sct_event_buffer* outputBuffer)override;
    virtual void pushAxis(const axis_ref* axis) override;
    virtual void pushPlane(S_plane* axis) override;


    virtual void processEventStart();
    virtual void processHit(const plane_hit&  p1, const plane_hit&  p2) = 0;
    virtual void processHit_1(const plane_hit&  p1) {}
    virtual void processHit_2(const plane_hit&  p2) {}
    virtual void processEventEnd();
    virtual returnProcessEvent  ProcessCurrentEvent() override final;
    virtual void fill() override final;
    void pushHit(Double_t x, Double_t y);
    void pushHit(Double_t x, Double_t y, Double_t ID);
    void pushHit(const sct_corr::plane_hit& hit, Double_t ID);
    virtual Long64_t Draw(const char* options, const char* cuts = "", const char* axis = "y:x") override;
    virtual Long64_t Draw(const S_DrawOption& opt) override;
    virtual  sct_type::collectionName_t getOutputName() const;
    rootEvent_X_Y_hits m_outPutEvent;

    std::shared_ptr<treeCollection_ouput> m_outTree;
    Int_t m_current = 0;
    const plane_hit* m_HitA=nullptr, *m_HitB=nullptr;

    sct_corr::plane *m_x = nullptr, *m_y = nullptr;

    Int_t m_size_x = 0, m_size_y = 0;
  };


  class plane_distance :public plotPlaneVsPlane {
  public:
    plane_distance(const s_plot_prob& plot_prob="");
    virtual void processEventStart();
    virtual void processHit(const plane_hit&  p1, const plane_hit&  p2);
    static double hit_abs(const plane_hit& h);
    double m__id = 0, m_status = 0;
    virtual const char* getType() const override;
    virtual s_plane_collection getOutputcollection();
  };
}

#endif // PLaneVsPlane_plots_h__
