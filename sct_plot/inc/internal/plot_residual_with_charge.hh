#ifndef plot_residual_with_charge_h__
#define plot_residual_with_charge_h__


#include "internal/plotsBase.hh"
#include "plane_def.h"
#include "sct_events/rootEvent_X_Y_hits.hh"
#include "sct_events/rootEvent_Alibava.hh"
#include <memory>
#include "planeAlibava_hit.hh"
#include "plane_hit.hh"
#include "plane.hh"

namespace sct_corr {
class treeCollection_ouput;

class plot_residual_with_charge :public plot {

public:
  plot_residual_with_charge(bool full_charge,const s_plot_prob& = "");
  virtual Long64_t Draw(const char* options, const char* cuts = "", const char* axis = "y:x") override;
  virtual Long64_t Draw(const S_DrawOption& opt) override;

  virtual returnProcessEvent ProcessCurrentEvent() override final;
  virtual void  fill() override final;
  void processEventStart();
  void processEventEnd();
  void processHit(const Alibava_hit& p1, const plane_hit& p2);
  void processHit_1(const Alibava_hit& p1);
  void processHit_2(const plane_hit& p1);

  void pushHit(Double_t x, Double_t y,Double_t charge,Double_t ID =0);
  virtual s_plane_collection getOutputcollection();
  virtual bool MakeReadyForData(sct_event_buffer* outputBuffer) override;
  virtual void pushAxis(const axis_ref* axis);
  virtual void pushPlane(S_plane* axis);
  virtual sct_type::collectionName_t getOutputName()  const override;
  virtual const char* getType() const override;
protected:
  const Alibava_hit* m_HitA = nullptr;
  const plane_hit *m_HitB = nullptr;
  Alibava_hit m_max_charge_hit,m_entire_cluster;
  sct_corr::plane *m_x = nullptr, *m_y = nullptr;

  rootEvent_Alibava m_outputEvent;
  bool m_full_cluster = true;
  std::shared_ptr<treeCollection_ouput> m_outTree;
  Int_t m_current = 0;
};
}


#endif // plot_residual_with_charge_h__
