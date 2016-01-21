#include "internal/plot_residual_with_charge.hh"

#include "treeCollection.h"
#include "sct_events/rootEvent_Alibava.hh"
#include "internal/planeAlibava_hit.hh"
#include "internal/plane_hit.hh"
#include "s_plane.h"
#include "sct_processors.h"
#include "s_plot_collection.h"


sct_corr::plot_residual_with_charge::plot_residual_with_charge(bool full_cahrge_,const s_plot_prob& plot_prob/*= ""*/):plot(plot_prob),m_full_cluster(full_cahrge_) {

}

Long64_t sct_corr::plot_residual_with_charge::Draw(const char* options, const char* cuts /*= ""*/, const char* axis /*= "y:x"*/) {
  return Draw(S_DrawOption(options, cuts, axis));
}


sct_corr::returnProcessEvent sct_corr::plot_residual_with_charge::ProcessCurrentEvent() {
  ++m_current;
  m_outputEvent.reset();

  bool first = true;
  processEventStart();
  while (m_y->next()) {
    while (m_x->next()) {
      processHit(*m_HitA, *m_HitB);
    }
    processHit_2(*m_HitB);
  }
  while (m_x->next()) {
    processHit_1(*m_HitA);

  }
  processEventEnd();
  return FILL_OK;
  
}

void sct_corr::plot_residual_with_charge::fill() {
  m_outTree->fill();
}

void sct_corr::plot_residual_with_charge::processEventStart() {
  m_max_charge_hit = Alibava_hit(0,0,0);
  m_entire_cluster = Alibava_hit(0, 0, 0);
}

void sct_corr::plot_residual_with_charge::processEventEnd() {
  if (m_max_charge_hit.charge==0)
  {
    return;
  }
  if (m_full_cluster)
  {
  pushHit(m_max_charge_hit.x, m_max_charge_hit.y, m_entire_cluster.charge, 0);
  }
  else {

  pushHit(m_max_charge_hit.x, m_max_charge_hit.y, m_max_charge_hit.charge , 0);
  }
}

void sct_corr::plot_residual_with_charge::processHit(const Alibava_hit& p1, const plane_hit& p2) {



  if (p1.charge > m_max_charge_hit.charge) {

    m_max_charge_hit.x = p1.x - p2.x;
    m_max_charge_hit.y = p1.y - p2.y;
    m_max_charge_hit.charge = p1.charge;
  }
 
}



void sct_corr::plot_residual_with_charge::processHit_1(const Alibava_hit& p1) {
  if (p1.charge<70)
  {
    return;
  }
  m_entire_cluster.x = (m_entire_cluster.x*m_entire_cluster.charge + p1.x*p1.charge) / (p1.charge + m_entire_cluster.charge);
  m_entire_cluster.y = (m_entire_cluster.y*m_entire_cluster.charge + p1.y*p1.charge) / (p1.charge + m_entire_cluster.charge);
  m_entire_cluster.charge = (p1.charge + m_entire_cluster.charge);
}

void sct_corr::plot_residual_with_charge::processHit_2(const plane_hit& p1) {

}

void sct_corr::plot_residual_with_charge::pushHit(Double_t x, Double_t y, Double_t charge, Double_t ID /*=0*/) {
  m_outputEvent.push_Hit(Alibava_hit(x, y, charge), ID);
}

s_plane_collection sct_corr::plot_residual_with_charge::getOutputcollection() {
  s_plane_collection ret;
  ret.push_back(plane_def(getOutputName(), sct_type::ID_t(0)));
  return ret;
}

bool sct_corr::plot_residual_with_charge::MakeReadyForData(sct_event_buffer* outputBuffer) {
  if (!m_x) {
    std::cout << "[plotPlaneVsPlane] first plane not set \n";
    return false;
  }

  if (!m_y) {
    std::cout << "[plotPlaneVsPlane] second plane not set \n";
    return false;
  }
  m_outputEvent = rootEvent_Alibava(getName());
  m_outTree = std::make_shared<treeCollection_ouput>(m_outputEvent, outputBuffer, getSave2disk());
  
  auto a_hit = m_x->getHit();
  auto b_hit = dynamic_cast<const Alibava_hit*>(a_hit);
  m_HitA = dynamic_cast<const sct_corr::Alibava_hit*>(m_x->getHit());
  if (!m_HitA)
  {
    std::cout << "[plotPlaneVsPlane] first plane not set \n";
    return false;
  }
  m_HitB = m_y->getHit();
  return true;
}

void sct_corr::plot_residual_with_charge::pushAxis(const axis_ref* axis) {
  std::cout << "[plotPlaneVsPlaneplot_residual_with_charge] axis not supported as input " << std::endl;
}

void sct_corr::plot_residual_with_charge::pushPlane(S_plane* plane_) {
  if (!m_x) {
    m_x = plane_->getPlane();
    return;
  }

  if (!m_y) {
    m_y = plane_->getPlane();
    return;
  }


  std::cout << "[plotPlaneVsPlane] this class only supports two input planes \n";
}

 sct_type::collectionName_t sct_corr::plot_residual_with_charge::getOutputName() const {
  if (m_outTree) {
    return m_outTree->m_name;
  }
  return getName();
}

const char* sct_corr::plot_residual_with_charge::getType() const {
  return "plot_residual_with_charge";
}

Long64_t sct_corr::plot_residual_with_charge::Draw(const S_DrawOption& opt) {
  return opt.Draw(m_outTree->getTTree());
}




sct_corr::plane_def_Alibava sct_corr::processor::residual_with_charge(
  const sct_corr::plane_def_Alibava& hits_A,
  const plane_def& hits_B, 
  bool full_cluster ,
  const s_plot_prob& plot_prob/*= "" */
  ) {

  auto pl = getPlotCollectionIfUnique(hits_A, hits_B);
  if (!pl) {
    std::cout << "[processor::moduloHitMap] referencing to different plot collection\n";
    return sct_corr::plane_def_Alibava(sct_type::collectionName_t("error"), sct_type::ID_t(0));
  }
 auto  ret= pl->addPlot(S_plot(new sct_corr::plot_residual_with_charge(full_cluster,plot_prob)), hits_A, hits_B)();
  
 sct_corr::plane_def_Alibava ret1(ret.getName(), sct_type::ID_t(ret.getID()), ret.getLayer());
 ret1.set_plot_collection(ret.get_plot());
 return ret1;
}
