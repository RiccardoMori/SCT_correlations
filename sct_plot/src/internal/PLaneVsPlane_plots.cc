#include "internal/PLaneVsPlane_plots.hh"
#include "treeCollection.h"

#include <memory>
#include <iostream>
#include "internal/plane.hh"
#include "plane_def.h"
#include "TMath.h"
namespace sct_corr {
plotPlaneVsPlane::plotPlaneVsPlane(const s_plot_prob& plot_prob) :plot(plot_prob) {

}

bool plotPlaneVsPlane::MakeReadyForData(sct_event_buffer* outputBuffer) {
  if (!m_x) {
    return false;
    std::cout << "[plotPlaneVsPlane] first plane not set \n";
  }

  if (!m_y) {
    return false;
    std::cout << "[plotPlaneVsPlane] second plane not set \n";
  }
  m_outPutEvent = rootEvent_X_Y_hits(getName());
  m_outTree = std::make_shared<treeCollection_ouput>(m_outPutEvent, outputBuffer, getSave2disk());
  m_HitA = m_x->getHit();
  m_HitB = m_y->getHit();
  return true;
}

void plotPlaneVsPlane::pushAxis(const axis_ref* axis) {
  std::cout << "[plotPlaneVsPlane] axis not supported as input " << std::endl;

}

void plotPlaneVsPlane::pushPlane(S_plane* plane_) {
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

void plotPlaneVsPlane::processEventStart() {

}

void plotPlaneVsPlane::processEventEnd() {

}

returnProcessEvent  plotPlaneVsPlane::ProcessCurrentEvent() {
  ++m_current;
  m_size_x = 0, m_size_y = 0;
  m_outPutEvent.reset();

  bool first = true;
  processEventStart();
  while (m_y->next()) {
    ++m_size_y;
    while (m_x->next()) {
      if (first) {
        ++m_size_x;
      }
      processHit(*m_HitA, *m_HitB);
    }
    first = false;
    processHit_2(*m_HitB);
  }

  while (m_x->next()) {
    processHit_1(*m_HitA);

  }


  processEventEnd();
  return FILL_OK;
}

void plotPlaneVsPlane::fill()
{
  m_outTree->fill();
}

void plotPlaneVsPlane::pushHit(Double_t x, Double_t y) {
  pushHit(x, y, 0);
}

void plotPlaneVsPlane::pushHit(Double_t x, Double_t y, Double_t ID) {
  m_outPutEvent.push_Hit(sct_corr::plane_hit(x, y), ID);
}
void plotPlaneVsPlane::pushHit(const sct_corr::plane_hit& hit, Double_t ID) {
  m_outPutEvent.push_Hit(hit, ID);
}

Long64_t plotPlaneVsPlane::Draw(const char* options, const char* cuts /*= ""*/, const char* axis /*= "y:x"*/) {
  return Draw(S_DrawOption(options, cuts, axis));
}

Long64_t plotPlaneVsPlane::Draw(const S_DrawOption& opt) {
  return opt.Draw(m_outTree->getTTree());

}

sct_type::collectionName_t plotPlaneVsPlane::getOutputName() const {
  if (m_outTree) {
    return m_outTree->m_name;
  }
  return getName();
}


plane_distance::plane_distance(const s_plot_prob& plot_prob) : plotPlaneVsPlane(plot_prob) {

}

void plane_distance::processEventStart() {
  m__id = 0;
  m_status = 0;
}

void plane_distance::processHit(const plane_hit& p1, const plane_hit& p2) {
  if (m__id == 0) {

    if (hit_abs(p2) < 0.1) {
      if (TMath::Abs(p1.y) < 0.2) {
        m_status = 0;
      } else {
        m_status = 1;
      }
    } else {
      if (TMath::Abs(p1.y) < 0.2) {
        m_status = 10;
      } else {
        m_status = 11;

      }
    }
  }
  pushHit(p1.x, p1.y, m__id * 100 + m_status);
  ++m__id;
  //pushHit(hit_abs(p1), hit_abs(p2));
}

double plane_distance::hit_abs(const plane_hit& h) {
  return TMath::Sqrt(h.x*h.x + h.y * h.y);
}

const char* plane_distance::getType() const {
  return "plane_distance__";
}

s_plane_collection plane_distance::getOutputcollection() {
  s_plane_collection ret;
  ret.m_planes.push_back(std::make_pair(std::string("A_and_b"), plane_def(getOutputName(), sct_type::ID_t(0))));
  return ret;
}

}

sct_corr::S_plot sct_plot::plane_distance(const s_plot_prob& plot_prob) {
  return sct_corr::S_plot(new sct_corr::plane_distance(plot_prob));
}