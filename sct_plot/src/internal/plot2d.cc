#include <memory>
#include "treeCollection.h"
#include "sct_plots.h"
#include "internal/plot2d.hh"
#include "internal/plane.hh"
#include "sct_plots.h"
#include <iostream>

namespace sct_corr{

  plot2d::plot2d(const s_plot_prob& plot_prob) :plot(plot_prob)
  {

  }

  Long64_t plot2d::Draw(const char* options, const char* cuts /*= ""*/, const char* axis /*= "y:x"*/)
  {
    return Draw(S_DrawOption(options, cuts, axis));
  }

  Long64_t plot2d::Draw(const S_DrawOption& opt)
  {
    return opt.Draw(m_outTree->getTTree());
    
  }

  sct_corr::returnProcessEvent plot2d::ProcessCurrentEvent() {
    ++m_current;
    m_outputEvent.reset();


    ProcessEvent();
    return FILL_OK;
  }

  void plot2d::fill()
{
    m_outTree->fill();
  }

  void plot2d::pushHit(Double_t x, Double_t y)
  {
    
    pushHit(x, y, 0);
  }

  void plot2d::pushHit(Double_t x, Double_t y, Double_t ID)
  {
    m_outputEvent.push_Hit(sct_corr::plane_hit(x, y), ID);
  }

  void plot2d::pushHit(const plane_hit& hit, Double_t ID) {
    m_outputEvent.push_Hit(hit, ID);
  }

  bool plot2d::MakeReadyForData(sct_event_buffer* outputBuffer)
  {
    if (m_x &&m_y)
    {
      m_outputEvent = rootEvent_X_Y_hits(getName());
      m_outTree = std::make_shared<treeCollection_ouput>(m_outputEvent, outputBuffer, getSave2disk());
      return true;
    }

    return false;
  }

  void plot2d::pushAxis(const axis_ref* axis)
  {
    if (!m_x)
    {
      m_x = axis;
      return;
    }

    if (!m_y)
    {
      m_y = axis;
      return;
    }


    std::cout << "[plot2d::pushAxis(axis_ref* axis)] only two axis are supported\n";
  }

  void plot2d::pushPlane(S_plane* plane_)
  {
    if (!m_x&&!m_y)
    {
      m_x = plane_->getPlane()->getAxis(x_axis_def);
      m_y = plane_->getPlane()->getAxis(y_axis_def);
      return;
    }
    std::cout << "[plot2d::pushPlane(S_plane* plane)] only one plane is supported \n";
  }

  sct_type::collectionName_t plot2d::getOutputName() const
  {
    return m_outTree->m_name;
  }

  Int_t plot2d::get_current_event() const {
    return m_current;
  }

}