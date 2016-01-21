#include "internal/plotsBase.hh"
#include <iostream>
#include "treeCollection.h"
#include <vector>
#include "sct_events/rootEvent_X_Y_hits.hh"
#include "internal/plane.hh"

namespace sct_corr{
  class plane_merger :public plot{
  public:
    plane_merger( const s_plot_prob& = "");
    virtual const char* getType() const  override;
    virtual bool MakeReadyForData(sct_event_buffer* outputBuffer) override;
    virtual void pushAxis(const axis_ref* axis) override;
    virtual void pushPlane(S_plane* axis) override;
    virtual returnProcessEvent ProcessCurrentEvent() override;
    virtual void fill() override;
    virtual Long64_t Draw(const char* options, const char* cuts = "", const char* axis = "y:x")  override;
    virtual Long64_t Draw(const S_DrawOption&) ;
    virtual s_plane_collection getOutputcollection();

    virtual sct_type::collectionName_t getOutputName()  const override;
    void processEvent(double x, double y, double id_);
    void processPlane(sct_corr::plane* plane_, double id_);
  private:
    std::vector<S_plane*> m_planes;
    rootEvent_X_Y_hits m_outputEvent;

    std::shared_ptr<treeCollection_ouput> m_outTree;
    Int_t m_current = 0;
  };

  plane_merger::plane_merger(const s_plot_prob& plot_prob/*= ""*/) :plot(plot_prob)
  {

  }

  const char* plane_merger::getType() const
  {
    return "plane_merger";
  }

  bool plane_merger::MakeReadyForData(sct_event_buffer* outputBuffer)
  {
    m_current = 0;
    if (m_planes.empty())
    {
      return false;
    }
    m_outputEvent = rootEvent_X_Y_hits(getName());
    m_outTree = std::make_shared<treeCollection_ouput>(m_outputEvent, outputBuffer, getSave2disk());
    return true;
  }

  void plane_merger::pushAxis(const axis_ref* axis)
  {
    std::cout << "[plane_merger] axis are not supported " << std::endl; 
  }

  void plane_merger::pushPlane(S_plane* plane_)
  {
    if (!plane_)
    {
      std::cout << "[plane_merger] plane is nullptr\n";
    }

    m_planes.push_back(plane_);
  }

  sct_corr::returnProcessEvent plane_merger::ProcessCurrentEvent() {

    ++m_current;
    m_outputEvent.reset();
    double i = 0;
    for (auto& e : m_planes) {
      processPlane(e->getPlane(), i++);

    }

    return FILL_OK;
  }

  void plane_merger::fill(){

  }

  void plane_merger::processEvent(double x, double y, double id_)
  {
    m_outputEvent.push_Hit(sct_corr::plane_hit(x, y), id_);
  }

  void plane_merger::processPlane(sct_corr::plane* plane_, double id_)
  {
    while (plane_->next())
    {
      auto h = plane_->getHit();
      processEvent(h->x, h->y, id_);
    }
  }


  Long64_t plane_merger::Draw(const S_DrawOption& opt)
  {
    return opt.Draw(m_outTree->getTTree());

  }

  Long64_t plane_merger::Draw(const char* options, const char* cuts /*= ""*/, const char* axis /*= "y:x"*/)
  {
    return Draw(S_DrawOption(options, cuts, axis));
  }

  s_plane_collection plane_merger::getOutputcollection()
  {
    s_plane_collection ret;
    double i = 0;
    for (auto& e : m_planes)
    {
      ret = ret + plane_def(e->getName(), sct_type::ID_t(i++));
    }
    return ret;
  }

  sct_type::collectionName_t plane_merger::getOutputName() const
  {
    if (m_outTree)
    {
      return m_outTree->m_name;
    }
    return getName();
  }
}

