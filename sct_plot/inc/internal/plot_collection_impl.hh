#ifndef plot_collection_impl_hh___
#define plot_collection_impl_hh___
#include <memory>
#include <map>

#include "internal/platform.hh"
#include <vector>
#include "sct_types.h"

namespace sct_corr {
class  DllExport plot_collection_impl :public plot_collection {
public:
  plot_collection_impl();
  virtual ~plot_collection_impl();

  virtual void addFile(TFile* file) override;
  virtual void setOutputFile(TFile* file) override;
  virtual void reset() override;
  virtual s_plane_collection addPlot(S_plot plot_def, const S_Axis& x_axis, const S_Axis& y_axis) override;

  virtual s_plane_collection addPlot(S_plot plot_def, const sct_corr::plane_def& p1, const sct_corr::plane_def & p2) override;
  virtual s_plane_collection addPlot(S_plot plot_def, const  sct_corr::plane_def& p1) override;
  virtual s_plane_collection addPlot(S_plot plot_def, const  s_plane_collection& p1) override;

  virtual void loop(Long64_t last = -1, Long64_t start = 0) override;


  virtual Long64_t Draw(const char* name, const S_DrawOption& option) override;
  virtual Long64_t Draw(const sct_corr::plane_def& name, const S_DrawOption& option) override;




  bool collectionExist(const sct_type::collectionName_t& name)  const override;

  void set_self_weak_pointer(std::weak_ptr<plot_collection> self_);
private:
  Long64_t getMaxEntriesFromTree(Long64_t last);
  s_plane_collection addPlot_internal(S_plot plot_def);
  const sct_corr::axis_ref* getAxis_ref(const S_Axis & axis);
  sct_corr::rootEventBase* getCollection(const sct_type::collectionName_t& name);
  TTree*                    getTTree(const sct_type::collectionName_t& name) const;
  S_plane* getPlane(const sct_type::ID_t& id_, sct_corr::rootEventBase* coll);
  S_plane* pushPlane(const sct_corr::plane_def& pl);

  std::shared_ptr<sct_corr::sct_event_buffer> m_eventBuffer;

  std::vector<std::shared_ptr<S_plane>> m_planes;
  std::vector<std::pair<std::string, S_plot>> m_plots;
  std::vector< std::pair<sct_type::collectionName_t, sct_corr::treeCollection*>> m_trees;
  std::vector<TFile*> m_file;
  std::weak_ptr<plot_collection> m_self;

};
}

#endif