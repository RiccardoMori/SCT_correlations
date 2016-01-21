
#include "s_plot_collection.h"
#include "sct_plots_internal.h"
#include "treeCollection.h"
#include <iostream>
#include "TCanvas.h"

#include "sct_event_buffer.h"

#include "internal/plotsBase.hh"
#include "s_DrawOption.h"
#include "plane_def.h"
#include "TTree.h"
#include "plane_def.h"
#include <ctime>
#include "internal/plot_collection_impl.hh"
#include "internal/exceptions.hh"
#include "sct_events/rootEventBase.hh"


namespace sct_corr {






sct_corr::plot_collection::plot_collection()  {

}
plot_collection::~plot_collection() {

}




plot_collection_impl::plot_collection_impl() : m_eventBuffer(std::make_shared<sct_corr::sct_event_buffer>()) {

}

plot_collection_impl::~plot_collection_impl() {

}

void plot_collection_impl::addFile(TFile* file) {
  m_file.push_back(file);
}

void plot_collection_impl::setOutputFile(TFile* file) {
  if (file) {
    addFile(file);

  }
  m_eventBuffer->setOutputFile(file);
}

void plot_collection_impl::reset() {
  m_eventBuffer->reset();
  m_planes.clear();
  m_plots.clear();
  m_trees.clear();
}

s_plane_collection plot_collection_impl::addPlot(S_plot plot_def, const S_Axis& x_axis, const S_Axis& y_axis) {
  plot_def->pushAxis(getAxis_ref(x_axis));
  plot_def->pushAxis(getAxis_ref(y_axis));
  return addPlot_internal(std::move(plot_def));
}


s_plane_collection plot_collection_impl::addPlot_internal(S_plot plot_def) {
  std::string coll_name = Un_necessary_CONVERSION(plot_def->getName());
  m_plots.push_back(std::make_pair(coll_name, std::move(plot_def)));
  if (!m_plots.back().second->MakeReadyForData(m_eventBuffer.get())) {
    SCT_THROW(std::string("unable to create plot ") + std::string(plot_def->getType()) + std::string(":") + necessary_CONVERSION(plot_def->getName()));
  }
  auto ret = m_plots.back().second->getOutputcollection();
  ret.set_plot_collection(m_self);
  return ret;
}

s_plane_collection plot_collection_impl::addPlot(S_plot plot_def, const plane_def& p1) {
  return addPlot(std::move(plot_def), s_plane_collection(p1));
}

s_plane_collection plot_collection_impl::addPlot(S_plot plot_def, const plane_def& p1, const plane_def & p2) {
  return addPlot(std::move(plot_def), p1 + p2);
}


s_plane_collection plot_collection_impl::addPlot(S_plot plot_def, const s_plane_collection& p1) {
  for (auto &e : p1.m_planes) {
    auto p1_pointer = pushPlane(e.second);
    if (!p1_pointer) {
      SCT_THROW("planes not set correctly");
    }
    plot_def->pushPlane(p1_pointer);
  }

  return addPlot_internal(plot_def);
}




Long64_t plot_collection_impl::Draw(const char* name, const S_DrawOption& option) {
  for (auto &e : m_plots) {
    if (e.first == name) {
      return  e.second->Draw(option);
    }
  }
  auto tree_ = getTTree(sct_type::collectionName_t(name));
  if (tree_) {
    return option.Draw(tree_);
  }
  return -1;
}



Long64_t plot_collection_impl::Draw(const sct_corr::plane_def& name, const S_DrawOption& option) {

  S_DrawOption local(option);

  TCut dummy = ("ID == " + std::to_string( necessary_CONVERSION(name.getID()))).c_str(); //external library 
  local.cut_add(dummy);

  return Draw(Un_necessary_CONVERSION(name.getName()).c_str(), local);

}
Long64_t plot_collection_impl::getMaxEntriesFromTree(Long64_t last) {
  bool external_last_event = true;
  if (last == -1) {
    last = kMaxInt;
    external_last_event = false;
  }
  for (auto & e : m_trees) {
    auto l = e.second->GetEntries();
    if (last < kMaxInt && l != last) {
      if (l < kMaxInt && external_last_event == false) {
        std::cout << "trees must have the same length" << std::endl;
      }


    }
    if (l < last) {
      last = l;
    }
  }

  return last;
}

void plot_collection_impl::loop(Long64_t last /*= -1*/, Long64_t start /*= 0*/) {


  last = getMaxEntriesFromTree(last);
  for (Long64_t i = start; i < last; ++i) {
#ifdef _DEBUG

    if (i % 1000 == 0) {
      std::cout << "processing event: " << i << " of " << last << "\n";
    }
#endif // _DEBUG

    for (auto& e : m_trees) {
      e.second->GetEntry(i);
    }
    returnProcessEvent ret = FILL_OK;
    for (auto& current_plot : m_plots) {
      auto current_ret = current_plot.second->ProcessCurrentEvent();
      if (current_ret == FILL_DONE) {
        std::cout << "run terminated by plot: " << current_plot.first << std::endl;
        return;
      } else if (current_ret == FILL_SKIP) {
        ret = FILL_SKIP;
      }
    }

    if (ret==FILL_OK) {
      for (auto& current_plot : m_plots) {
        current_plot.second->fill();
      }
    }
  }
}

bool plot_collection_impl::collectionExist(const sct_type::collectionName_t& name) const {


  if (getTTree(name)) {
    return true;
  }
  return false;
}

void plot_collection_impl::set_self_weak_pointer(std::weak_ptr<sct_corr::plot_collection> self_) {
  m_self = self_;
}

const sct_corr::axis_ref* plot_collection_impl::getAxis_ref(const S_Axis & axis) {

  return getPlane(axis.m_planeID, getCollection(axis.m_collectionName))->getAxis(axis.m_axis);
}

sct_corr::rootEventBase* plot_collection_impl::getCollection(const sct_type::collectionName_t& name) {
  for (auto&e : m_trees) {
    if (Un_necessary_CONVERSION(e.first) == Un_necessary_CONVERSION(name)) { // mising opperator == 
      return &(e.second->m_rootBuffer);
    }
  }

  if (m_eventBuffer->IsCollection(name)) {
    
    return m_eventBuffer->get(name);
  }

  if (m_file.empty()) {
    SCT_THROW("file empty");
  }
  TTree *collection = getTTree(name);



  if (!collection) {
    SCT_THROW("collection not found! Collection name: \"" + necessary_CONVERSION(name) + "\"");
  }


  sct_corr::treeCollection* tree_pointer = new sct_corr::treeCollection(collection);

  m_trees.push_back(std::make_pair(name, tree_pointer));
  return &(tree_pointer->m_rootBuffer);

}

TTree* plot_collection_impl::getTTree(const sct_type::collectionName_t& name) const {
  TTree *collection = NULL;

  for (auto& e : m_file) {
    e->GetObject(necessary_CONVERSION(name).c_str(), collection); ////external library 
    if (collection) {
      return collection;
    }
  }
  SCT_THROW("unable to find colllection. Collection name: " + necessary_CONVERSION(name));
  return nullptr;
}

S_plane* plot_collection_impl::getPlane(const sct_type::ID_t&  ID, sct_corr::rootEventBase* coll) {
  if (!coll) {
    SCT_THROW("unable to get plane, Root Event base is empty.");
    return nullptr;
  }

  m_planes.push_back(std::make_shared<S_plane>(plane_def(coll->getName(), ID), coll));
  return m_planes.back().get();
}

S_plane* plot_collection_impl::pushPlane(const plane_def& pl) {

  auto pl_pointer = getCollection(pl.getName());
  if (!pl_pointer) {
    SCT_THROW("Unable to push plane. collection not found. collection name: " + necessary_CONVERSION(pl.getName()));
    return nullptr;
  }


  m_planes.push_back(std::make_shared<S_plane>(pl, pl_pointer));

  return m_planes.back().get();
}

std::shared_ptr<sct_corr::plot_collection> create_plot_collection() {
  auto  ret = std::make_shared<sct_corr::plot_collection_impl>();

  ret->set_self_weak_pointer(ret);

  return std::dynamic_pointer_cast<sct_corr::plot_collection>(ret);
}





}