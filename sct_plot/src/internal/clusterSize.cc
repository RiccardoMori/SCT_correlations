
#include "internal/plot_hit2d.hh"
#include "cluster.h"
#include "sct_plots.h"
namespace sct_corr {

class Cluster_strip :public plot_hit2d {
public:
  Cluster_strip(double clusterSize_, axis_def search_axis, const s_plot_prob& = "");
  virtual void processEventStart() override;
  virtual const char* getType() const override;
  virtual void processHit(double x, double y) override;
  virtual void processEventEnd() override;
  clusterMaker<Double_t> m_cl;
  virtual s_plane_collection getOutputcollection();
  axis_def m_ax;
};

Cluster_strip::Cluster_strip(double clusterSize_, axis_def search_axis, const s_plot_prob& plot_prob) : plot_hit2d(plot_prob), m_ax(search_axis) {
  m_cl = clusterMaker<Double_t>(clusterSize_);
}




void Cluster_strip::processEventStart() {
  m_cl.reset();
}

const char* Cluster_strip::getType() const {
  return "cluster_size_strip__";
}

void Cluster_strip::processHit(double x, double y) {
  m_cl.push_pixel(pixelHit<Double_t>(x, y));
}




void Cluster_strip::processEventEnd() {

  if (m_ax == x_axis_def) {
    for (size_t i = 0; i < m_cl.NumOfCluster(); ++i) {
      pushHit(m_cl.getCluster(i).x(), m_cl.getCluster(i).size(), 0);
    }

  }
  else if (m_ax==y_axis_def)
  {
    for (size_t i = 0; i < m_cl.NumOfCluster(); ++i) {
      pushHit(m_cl.getCluster(i).y(), m_cl.getCluster(i).size(), 0);
    }


  }

}

s_plane_collection Cluster_strip::getOutputcollection() {
  s_plane_collection ret;
  ret.m_planes.push_back(std::make_pair(std::string("clusterPos_size"), plane_def(getOutputName(), sct_type::ID_t(0))));
  return ret;
}

class clusterSize :public plot_hit2d {


public:
  clusterSize(double clusterSize_, const s_plot_prob& = "");
  virtual void processEventStart() override;
  virtual const char* getType() const override;
  virtual void processHit(double x, double y) override;
  virtual void processEventEnd() override;
  clusterMaker<Double_t> m_cl;
  virtual s_plane_collection getOutputcollection();
};



clusterSize::clusterSize(double clusterSize_, const s_plot_prob& plot_prob) : plot_hit2d(plot_prob) {
  m_cl = clusterMaker<Double_t>(clusterSize_);
}

void clusterSize::processEventStart() {
  m_cl.reset();
}

const char* clusterSize::getType() const {
  return "clusterSize__";
}

void clusterSize::processHit(double x, double y) {
  m_cl.push_pixel(pixelHit<Double_t>(x, y));
}

void clusterSize::processEventEnd() {
  for (size_t i = 0; i < m_cl.NumOfCluster(); ++i) {
    pushHit(m_cl.getCluster(i).x(), m_cl.getCluster(i).y(), 0);
    pushHit(m_cl.getCluster(i).size(), 1,1);
  }
}

s_plane_collection clusterSize::getOutputcollection() {

  s_plane_collection ret;
  ret.m_planes.push_back(std::make_pair(std::string("clusterPos"), plane_def(getOutputName(), sct_type::ID_t(0))));
  ret.m_planes.push_back(std::make_pair(std::string("clusterSize"), plane_def(getOutputName(), sct_type::ID_t(1))));
  return ret;
}
}


sct_corr::S_plot sct_plot::clustering_strip(axis_def searchAxis, Double_t Pixel_distance /*= 2*/, const s_plot_prob& plot_prob/*= ""*/) {

  return sct_corr::S_plot(new sct_corr::Cluster_strip(Pixel_distance, searchAxis, plot_prob));

}

sct_corr::S_plot sct_plot::clustering(Double_t Pixel_distance/*=2*/, const s_plot_prob& plot_prob) {


  return sct_corr::S_plot(new sct_corr::clusterSize(Pixel_distance, plot_prob));
}

