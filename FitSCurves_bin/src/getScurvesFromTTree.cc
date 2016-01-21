#include "getScurvesFromTTree.hh"
#include "treeCollection.h"


class linear_trans {
public:
  linear_trans(double slope, double constant) :m_slope(slope), m_const(constant) {}
  double operator()(double in) {
    return m_slope*in + m_const;
  }

private:
  double m_slope, m_const;
}g_threshold_trans(1, 0);

void setThreshold_linear_trans(double slope, double const__) {
  g_threshold_trans = linear_trans(slope, const__);
}


getScurvesFromTTree::getScurvesFromTTree(TTree* tree) : m_tree(tree) {
  tree_pointer = new sct_corr::treeCollection(m_tree);
  tree->SetBranchAddress("Threshold", &m_Threshold_out);
  tree_pointer->GetEntry(0);
  m_max_index = tree_pointer->m_rootBuffer.getData(x_axis_def)->size();

  m_Ocuupancy.resize(m_max_index);
  m_Occuoancy_error.resize(m_max_index);
  m_total_hits.resize(m_max_index);
  x_pos.resize(m_max_index);

  auto sorted_byThreshold = get_sort_by_threhsold();
  m_Threshold.resize(sorted_byThreshold.size());
  size_t i = 0;
  for (auto& e : sorted_byThreshold) {

    tree_pointer->GetEntry(e);
    for (size_t index = 0; index < m_max_index; ++index) {
      m_Threshold[i] = g_threshold_trans(m_Threshold_out);
      m_Ocuupancy[index].push_back(tree_pointer->m_rootBuffer.getData(Occupancy_axis_def)->at(index));
      auto Occupancy_error_axis = tree_pointer->m_rootBuffer.getData(Occupancy_error_axis_def);
      if (Occupancy_error_axis) {
        m_Occuoancy_error[index].push_back(Occupancy_error_axis->at(index));

      } else {
        m_Occuoancy_error[index].push_back(0);
      }
      m_total_hits[index].push_back(tree_pointer->m_rootBuffer.getData(NumOfEvents_axis_def)->at(index));
      x_pos[index] = tree_pointer->m_rootBuffer.getData(x_axis_def)->at(index);


    }
    ++i;
  }
}

std::vector<int> getScurvesFromTTree::get_sort_by_threhsold() const {
  std::map<double, int> dummy;

  for (Int_t i = 0; i < tree_pointer->GetEntries(); ++i) {
    tree_pointer->GetEntry(i);
    dummy[m_Threshold_out] = i;
  }
  std::vector<int> ret;
  for (auto& e : dummy) {
    ret.push_back(e.second);
  }
  return ret;
}

size_t getScurvesFromTTree::size() const {
  return m_max_index;
}

Long_t getScurvesFromTTree::set_pos(size_t index) {
  m_index = index;

  return 1;
}

double getScurvesFromTTree::getX() const {
  return x_pos.at(m_index);
}

const std::vector<double> * getScurvesFromTTree::get_Threshold() const {
  return &m_Threshold;
}

const std::vector<double> * getScurvesFromTTree::get_Ocuupancy() const {
  return &m_Ocuupancy.at(m_index);
}

const std::vector<double> * getScurvesFromTTree::get_Occuoancy_error() const {
  return &m_Occuoancy_error.at(m_index);
}

const std::vector<double> * getScurvesFromTTree::get_total_hits() const {
  return &m_total_hits.at(m_index);
}
