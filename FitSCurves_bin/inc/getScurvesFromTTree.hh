#ifndef getScurvesFromTTree_h__
#define getScurvesFromTTree_h__
#include <vector>
#include "TTree.h"
namespace sct_corr {


  class treeCollection;
}
void setThreshold_linear_trans(double slope, double const__);
class getScurvesFromTTree {
public:
  getScurvesFromTTree(TTree* tree);
  std::vector<int> get_sort_by_threhsold() const;
  size_t size() const;
  Long_t set_pos(size_t index);
  double getX() const;
  const std::vector<double> * get_Threshold() const;
  const std::vector<double> * get_Ocuupancy() const;
  const std::vector<double> * get_Occuoancy_error() const;
  const std::vector<double> * get_total_hits() const;
private:
  double m_Threshold_out = 0;
  std::vector<double> x_pos, m_Threshold;
  size_t m_index, m_max_index;
  TTree* m_tree = nullptr;
  sct_corr::treeCollection* tree_pointer = nullptr;
  std::vector<std::vector<double>> m_Ocuupancy, m_Occuoancy_error, m_total_hits;
};
#endif // getScurvesFromTTree_h__
