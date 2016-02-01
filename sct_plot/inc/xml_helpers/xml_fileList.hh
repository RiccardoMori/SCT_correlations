#ifndef xml_fileList_h__
#define xml_fileList_h__


#include "rapidxml.hpp"
#include "rapidxml_utils.hpp"
#include "rapidxml.hpp"
#include "xml_helpers/xml_util.hh"

namespace xmlImputFiles {



  using xml_n = rapidxml::xml_node < char > ;
  template<typename T>
  class MinMaxRange {
  public:
    MinMaxRange(T min_, T max_) :m_min(min_), m_max(max_),m_stepSize(1) {}
    MinMaxRange(T min_, T step_, T max_) :m_min(min_), m_max(max_), m_stepSize(step_) {}

    T getMin() const {
      return m_min;
    }
    T getMax() const {
      return m_max;
    }
    T getStep() const {
      return m_stepSize;
    }
  private:
    T m_min, m_max,m_stepSize;

  };
  using intRange = MinMaxRange < int > ;
  class xml_globalConfig {
  public:
    xml_globalConfig(xml_n* xIn) :m_node(xIn) {
      // std::cout << xIn->name << std::endl;
    }
    xml_globalConfig() {}
    std::string CollectionName() const {
      return xml_util::getAtribute_string(m_node->first_node("CollectionName"), "value");
    }
    int NumberOfBins() const {
      return xml_util::getAtribute<int>(m_node->first_node("NumberOfBins"), "value");
    }
    int NumberOfStrips() const {
      return xml_util::getAtribute<int>(m_node->first_node("NumberOfStrips"), "value");
    }
    const S_Cut& cut() const {
      auto ymin_ = xml_util::getAtribute<Double_t>(m_node->first_node("YCut"), "min");
      auto ymax_ = xml_util::getAtribute<Double_t>(m_node->first_node("YCut"), "max");
      auto xmin_ = xml_util::getAtribute<Double_t>(m_node->first_node("XCut"), "min");
      auto xmax_ = xml_util::getAtribute<Double_t>(m_node->first_node("XCut"), "max");
      m_cut = S_YCut(ymin_, ymax_) + S_XCut(xmin_, xmax_);
      return  m_cut;
    }
    int Device() const {
      return xml_util::getAtribute<int>(m_node->first_node("Device"), "value");
    }

    intRange AvtiveStrips() const {
      auto min_ = xml_util::getAtribute<int>(m_node->first_node("AvtiveStrips"), "min");
      auto max_ = xml_util::getAtribute<int>(m_node->first_node("AvtiveStrips"), "max");
      return intRange(min_, max_);
    }
    double Rotation()const {
      return xml_util::getAtribute<double>(m_node->first_node("Rotation"), "value");
    }
    std::string Position_name()const {
      return xml_util::getAtribute_string(m_node->first_node("Position"), "name");
    }
    double Position_value()const {
      return xml_util::getAtribute<double>(m_node->first_node("Position"), "value");
    }
    double residual_cut() const {
      return xml_util::getAtribute<double>(m_node->first_node("residual_cut"), "value");
    }
    std::string gearFile() const {
      return xml_util::getAtribute_string(m_node->first_node("gearFile"), "name");
    }
  private:
    xml_n* m_node = nullptr;
    mutable  S_CutCoollection m_cut = S_YCut(-10000000, 100000000000)+S_XCut(-10000,10000);
  };
  class xml_file {
  public:
    xml_file(xml_n* xIn) :
      m_node(xIn) {

    }
    static const char* NodeName() {
      return "file";
    }
    xml_file() {}
    std::string name() const {
      return xml_util::getAtribute_string(m_node->first_node("name"), "value");
    }
    double threshold() const {
      return xml_util::getAtribute<double>(m_node->first_node("threshold"), "value");
    }
    double HV() const {
      return xml_util::getAtribute<double>(m_node->first_node("HV"), "value");
    }
    int runNumber() const {
      return xml_util::getAtribute<int>(m_node->first_node("runNumber"), "value");
    }
  private:
    xml_n* m_node = nullptr;
  };
  class XML_imput_file {
  public:
    XML_imput_file(const char* name) :m_file(name) {
      m_doc.parse<0>(m_file.data());    // 0 means default parse flags
      //      auto p = m_doc.next_sibling("");
      m_globalConfig = xml_globalConfig(m_doc.first_node("RunCollection")->first_node("globalConfig"));

      m_files = xml_util::getVectorOfT<xml_file>(m_doc.first_node("RunCollection")->first_node("fileList"));

    }

    const xml_globalConfig& globalConfig() const {
      return m_globalConfig;
    }
    const std::vector<xml_file>& fileList() const {
      return m_files;

    }


  private:
    rapidxml::xml_document<> m_doc;
    rapidxml::file<> m_file;
    xml_globalConfig m_globalConfig;
    std::vector<xml_file> m_files;
  };
}
#endif // xml_fileList_h__