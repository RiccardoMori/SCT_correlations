#ifndef setup_description_h__
#define setup_description_h__
#include <string>
#include <vector>
#include "rapidxml.hpp"
#include "rapidxml_utils.hpp"
#include "rapidxml.hpp"
#include "xml_helpers/xml_util.hh"
#include "TMath.h"
#include "internal/exceptions.hh"
#define  sct_corr_XML_ERROR_DEFAULT_VALUE  double(-666666.666666)

namespace sct_corr {
using xml_n = rapidxml::xml_node < char >;
struct Xladder {
public:
  Xladder(xml_n* n) :

    ID(xml_util::getAtribute(n, "ID", sct_corr_XML_ERROR_DEFAULT_VALUE)),
    positionX(xml_util::getAtribute(n, "positionX", sct_corr_XML_ERROR_DEFAULT_VALUE)),
    positionY(xml_util::getAtribute(n, "positionY", sct_corr_XML_ERROR_DEFAULT_VALUE)),
    positionZ(xml_util::getAtribute(n, "positionZ", sct_corr_XML_ERROR_DEFAULT_VALUE)),
    radLength(xml_util::getAtribute(n, "radLength", sct_corr_XML_ERROR_DEFAULT_VALUE)),
    rotationXY(xml_util::getAtribute(n, "rotationXY", sct_corr_XML_ERROR_DEFAULT_VALUE)),
    rotationZX(xml_util::getAtribute(n, "rotationZX", sct_corr_XML_ERROR_DEFAULT_VALUE)),
    rotationZY(xml_util::getAtribute(n, "rotationZY", sct_corr_XML_ERROR_DEFAULT_VALUE)),


    sizeX(xml_util::getAtribute(n, "sizeX", sct_corr_XML_ERROR_DEFAULT_VALUE)),
    sizeY(xml_util::getAtribute(n, "sizeY", sct_corr_XML_ERROR_DEFAULT_VALUE)),
    thickness(xml_util::getAtribute(n, "thickness", sct_corr_XML_ERROR_DEFAULT_VALUE)) {
  }
  double ID = sct_corr_XML_ERROR_DEFAULT_VALUE,
    positionX = sct_corr_XML_ERROR_DEFAULT_VALUE,
    positionY = sct_corr_XML_ERROR_DEFAULT_VALUE,
    positionZ = sct_corr_XML_ERROR_DEFAULT_VALUE,
    radLength = sct_corr_XML_ERROR_DEFAULT_VALUE,
    rotationXY = sct_corr_XML_ERROR_DEFAULT_VALUE,
    rotationZX = sct_corr_XML_ERROR_DEFAULT_VALUE,
    rotationZY = sct_corr_XML_ERROR_DEFAULT_VALUE,
    sizeX = sct_corr_XML_ERROR_DEFAULT_VALUE,
    sizeY = sct_corr_XML_ERROR_DEFAULT_VALUE,
    thickness = sct_corr_XML_ERROR_DEFAULT_VALUE;


};
struct Xsensitive {
public:
  Xsensitive(xml_n* n) :

    ID(xml_util::getAtribute(n, "ID", sct_corr_XML_ERROR_DEFAULT_VALUE)),
    npixelX(xml_util::getAtribute(n, "npixelX", sct_corr_XML_ERROR_DEFAULT_VALUE)),
    npixelY(xml_util::getAtribute(n, "npixelY", sct_corr_XML_ERROR_DEFAULT_VALUE)),


    pitchX(xml_util::getAtribute(n, "pitchX", sct_corr_XML_ERROR_DEFAULT_VALUE)),
    pitchY(xml_util::getAtribute(n, "pitchY", sct_corr_XML_ERROR_DEFAULT_VALUE)),
    positionX(xml_util::getAtribute(n, "positionX", sct_corr_XML_ERROR_DEFAULT_VALUE)),
    positionY(xml_util::getAtribute(n, "positionY", sct_corr_XML_ERROR_DEFAULT_VALUE)),
    positionZ(xml_util::getAtribute(n, "positionZ", sct_corr_XML_ERROR_DEFAULT_VALUE)),
    radLength(xml_util::getAtribute(n, "radLength", sct_corr_XML_ERROR_DEFAULT_VALUE)),

    resolution(xml_util::getAtribute(n, "resolution", sct_corr_XML_ERROR_DEFAULT_VALUE)),
    rotation1(xml_util::getAtribute(n, "rotation1", sct_corr_XML_ERROR_DEFAULT_VALUE)),
    rotation2(xml_util::getAtribute(n, "rotation2", sct_corr_XML_ERROR_DEFAULT_VALUE)),
    rotation3(xml_util::getAtribute(n, "rotation3", sct_corr_XML_ERROR_DEFAULT_VALUE)),


    rotation4(xml_util::getAtribute(n, "rotation4", sct_corr_XML_ERROR_DEFAULT_VALUE)),
    sizeX(xml_util::getAtribute(n, "sizeX", sct_corr_XML_ERROR_DEFAULT_VALUE)),
    sizeY(xml_util::getAtribute(n, "sizeY", sct_corr_XML_ERROR_DEFAULT_VALUE)),
    thickness(xml_util::getAtribute(n, "thickness", sct_corr_XML_ERROR_DEFAULT_VALUE))

  {
  }
  double ID = sct_corr_XML_ERROR_DEFAULT_VALUE,
    npixelX = sct_corr_XML_ERROR_DEFAULT_VALUE,
    npixelY = sct_corr_XML_ERROR_DEFAULT_VALUE,
    pitchX = sct_corr_XML_ERROR_DEFAULT_VALUE,
    pitchY = sct_corr_XML_ERROR_DEFAULT_VALUE,
    positionX = sct_corr_XML_ERROR_DEFAULT_VALUE,
    positionY = sct_corr_XML_ERROR_DEFAULT_VALUE,
    positionZ = sct_corr_XML_ERROR_DEFAULT_VALUE,
    radLength = sct_corr_XML_ERROR_DEFAULT_VALUE,
    resolution = sct_corr_XML_ERROR_DEFAULT_VALUE,
    rotation1 = sct_corr_XML_ERROR_DEFAULT_VALUE,
    rotation2 = sct_corr_XML_ERROR_DEFAULT_VALUE,
    rotation3 = sct_corr_XML_ERROR_DEFAULT_VALUE,
    rotation4 = sct_corr_XML_ERROR_DEFAULT_VALUE,
    sizeX = sct_corr_XML_ERROR_DEFAULT_VALUE,
    sizeY = sct_corr_XML_ERROR_DEFAULT_VALUE,
    thickness = sct_corr_XML_ERROR_DEFAULT_VALUE;

};
struct Xlayer {
public:
  Xlayer(xml_n* n) :ladder(n->first_node("ladder")), sensitive(n->first_node("sensitive")) {}
  Xladder ladder;
  Xsensitive  sensitive;
  static const char* NodeName() {
    return "layer";
  }
};
struct Xdetector {
public:
  Xdetector(xml_n* n) {
    geartype = xml_util::getAtribute(n, "geartype", "notSet");
    name = xml_util::getAtribute(n, "name", "notSet");

    siplanesID = xml_util::getAtribute(n->first_node("siplanesID"), "ID", siplanesID);
    siplanesType = xml_util::getAtribute(n->first_node("siplanesType"), "type", "notSet");

    siplanesNumber = xml_util::getAtribute(n->first_node("siplanesNumber"), "number", siplanesNumber);

    layer = xml_util::getVectorOfT<Xlayer>(n->first_node("layers"));
  }
  std::string geartype, name;
  double siplanesID = 0;
  std::string siplanesType;
  double siplanesNumber = 0;
  Xlayer* layer_by_ID(double ID) {
    for (auto& e : layer) {
      if (e.ladder.ID == ID) {
        return &e;
      }
    }
    SCT_THROW("Layer not Found. searching for Layer " + std::to_string(ID));
    return nullptr;
  }
  const Xlayer* layer_by_ID(double ID) const {
    for (auto& e : layer) {
      if (e.ladder.ID == ID) {
        return &e;
      }
    }
    SCT_THROW("Layer not Found. searching for Layer " + std::to_string(ID));
    return nullptr;
  }
  std::vector<Xlayer> layer;
};
struct XBField {
public:
  XBField(xml_n* n) :type(xml_util::getAtribute(n, "type", "notSet")), x(xml_util::getAtribute(n, "x", (double)0.0)), y(xml_util::getAtribute(n, "y", (double)0.0)), z(xml_util::getAtribute(n, "z", (double)0.0)) {}
  std::string type;
  double x, y, z;
};

struct Xgear {

public:
  Xgear(xml_n* n) :detectorName(xml_util::getAtribute(n->first_node("global"), "detectorName", "notSet")), bfield(n->first_node("BField")), detector(n->first_node("detectors")->first_node("detector")) {}


  std::string detectorName;
  XBField bfield;
  Xdetector detector;
};

inline Xgear load_gear(const char* name) {
  rapidxml::file<> m_file(name);
  rapidxml::xml_document<> m_doc;
  m_doc.parse<0>(m_file.data());

  return sct_corr::Xgear(m_doc.first_node("gear"));
}
}

#endif // setup_description_h__
