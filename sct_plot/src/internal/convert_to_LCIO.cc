#include "sct_plots.h"
#ifdef USE_LCIO
#include "internal/plotsBase.hh"
#include <iostream>



#include "lcio.h"
#include "IMPL/LCEventImpl.h"
#include "IMPL/TrackerRawDataImpl.h"
#include "IMPL/LCRunHeaderImpl.h"
#include "IMPL/LCCollectionVec.h"
#include "EUTELESCOPE.h"
#include "EUTelRunHeaderImpl.h"


#include "EUTelSetupDescription.h"
#include "EUTelTrackerDataInterfacerImpl.h"
#  include "IMPL/LCEventImpl.h"
#  include "IMPL/TrackerRawDataImpl.h"
#  include "IMPL/TrackerDataImpl.h"
#  include "IMPL/LCCollectionVec.h"
#  include "IMPL/LCGenericObjectImpl.h"
#  include "UTIL/CellIDEncoder.h"
#  include "lcio.h"
#include "IMPL/TrackerHitImpl.h"
#include "EVENT/LCIO.h"
#include "UTIL/CellIDDecoder.h"
#include "CellIDReencoder.h"
#include "s_plot_prob.h"
#include "internal/plane_hit.hh"
#include "internal/plane.hh"
using namespace IMPL;

using lcio::TrackerDataImpl;
using eutelescope::EUTELESCOPE;


using eutelescope::EUTelTrackerDataInterfacerImpl;
using eutelescope::EUTelGenericSparsePixel;
namespace sct_corr{

  class convert_to_LCIO :public plot{
  public:
    convert_to_LCIO(const char *file_name, unsigned runNum, const s_plot_prob& = "");
    virtual const char* getType() const;
    virtual bool MakeReadyForData(sct_event_buffer* outputBuffer);
    virtual void pushAxis(const axis_ref* axis);
    virtual void pushPlane(S_plane* plane_);
    virtual returnProcessEvent ProcessCurrentEvent() override final;
    virtual void fill() override final;
    virtual Long64_t Draw(const char* options, const char* cuts = "", const char* axis = "y:x");
    virtual Long64_t Draw(const S_DrawOption&);

    std::unique_ptr<lcio::LCRunHeader>  GetLCRunHeader();
    void GetLCIORunHeader(lcio::LCRunHeader & header, S_plane* plane_) const;

    std::unique_ptr<lcio::LCEvent > ConvertToLCIO();
    void ConvertLCIOSubEvent(lcio::LCEvent & dest, S_plane* source);
    void ConvertPlaneToLCIOGenericPixel(S_plane* source, lcio::TrackerDataImpl& zsFrame);
    virtual sct_type::collectionName_t getOutputName()  const override;

    std::string m_filename;
    std::vector<S_plane*> m_planes;
    lcio::LCWriter *m_lcwriter = nullptr; /// The lcio writer
    bool m_fileopened = false; /// We have to keep track whether a file is open ourselves
    unsigned m_runNum = 0, m_event_nr = 0;

    bool Collection_createIfNotExist(lcio::LCCollectionVec** zsDataCollection, const lcio::LCEvent & lcioEvent, const char * name);
  };

  convert_to_LCIO::convert_to_LCIO(const char *file_name, unsigned runNum, const s_plot_prob& plot_prob) :
    plot(s_plot_prob(plot_prob).SaveToDisk()),
    m_filename(file_name),
    m_runNum(runNum)
  {
    m_lcwriter = lcio::LCFactory::getInstance()->createLCWriter();
  }

  const char* convert_to_LCIO::getType() const
  {
    return "conveter_to_LCIO";
  }

  bool convert_to_LCIO::MakeReadyForData(sct_event_buffer* outputBuffer)
  {
    m_event_nr = 0;
    if (m_fileopened) {
      m_lcwriter->close();
      m_fileopened = false;
    }

    // open a new file
    try {
      m_lcwriter->open(m_filename.c_str(), lcio::LCIO::WRITE_NEW);
      m_fileopened = true;
    }
    catch (const lcio::IOException & e) {
      std::cout << e.what() << std::endl;
      ///FIXME Error message to run control and logger
      return m_fileopened;
    }


    auto lcRunHeader = GetLCRunHeader();
    m_lcwriter->writeRunHeader(lcRunHeader.get());

    return m_fileopened;
  }

  void convert_to_LCIO::pushAxis(const axis_ref* axis)
  {
    std::cout << "[convert_to_LCIO] axis are not supported \n";
  }

  void convert_to_LCIO::pushPlane(S_plane* plane_)
  {
    m_planes.push_back(plane_);
  }

  sct_corr::returnProcessEvent convert_to_LCIO::ProcessCurrentEvent() {
    return FILL_OK;
  }

  void convert_to_LCIO::fill()
{
    //  std::cout << "EUDAQ_DEBUG: FileWriterLCIO::WriteEvent() processing event "
    //    <<  devent.GetRunNumber () <<"." << devent.GetEventNumber () << std::endl;

    auto lcevent = ConvertToLCIO();

    // only write non-empty events
    //std::cout << lcevent->getDetectorName() << std::endl;
    if (!lcevent->getCollectionNames()->empty()) {
      //  std::cout << " FileWriterLCIO::WriteEvent() : doing the actual writing : " <<std::flush;
      m_lcwriter->writeEvent(lcevent.get());
      //  std::cout << " done" <<std::endl;
    }
   
  }


  Long64_t convert_to_LCIO::Draw(const S_DrawOption&)
  {
    std::cout << "[convert_to_LCIO] drawing not supported \n";
    return -1;
  }

  Long64_t convert_to_LCIO::Draw(const char* options, const char* cuts /*= ""*/, const char* axis /*= "y:x"*/)
  {
    return Draw(S_DrawOption(options, cuts, axis));
  }

  std::unique_ptr<lcio::LCRunHeader> convert_to_LCIO::GetLCRunHeader()
  {
    auto ret = std::unique_ptr<lcio::LCRunHeader>(new IMPL::LCRunHeaderImpl());
    auto lcHeader = dynamic_cast<IMPL::LCRunHeaderImpl*>(ret.get());
    lcHeader->setRunNumber(m_runNum);
    lcHeader->setDetectorName("EUTelescope");
    eutelescope::EUTelRunHeaderImpl runHeader(lcHeader);
    runHeader.setDateTime();
    runHeader.setDataType(EUTELESCOPE::DAQDATA);
    runHeader.setDAQSWName(EUTELESCOPE::EUDAQ);


    runHeader.setGeoID(0);

    for (auto& e : m_planes) {

      GetLCIORunHeader(*lcHeader, e);
    }
    return ret;
  }

  void convert_to_LCIO::GetLCIORunHeader(lcio::LCRunHeader & header, S_plane* plane_) const
  {

  }

  std::unique_ptr<lcio::LCEvent >  convert_to_LCIO::ConvertToLCIO()
  {
    auto ret = std::unique_ptr<lcio::LCEvent >(new lcio::LCEventImpl());
    auto lcEv = dynamic_cast<lcio::LCEventImpl*>(ret.get());
    lcEv->setEventNumber(m_event_nr++);
    lcEv->setRunNumber(m_runNum);
    lcEv->setTimeStamp(0);

    for (auto& e : m_planes) {
      ConvertLCIOSubEvent(*lcEv, e);
    }

    return ret;
  }

  void convert_to_LCIO::ConvertLCIOSubEvent(lcio::LCEvent & lcEv, S_plane* source)
  {

    // If we get here it must be a data event


    lcEv.parameters().setValue(eutelescope::EUTELESCOPE::EVENTTYPE, eutelescope::kDE);


    LCCollectionVec* outputCollection = nullptr;
    try
    {
      outputCollection = static_cast<LCCollectionVec*> (lcEv.getCollection(necessary_CONVERSION(source->getName()))); //external library 
    }
    catch (...)
    {
      outputCollection = new LCCollectionVec(LCIO::TRACKERHIT);
    }



    std::string encoding = EUTELESCOPE::HITENCODING;
    lcio::CellIDDecoder<TrackerHitImpl> hitDecoder(encoding);
    lcio::UTIL::CellIDReencoder<TrackerHitImpl> cellReencoder(encoding, outputCollection);


    auto id = source->getID();
    double outputPos[3];
    while (source->next())
    {

      TrackerHitImpl* outputHit = new IMPL::TrackerHitImpl();

      //Call the local2masterHit/master2localHit function defined int EUTelGeometryTelescopeDescription

      auto hit = source->getPlane()->getHit();
      outputPos[0] = hit->x;
      outputPos[1] = hit->y;
      outputPos[2] = 1;



      //Fill the new outputHit with information
      outputHit->setPosition(outputPos);
      float cov[TRKHITNCOVMATRIX] = { 0., 0., 0., 0., 0., 0. };
      cov[0] = 1;
      cov[2] = 1;
      outputHit->setCovMatrix(cov);
      outputHit->setType(0);
      outputHit->setTime(0);
      outputHit->setCellID0(static_cast<int>(necessary_CONVERSION(id))); //external library 
      outputHit->setCellID1(static_cast<int>(necessary_CONVERSION(id))); //external library 



      cellReencoder.readValues(outputHit);


      cellReencoder.setCellID(outputHit);

      outputCollection->push_back(outputHit);
    }

    try
    {
      lcEv.addCollection(outputCollection, necessary_CONVERSION(source->getName()));//external library 
    }
    catch (...)
    {
      std::cout << "Problem with pushing collection onto event" << std::endl;
    }
  }



  sct_type::collectionName_t convert_to_LCIO::getOutputName() const
  {
    return sct_type::collectionName_t(m_filename.c_str());
  }

  bool convert_to_LCIO::Collection_createIfNotExist(lcio::LCCollectionVec** zsDataCollection, const lcio::LCEvent & lcioEvent, const char * name)
  {
    bool zsDataCollectionExists = false;
    try
    {
      *zsDataCollection = static_cast<lcio::LCCollectionVec*> (lcioEvent.getCollection(name));
      zsDataCollectionExists = true;
    }
    catch (lcio::DataNotAvailableException& )
    {
      *zsDataCollection = new lcio::LCCollectionVec(lcio::LCIO::TRACKERDATA);
    }

    return zsDataCollectionExists;
  }
}

sct_corr::S_plot sct_plot::save2LCIO(const char* filename, unsigned runnum, const s_plot_prob& plot_prob)
{
  return sct_corr::S_plot(new sct_corr::convert_to_LCIO(filename, runnum, plot_prob));

}
#else
#include "s_plot_prob.h"
sct_corr::S_plot sct_plot::save2LCIO(const char* filename, unsigned runnum, const s_plot_prob& plot_prob)
{
  return sct_corr::S_plot(nullptr);
}

#endif