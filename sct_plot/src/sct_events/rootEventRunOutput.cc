#include "sct_events/rootEventRunOutput.hh"
#include "TTree.h"
#include <string>
#include <iostream>


std::vector<sct_type::AxesName_t > get_axis_list(){
  return{
    getIDString(),
    axis2String(x_axis_def),
    axis2String(y_axis_def),
    axis2String(Occupancy_axis_def),
    axis2String(Occupancy_error_axis_def),
    axis2String(NumOfEvents_axis_def)
  };
    
}
namespace sct_corr{
    
    rootEventRunOutput::rootEventRunOutput(TTree* tree) : rootEventBase(tree),
    m_x(getData(x_axis_def)),
    m_y(getData(y_axis_def)),
    m_Occupancy(getData(Occupancy_axis_def)),
    m_NumOfEvents(getData(NumOfEvents_axis_def)),
    m_id(getData(getIDString()))
    {
        m_totalNumOfEvents = std::make_shared< double >(0);
        m_total_efficiency = std::make_shared< double >(0);
        m_error_efficiency = std::make_shared< double >(0);
        m_residual = std::make_shared< double >(0);
        m_offset = std::make_shared< double >(0);
        m_rotation= std::make_shared< double >(0);
        m_Threshold = std::make_shared< double >(0);
        m_RunNumber = std::make_shared< double >(0);
        m_HV = std::make_shared< double >(0);
        
        
        tree->SetBranchAddress("totalNumOfEvents", m_totalNumOfEvents.get());
        tree->SetBranchAddress("total_efficiency", m_total_efficiency.get());
        tree->SetBranchAddress("error_efficiency", m_error_efficiency.get());
        tree->SetBranchAddress("residual", m_residual.get());
        tree->SetBranchAddress("offset", m_offset.get());
        tree->SetBranchAddress("rotation", m_rotation.get());
        tree->SetBranchAddress("Threshold", m_Threshold.get());
        tree->SetBranchAddress("RunNumber", m_RunNumber.get());
        tree->SetBranchAddress("HV", m_HV.get());
    }
    
    rootEventRunOutput::rootEventRunOutput(const sct_type::collectionName_t& collectionName) :rootEventBase(collectionName, get_axis_list()),
    m_x(getData(x_axis_def)),
    m_y(getData(y_axis_def)),
    m_Occupancy(getData(Occupancy_axis_def)),
    m_NumOfEvents(getData(NumOfEvents_axis_def)),
    m_id(getData(getIDString()))
    {
        m_totalNumOfEvents = std::make_shared< double >(0);
        m_total_efficiency = std::make_shared< double >(0);
        m_error_efficiency = std::make_shared< double >(0);
        m_residual = std::make_shared< double >(0);
        m_offset = std::make_shared< double >(0);
        m_rotation = std::make_shared< double >(0);
        m_Threshold = std::make_shared< double >(0);
        m_RunNumber = std::make_shared< double >(0);
        m_HV = std::make_shared< double >(0);
    }
    
    std::shared_ptr<plane> rootEventRunOutput::createPlane(const sct_type::ID_t& ID)
    {
        return nullptr;
    }
    
    void rootEventRunOutput::Save2Tree(TTree* outputTree)
    {
        rootEventBase::Save2Tree(outputTree);
        outputTree->Branch("totalNumOfEvents", m_totalNumOfEvents.get());
        outputTree->Branch("total_efficiency", m_total_efficiency.get());
        outputTree->Branch("error_efficiency", m_error_efficiency.get());
        outputTree->Branch("residual", m_residual.get());
        outputTree->Branch("rotation", m_rotation.get());
        outputTree->Branch("offset", m_offset.get());
        outputTree->Branch("Threshold", m_Threshold.get());
        outputTree->Branch("RunNumber", m_RunNumber.get());
        outputTree->Branch("HV", m_HV.get());
    }
    
    void rootEventRunOutput::set_TotalNumOfEvents(double numOfEvents)
    {
 
        *m_totalNumOfEvents = numOfEvents;
    }
    
    void rootEventRunOutput::set_Total_efficiency(double effi)
    {
        *m_total_efficiency = effi;
    }
    
    void rootEventRunOutput::set_Error_efficiency(double error_effi)
    {
        *m_error_efficiency = error_effi;
    }
    
    void rootEventRunOutput::set_Threshold(double thr)
    {
        *m_Threshold = thr;
    }
    
    void rootEventRunOutput::set_rotation(double rotation_) {
      *m_rotation = rotation_;
    }

    void rootEventRunOutput::set_HV(double HV__)
    {
        *m_HV = HV__;
    }
    
    void rootEventRunOutput::set_residual(double res)
    {
        *m_residual = res;
    }
    
    void rootEventRunOutput::set_offset(double offset)
    {
        *m_offset = offset;
    }
    
    void rootEventRunOutput::set_RunNumber(double RunNum)
    {
        *m_RunNumber = RunNum;
    }
    
    void rootEventRunOutput::reset()
    {
        rootEventBase::reset();
        *m_totalNumOfEvents = 0;
        *m_total_efficiency = 0;
        *m_error_efficiency = 0;
        *m_residual = 0;
        *m_offset = 0;
        *m_Threshold = 0;
        *m_RunNumber = 0;
        *m_HV = 0;
        
    }
    
}