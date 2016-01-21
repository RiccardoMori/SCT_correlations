#ifndef s_cuts_h__
#define s_cuts_h__

#include "internal/platform.hh"
#include <memory>
#include "Rtypes.h"

class DllExport S_Cut
{
public:

  virtual ~S_Cut(){}

#ifndef __CINT__
  virtual std::shared_ptr<S_Cut> copy()const { 
    return  std::shared_ptr<S_Cut>(new S_Cut(*this));
  }
#endif // __CINT__

  virtual bool isOutOfRange(Double_t BinContent, Double_t x, Double_t y) const { return false; }
  virtual bool isOutOfRange(Double_t BinContent, Double_t x)  const { return false; }
  virtual bool isOutOfRange(Double_t BinContent)  const { return false; }
 
};
class DllExport S_NoCut:public S_Cut
{
public:

  virtual ~S_NoCut(){}

#ifndef __CINT__
  virtual std::shared_ptr<S_Cut> copy()const {
    return  std::shared_ptr<S_Cut>(new S_NoCut(*this));
  }
#endif // __CINT__

  
};
class DllExport S_Cut_min_max : public S_Cut
{
public:
  virtual ~S_Cut_min_max(){}
  S_Cut_min_max(Double_t min_, Double_t max_);
  S_Cut_min_max(Double_t min_);

protected:
#ifndef __CINT__

  Double_t m_min, m_max;
  Bool_t m_cut_min, m_cut_max;
  bool isOutOfRange_intern(Double_t val)  const;
#endif // !__CINT__

  
};

class DllExport  S_XCut :public S_Cut_min_max
{
public:
  S_XCut(Double_t min_, Double_t max_);
  S_XCut(Double_t min_);
#ifndef __CINT__
  virtual std::shared_ptr<S_Cut> copy()const {
    return  std::shared_ptr<S_Cut>(new S_XCut(*this));
  }
#endif // __CINT__
  virtual bool isOutOfRange(Double_t BinContent, Double_t x, Double_t y) const;
  virtual bool isOutOfRange(Double_t BinContent, Double_t x) const;
  virtual bool isOutOfRange(Double_t BinContent)  const;

 
};
class DllExport S_YCut :public S_Cut_min_max
{
public:
  S_YCut(Double_t min_, Double_t max_);
  S_YCut(Double_t min_);
#ifndef __CINT__
  virtual std::shared_ptr<S_Cut> copy()const {
    return  std::shared_ptr<S_Cut>(new S_YCut(*this));
  }
#endif // __CINT__
  virtual bool isOutOfRange(Double_t BinContent, Double_t x, Double_t y) const;
  virtual bool isOutOfRange(Double_t BinContent, Double_t x) const;
  virtual bool isOutOfRange(Double_t BinContent) const;
 
};
class  DllExport S_Cut_BinContent :public S_Cut_min_max
{
public:
  S_Cut_BinContent(Double_t min_, Double_t max_);
  S_Cut_BinContent(Double_t min_);
#ifndef __CINT__
  virtual std::shared_ptr<S_Cut> copy()const override{
    return  std::shared_ptr<S_Cut>(new S_Cut_BinContent(*this));
  }
#endif // __CINT__
  virtual bool isOutOfRange(Double_t BinContent, Double_t x, Double_t y) const override;
  virtual bool isOutOfRange(Double_t BinContent, Double_t x) const override;
  virtual bool isOutOfRange(Double_t BinContent) const override;

 
};

class  DllExport S_CutCoollection :public S_Cut
{
public:

#ifndef __CINT__
  virtual std::shared_ptr<S_Cut> copy()const {
    return  std::shared_ptr<S_Cut>(new S_CutCoollection(*this));
  }
#endif // __CINT__
  S_CutCoollection();
  virtual bool isOutOfRange(Double_t BinContent, Double_t x, Double_t y) const;
  virtual bool isOutOfRange(Double_t BinContent, Double_t x) const;
  virtual bool isOutOfRange(Double_t BinContent) const;
  void add_Cut(const S_Cut& cut_);
  void add_Cut(const S_CutCoollection& cut_);
#ifndef __CINT__
  std::vector<std::shared_ptr<S_Cut> > m_cuts;
#endif // !__CINT__

};
DllExport S_CutCoollection operator+(const S_Cut& cut_a, const S_Cut& cut_b);



#endif // s_cuts_h__
