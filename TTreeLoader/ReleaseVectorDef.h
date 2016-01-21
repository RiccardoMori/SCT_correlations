#ifndef ReleaseVectorDef_h__
#define ReleaseVectorDef_h__
#include <memory>

#define DllImport   __declspec( dllimport )
#define DllExport   __declspec( dllexport )
class TTree;

class DllExport ReleaseVectorDef{
public:
  ReleaseVectorDef(const char* name);
  bool load_from_TTree(TTree* tree);
  bool push2TTree(TTree* tree);
  double at(size_t i) const;
  size_t size() const;
  void clear();
  void push_back(double val);

private:
  class vectorImpl;
  std::shared_ptr<vectorImpl> m_vecImpl;
};
#endif // ReleaseVectorDef_h__
